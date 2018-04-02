/**
 *	@file MultShootEngine.cpp
 *	@brief Engine object that applies differential corrections to arcsets
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of Astrohelion
 *
 *  Astrohelion is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Astrohelion is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Astrohelion.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <cmath>
#include <vector>

#include "MultShootEngine.hpp"

#include "AsciiOutput.hpp"
#include "Calculations.hpp"
#include "EigenDefs.hpp"
#include "Exceptions.hpp" 
#include "MultShootData.hpp"
#include "Node.hpp"
#include "Arcset_cr3bp.hpp"
#include "SimEngine.hpp"
#include "Utilities.hpp"

using Eigen::MatrixXd;
using Eigen::JacobiSVD;

namespace astrohelion{
//-----------------------------------------------------------------------------
//      *structors
//-----------------------------------------------------------------------------

/**
 *  @brief Default constructor
 */
MultShootEngine::MultShootEngine(){}

/**
 *	@brief Copy constructor - create this engine by copying the input engine
 *	@param e input correction engine
 */
MultShootEngine::MultShootEngine(const MultShootEngine &e){
	copyMe(e);
}//=======================================================

/**
 *	@brief Destructor
 */
MultShootEngine::~MultShootEngine(){}

/**
 *	@brief Copy all engine variables
 *	@param e an engine reference
 */
void MultShootEngine::copyMe(const MultShootEngine &e){
	Engine::copyBaseEngine(e);
	verbosity = e.verbosity;
	maxIts = e.maxIts;
	tolF = e.tolF;
	tolX = e.tolX;
	tolA = e.tolA;
	bInitLUSolver = e.bInitLUSolver;
	bInitQRSolver = e.bInitQRSolver;
	bFindEvent = e.bFindEvent;
	bIgnoreCrash = e.bIgnoreCrash;
	bIgnoreDiverge = e.bIgnoreDiverge;
	bFullFinalProp = e.bFullFinalProp;
	bLUFailed = e.bLUFailed;
	bSaveEachIt = e.bSaveEachIt;
	tofTp = e.tofTp;
	bLineSearchAttenFactor = e.bLineSearchAttenFactor;
	ls_alpha = e.ls_alpha;
	ls_maxStepSize = e.ls_maxStepSize;
}//====================================================

//-----------------------------------------------------------------------------
//      Operator Functions
//-----------------------------------------------------------------------------

/**
 *	@brief Copy operator; make a copy of the input correction engine.
 *
 *	@param e
 *	@return this correction engine
 */
MultShootEngine& MultShootEngine::operator =(const MultShootEngine &e){
	copyMe(e);
	return *this;
}//====================================

//-----------------------------------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------------------------------

/**
 *  @brief Retrieve whether or not the engine is locating an event crossing
 *	@return whether or not the algorithm will optimize the process to find an 
 *	event
 */
bool MultShootEngine::isFindingEvent() const { return bFindEvent; }

/**
 *  @brief Retreive whether or not the engine will use a full, variable-step
 *  propagation for the final propagation.
 *  @details By default, this setting is TRUE. For lower computation time, 
 *  set to false via setFullFinalProp().
 *  
 *  @return whether or not the engine will use a full, variable-step
 *  propagation for the final propagation
 */
bool MultShootEngine::doesFullFinalProp() const { return bFullFinalProp; }

/**
 *  @brief Retreive whether or not the engine conducts a line search to choose
 *  the Newton step size.
 *  @return whether or not the engine conducts a line search to choose
 *  the Newton step size.
 */
bool MultShootEngine::doesLineSearch() const { return bLineSearchAttenFactor; }

/**
 *  @brief Retrieve the maximum permitted error magnitude
 *  @details If the error rises above this value, the corrections processes
 *  is deemed diverged and the process will quit.
 *  @return the maximum permitted error magnitude
 */
double MultShootEngine::getMaxErr() const{ return maxErr; }

/**
 *  @brief Retrieve the maximum number of iterations to attempt
 *	@return the maximum number of iterations to attempt before giving up
 */
int MultShootEngine::getMaxIts() const { return maxIts; }

/**
 *  @brief Retrieve the enumerated type describing how time-of-flight values
 *  are encoded (if at all) in the free variable vector
 *  
 *  @return the time-of-flight type
 */
MSTOF_tp MultShootEngine::getTOFType() const{ return tofTp; }

/**
 *  @brief Retrieve the minimum error tolerance
 *	@return the minimum error tolerance (non-dimensional units); errors
 *	less than this value are considered negligible
 */
double MultShootEngine::getTol() const { return tolF; }

/**
 *  @brief Set whether or not the engine conducts a line search to choose 
 *  the Newton step size.
 *  @details Although the Newton step direction is guaranteed to point toward a 
 *  decreasing constraint vector, the full step may be too large. By leveraging
 *  a line search, the step size is chosen such that the magnitude of the 
 *  constraint vector decreases. WARNING: This adds many evaluations of 
 *  expensive functions and will greatly decrease the speed. The line search is 
 *  only recommended for particularly stubborn correction processes.
 * 
 *  @param b 
 */
void MultShootEngine::setDoLineSearch(bool b){ bLineSearchAttenFactor = b; }

/**
 *  @brief Set whether or not the engine will use a full, variable-step
 *  propagation for the final propagation
 *  @details By default, this setting is TRUE. For lower computation time, 
 *  set to false. When set to false, only two points per segment are stored
 *  (the initial and final states), i.e., no intermediate points between nodes
 *  are computed
 * 
 *  @param b whether or not the engine will use a full, variable-step
 *  propagation for the final propagation
 */
void MultShootEngine::setFullFinalProp(bool b){ bFullFinalProp = b; }

/**
 * @brief Tell the corrector to ignore crash events (or to not to).
 * @details By default, the corrector does monitor crashes and will 
 * run into issues if the trajectory being corrected passes through a 
 * primary.
 * 
 * @param b whether or not to ignore crashes (default is false)
 */
void MultShootEngine::setIgnoreCrash(bool b){ bIgnoreCrash = b; }

/**
 *  @brief Tell the corrector to ignore divergence and return the partially
 *  corrected iteration data instead of throwing an exception when divergence
 *  occurs.
 * 
 *  @param b Whether or not to ignore divergance
 */
void MultShootEngine::setIgnoreDiverge(bool b){ bIgnoreDiverge = b;}

/**
 *  @brief Set the maximum permitted error magnitude
 *  @details If the error rises above this value, the corrections processes
 *  is deemed diverged and the process will quit.
 *  
 *  @param e the maximum permitted error magnitude
 */
void MultShootEngine::setMaxErr(double e){ maxErr = e; }

/**
 *	@brief Set maximum iterations
 *	@param i the maximum number of iterations to attempt before giving up
 */
void MultShootEngine::setMaxIts(int i){ maxIts = i; }

/**
 * @brief Set the flag to tell the engine whether or not to save each
 * updated arcset to file
 * @details Files are saved with the format "multShoot_it####.mat"
 * 
 * @param bSave whether or not to save each updated arcset to file
 */
void MultShootEngine::setSaveEachIt(bool bSave){ bSaveEachIt = bSave; }

/**
 *  @brief Set the way times-of-flight are encoded (if at all) in the
 *  free variable vector
 * 
 *  @param tp Describes how times-of-flight are encoded in the free
 *  variable vector
 */
void MultShootEngine::setTOFType(MSTOF_tp tp){ tofTp = tp; }

/**
 *	@brief Set the error tolerance
 *	@param d errors below this value will be considered negligible
 */
void MultShootEngine::setTol(double d){
	tolF = d;

	if(tolF > 1)
		astrohelion::printWarn("MultShootEngine::setTol: "
			"tolerance is greater than 1... just FYI\n");
}//====================================================

/**
 *	@brief Set the findEven flag
 *	@param b whether or not the algorithm will be looking for an event
 */
void MultShootEngine::setFindEvent(bool b){ bFindEvent = b; }

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Correct a generic arcset using multiple shooting
 *	@details This algorithm employs multiple shooting to correct a set of nodes
 *	subject to a set of constraints. The nodes and constraints are all stored in the 
 *	input arcset object. 
 *	
 *	Numerical integration is used to generate state transition matrices between
 *	nodes. In order to optimize performance and accuracy, only two steps are 
 *	computed by the simulation engine, and the step size is fixed to force the 
 *	usage of the Adams-Bashforth Adams-Moulton method.
 *	
 *	@param pArcIn pointer to the arcset that needs to be corrected
 *	@param pArcOut pointer to the arcset object that will contain the results of
 *	the shooting process
 *	@return the iteration data object for this corrections process
 *	@throws DivergeException if the corrections process does not converge
 *	@throws Exception
 *	* if the input and output arcsets contain different system data objects
 *	* if the dynamic model associated with the input
 *	arcset does not support one or more of the arcset's constraints
 *	* if the input arcset contains more than one delta-v constraint
 *	* if the input arcset contains more than one TOF constraint
 */
void MultShootEngine::multShoot(const Arcset *pArcIn, Arcset *pArcOut,
	MultShootData *pData){

	if(pArcOut != nullptr && *(pArcIn->getSysData()) != *(pArcOut->getSysData()))
		throw Exception("MultShootEngine::multShoot: "
			"Input and Output arcsets must use the same system data object");

	if(!bIsClean)
		cleanEngine();

	bIsClean = false;
	bool bDataOnStack = pData == nullptr;

	if(bDataOnStack)
		pData = new MultShootData(pArcIn);
	else
		*pData = MultShootData(pArcIn);

	// Create structure to store iteration data for easy sharing
	// MultShootData it(pArcIn);
	pData->pArcOut = pArcOut;
	pData->tofTp = tofTp;
	
	if(verbosity >= Verbosity_tp::ALL_MSG){
		printf("Multiple Shooting Algorithm:\n");
		printf("  it.numNodes = %u\n", pData->numNodes);
		printf("  sysType = %s\n", pArcIn->getSysData()->getTypeStr().c_str());
		printf("  TOF Type = %s\n", MSTOF_tp_cStr(tofTp));
		printf("  Do full final prop? %s\n", bFullFinalProp ? "YES" : "NO");
		printf("  Find event? %s\n", bFindEvent ? "YES" : "NO");
		printf("  Ignore crash? %s\n", bIgnoreCrash ? "YES" : "NO");
		printf("  Ignore diverge? %s\n", bIgnoreDiverge ? "YES" : "NO");
		printf("  Do line search? %s\n", bLineSearchAttenFactor ? "YES" : "NO");
	}
	

	// Get the model associated with the arcset
	const DynamicsModel *pModel = pArcIn->getSysData()->getDynamicsModel();
	pModel->multShoot_initDesignVec(*pData);

	// Create constraints that enforce continuity between nodes; this process
	// does account for velocity discontinuities specified in the arcset
	pData->allCons.clear();
	pModel->multShoot_createContCons(*pData);

	// Add all node constraints
	for(unsigned int n = 0; n < pArcIn->getNumNodes(); n++){
		std::vector<Constraint> nodeCons = 
			pArcIn->getNodeRefByIx_const(n).getConstraints();
		pData->allCons.insert(pData->allCons.end(), nodeCons.begin(), nodeCons.end());
	}

	// Add all segment constraints
	for(unsigned int s = 0; s < pArcIn->getNumSegs(); s++){
		std::vector<Constraint> segCons = 
			pArcIn->getSegRefByIx_const(s).getConstraints();
		pData->allCons.insert(pData->allCons.end(), segCons.begin(), segCons.end());
	}

	// Add all arcset constraints
	std::vector<Constraint> arcCons = pArcIn->getArcConstraints();
	pData->allCons.insert(pData->allCons.end(), arcCons.begin(), arcCons.end());

	// Compute number of extra consraint functions to add
	pData->numSlack = 0;

	// Initialize vector to keep track of which row each constraint begins on
	// Also add slack variables to the design variable vector
	pData->conRows.assign(pData->allCons.size(), -1);	// Fill with -1
	int conRow = 0;
	bool foundDVCon = false;
	bool foundTOFCon = false;
	for(unsigned int c = 0; c < pData->allCons.size(); c++){
		int addToRows = 0;
		Constraint con = pData->allCons[c];

		if(!pModel->supportsCon(con.getType()))
			throw Exception("MultShootEngine::multShoot: "
				"The dynamic model does not support one of the constraints!");

		switch(con.getType()){
			case Constraint_tp::CONT_PV:
			case Constraint_tp::CONT_EX:
			case Constraint_tp::CONT_CTRL:
			case Constraint_tp::SEG_CONT_PV:
			case Constraint_tp::STATE:
			case Constraint_tp::CTRL:
			case Constraint_tp::MATCH_CUST:
			case Constraint_tp::ENDSEG_STATE:
				addToRows = con.countConstrainedStates();
				break;
			case Constraint_tp::MATCH_ALL:
				addToRows = pModel->getCoreStateSize();
				break;
			case Constraint_tp::SP:
				addToRows = 3;
				break;
			case Constraint_tp::SP_RANGE:
				addToRows = 1;
				pData->X.push_back(pModel->multShoot_getSlackVarVal(*pData, con));
				pData->slackAssignCon.push_back(c);
				pData->numSlack++;
				break;
			case Constraint_tp::SP_MAX_DIST:
				pData->X.push_back(pModel->multShoot_getSlackVarVal(*pData, con));
				pData->slackAssignCon.push_back(c);
				pData->numSlack++;
			case Constraint_tp::SP_DIST:
				addToRows = 1;
				break;
			case Constraint_tp::MAX_DIST:
			case Constraint_tp::MIN_DIST:
			case Constraint_tp::ENDSEG_MAX_DIST:
			case Constraint_tp::ENDSEG_MIN_DIST:
				pData->X.push_back(pModel->multShoot_getSlackVarVal(*pData, con));
				// remember where this slack variable is hiding
				pData->slackAssignCon.push_back(c);	
				pData->numSlack++;
				// do NOT break here, continue on to do stuff for DIST as well
			case Constraint_tp::DIST:
			case Constraint_tp::ENDSEG_DIST:
				addToRows = 1;
				break;
			case Constraint_tp::MAX_DELTA_V:
			case Constraint_tp::DELTA_V:
				if(!foundDVCon){
					addToRows = 1;
					foundDVCon = true;

					if(con.getType() == Constraint_tp::MAX_DELTA_V){
						/* Add a slack variable to the end of design vector and 
						 * keep track of which constraint pData is assigned to; 
						 * value of slack variable will be recomputed later
						 */
						pData->X.push_back(pModel->multShoot_getSlackVarVal(*pData, con));
						pData->numSlack++;
						pData->slackAssignCon.push_back(c);
					}
				}else{
					throw Exception("MultShootEngine::multShoot: "
						"You can only apply ONE delta-V constraint");
				}
				break;
			case Constraint_tp::ANGLE:
			case Constraint_tp::ENDSEG_ANGLE:
			case Constraint_tp::APSE:
			case Constraint_tp::ENDSEG_APSE:
			case Constraint_tp::JC:
			case Constraint_tp::ENDSEG_JC:
			case Constraint_tp::HLT:
			case Constraint_tp::PSEUDOARC:
			case Constraint_tp::SEG_CONT_EX:
				addToRows = 1;
				break;
			case Constraint_tp::TOF_TOTAL:
				if(to_underlying(tofTp) <= 0)
					astrohelion::printWarn("MultShootEngine::multShoot: "
						"Attempting to constraint TOF without variable time... "
						"won't work!\n");
				
				if(!foundTOFCon)
					addToRows = 1;
				else
					throw Exception("MultShootEngine::multShoot: "
						"You can only apply ONE TOF constraint");
				break;
			case Constraint_tp::EPOCH:
				if(to_underlying(tofTp) <= 0)
					printWarn("MultShootEngine::multShoot: "
						"Attempting to constrain Epoch without variable time..."
						" won't work!\n");

				addToRows = 1;
				break;
			case Constraint_tp::RM_STATE:
			case Constraint_tp::RM_EPOCH:
			case Constraint_tp::RM_CTRL:
				// These constraints are handled differently
				addToRows = 0;
				break;
			default:
				printWarn("MultShootEngine::multShoot: "
					"Unhandled constraint: %s\n", con.getTypeStr());
				break;
		}

		// Save the index of the first row for this constraint
		pData->conRows[c] = conRow;
		conRow += addToRows;	// remember we've added rows
		pData->totalCons += addToRows;
	}// END of loop through constraints

	// Determine the number of free/design variables based on the system type
	pData->totalFree = pData->X.size();

	// Save the initial free variable vector
	pData->X0 = pData->X;

	// Print debugging information
	printVerb(verbosity >= Verbosity_tp::ALL_MSG,
		"  # Free: %d\n  # Constraints: %d\n", pData->totalFree, pData->totalCons);
	printVerb(verbosity >= Verbosity_tp::ALL_MSG,
		"  -> # Slack Variables: %d\n", pData->numSlack);

	printVerb(verbosity >= Verbosity_tp::ALL_MSG, "ALL CONSTRAINTS:\n\n");
	if(verbosity >= Verbosity_tp::ALL_MSG){
		for(unsigned int n = 0; n < pData->allCons.size(); n++){
			pData->allCons[n].print();
		}
	}
	
	// Run the multiple shooting process
	try{
		multShoot(pData);
	}catch(const std::exception &e){
		// Free allocated data
		if(bDataOnStack){ delete pData; }
		throw e;	// Throw the error up a level
	}

	// Free allocated data
	if(bDataOnStack){ delete pData; }
}//==========================================================

/**
 *  @brief Run a multiple shooting algorithm given an MultShootData object
 * 
 *  @param it A completely formed MultShootData object that describes a 
 *  multiple shooting problem. These are created from Arcset and its
 *  derivative types by the other implementation of multShoot()
 *  @return A corrected MultShootData object
 *  @see multShoot(Arcset*)
 *  @throws DivergeException if the multiple shooting process does not converge
 */
void MultShootEngine::multShoot(MultShootData *pData){
	pData->count = 0;

	// create a simulation engine
	SimEngine simEngine;
	if(verbosity >= Verbosity_tp::DEBUG)
		simEngine.setVerbosity(static_cast<Verbosity_tp>(static_cast<int>(verbosity) - 1));
	else
		simEngine.setVerbosity(Verbosity_tp::NO_MSG);
	
	// Set both tolerances of simulation engine to be three orders of 
	// magnitude less corrector
	double simTol = tolF/1000 < 1e-15 ? 1e-15 : tolF/1000;
	simEngine.setAbsTol(simTol);
	simEngine.setRelTol(simTol);

	// Only need info about the final state, no need to generate lots of 
	// intermediate data points. This forces the integrator to use the 
	// Adams-Bashforth Adams-Moulton integration method, which is most similar 
	// to Matlab's ode113
	simEngine.setVarStepSize(false);
	simEngine.setNumSteps(2);

	if(bFindEvent || bIgnoreCrash){
		// don't use crash events when searching for an event
		simEngine.setMakeDefaultEvents(false);
	}

	// Define values for use in corrections loop
	double errF = 10*tolF, errX = 10*tolX, errF_infty = 10*tolF;
	unsigned int coreStateSize = 
		pData->pArcIn->getSysData()->getDynamicsModel()->getCoreStateSize();

	Eigen::VectorXd oldX(pData->totalFree, 1), newX(pData->totalFree, 1), 
		FX(pData->totalCons, 1);

	// Loop through iterations; error checking is done inside the loop
	while(pData->count < maxIts){
		if(pData->count > 0){
			// Solve for newX and copy into working vector X
			oldX = Eigen::Map<Eigen::VectorXd>(&(pData->X[0]), pData->totalFree, 1);
			
			try{
				solveUpdateEq(*pData, oldX, FX, newX);
			}catch(LinAlgException &e){
				throw e;	// Rethrow error
			}

			pData->X.clear();
			pData->X.insert(pData->X.begin(), newX.data(), newX.data()+pData->totalFree);
		}

		pData->FX.clear();					// Clear vectors each iteration
		pData->FX.assign(pData->totalCons, 0);	// Size the vectors and fill with zeros

		pData->DF_elements.clear();
		pData->DF_elements.reserve(pData->totalCons * coreStateSize);

		// Fill each trajectory object with a propagated arc
		propSegsFromFreeVars(*pData, simEngine);

		// Save the updated solution to file (for debugging)
		if(bSaveEachIt){
			pData->pArcIn->getSysData()->getDynamicsModel()->multShoot_createOutput(*pData);
			char filename[24];
			sprintf(filename, "multShoot_it%04d.mat", pData->count+1);
			pData->pArcOut->saveToMat(filename);
			pData->pArcOut->reset();	// Must be empty for next createOutput() call
		}
		// waitForUser();

		// Loop through all constraints and compute the constraint values, 
		// partials, and apply them to the FX and DF matrices
		for(unsigned int c = 0; c < pData->allCons.size(); c++){
			pData->pArcIn->getSysData()->getDynamicsModel()->\
				multShoot_applyConstraint(*pData, pData->allCons[c], c);
			printVerb(verbosity >= Verbosity_tp::DEBUG, 
				"* Applying %s constraint\n", pData->allCons[c].getTypeStr());
		}

		// Check to see what the error is; if it's too high, update X and continue another iteration
		errF_infty = std::abs(pData->FX[0]);
		for(unsigned int i = 1; i < pData->FX.size(); i++){
			if(std::abs(pData->FX[i]) > errF_infty)
				errF_infty = std::abs(pData->FX[i]);
		}

		FX = Eigen::Map<Eigen::VectorXd>(&(pData->FX[0]), pData->totalCons, 1);
		errF = FX.norm();

		if(verbosity >= Verbosity_tp::DEBUG)
			reportConMags(*pData);

		// Compute error: difference between subsequent free variable vectors
		if(pData->count > 0){
			Eigen::VectorXd diff = newX - oldX;
			errX = diff.norm();
		}

		pData->count++;
		printVerbColor((bFindEvent && verbosity >= Verbosity_tp::ALL_MSG) || 
			(!bFindEvent && verbosity > Verbosity_tp::NO_MSG),
			YELLOW, "It %02d : ||F(X)||_2 = %6.4e / %4.2e : "
			"||F(X)||_inf = %6.4e / %4.2e : "
			"||dX||_2 = %6.4e / %4.2e\n", pData->count, errF, tolF, 
			errF_infty, tolF, errX, tolX);

		// End the iterations if the constraint error grows too large or if
		// the constraint error or step size reaches the desired precision
		if(errF > maxErr || errF < tolF || errX < tolX)
			break;
	}// end of corrections loop

	if(!bIgnoreDiverge){
		if( (errF > maxErr) || (errF > tolF && errX > tolX)){
			char msg[256];
			if(errF > maxErr){
				sprintf(msg, "MultShootEngine: Diverged! "
					"Constraint Error = %e > maxErr = %e", errF, maxErr);
			}else{
				sprintf(msg, "MultShootEngine: Diverged! "
					"Constraint Error = %e > tolF = %e", errF, tolF);
			}

			throw DivergeException(msg);
		}
	}

	if(pData->pArcOut){
		try{
			// Save propagated data and free variable vector values to the 
			// output arcset
			if(bFullFinalProp){
				simEngine.setVarStepSize(true);
				propSegsFromFreeVars(*pData, simEngine, verbosity);
			}
			pData->pArcIn->getSysData()->getDynamicsModel()->multShoot_createOutput(*pData);
		}catch(Exception &e){
			printErr("MultShootEngine::multShoot: "
				"Unable to create output arcset\n  Err: %s\n", e.what());
			throw e;
		}
	}
}//=====================================================

/**
 *  @brief Propagate all segments along the arcset
 *  @details Initial conditions and other parameters for each integrated
 *  arc are obtained from the free variable vector or the input arcset
 * 
 *  @param it reference to the multiple shooting data structure
 *  @param sim reference to a simulalation engine initialized for multiple
 *  shooting propagations
 *  @param verbosity how verbose the output should be; the verbosity of the 
 *  input SimEngine object is set separately (i.e., via sim.setVerbosity()) and 
 *  is not affected by this input
 */
void MultShootEngine::propSegsFromFreeVars(MultShootData& it, SimEngine &sim, 
	Verbosity_tp verbosity){
	
	unsigned int coreStateSize = it.pArcIn->getSysData()->getDynamicsModel()->\
		getCoreStateSize();

	it.deltaVs.clear();
	it.deltaVs.assign(3*it.numNodes, 0);

	// initialize a vector of trajectory objects to store each propagated 
	// segment
	it.propSegs.clear();
	it.pArcIn->getSysData()->getDynamicsModel()->multShoot_initIterData(it);

	for(unsigned int s = 0; s < it.pArcIn->getNumSegs(); s++){
		// printf("Retrieving ICs for segment (ix %02d):\n", s);
		// Get simulation conditions from design vector via dynamic model 
		// implementation
		double t0 = 0, tof = 0;
		std::vector<double> ic(coreStateSize, 0);

		ControlLaw *pLaw = it.pArcIn->getSegRefByIx_const(s).getCtrlLaw();
		std::vector<double> ctrl0;
		if(pLaw){
			ctrl0.assign(pLaw->getNumStates(), 0);
		}

		it.pArcIn->getSysData()->getDynamicsModel()->multShoot_getSimICs(it, 
			it.pArcIn->getSegRefByIx_const(s).getID(),
			&(ic.front()), &(ctrl0.front()), &t0, &tof);

		sim.setRevTime(tof < 0);
		
		printVerbColor(verbosity >= Verbosity_tp::DEBUG, MAGENTA, 
			"Simulating segment %d:\n  t0 = %.4f\n  tof = %.4f\n", s, t0, tof);

		try{
			sim.runSim(ic, ctrl0, t0, tof, &(it.propSegs[s]), pLaw);
		}catch(DivergeException &e){
			if(verbosity >= Verbosity_tp::SOME_MSG){
				printf("%s", RED);
				printf("SimEngine integration diverged on segment %u...\n", s);
				printf("  > ic = [");
				for(unsigned int i = 0; i < ic.size(); i++){
					printf(" %.4f ", ic[i]);
				}
				printf("]\n  > t0 = %.4f\n  > tof = %.4f\n", t0, tof);
			}
		}catch(Exception &e){
			printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, RED, 
				"SimEngine Error on segment %u:\n%s\nEnding corrections.\n",
				s, e.what());
		}

		// if(verbosity >= Verbosity_tp::DEBUG){
		// 	it.propSegs[s].print();
		// }

		std::vector<double> lastState = it.propSegs[s].getStateByIx(-1);
		int termID = it.pArcIn->getSegRefByIx_const(s).getTerminus();
		if(termID != Linkable::INVALID_ID){

			// Get the state of the terminal node
			MSVarMap_Obj stateVar = it.getVarMap_obj(MSVar_tp::STATE, termID);
			std::vector<double> state;
			if(stateVar.row0 == -1){
				state = it.pArcIn->getState(termID);
			}else{
				state = std::vector<double>(it.X.begin() + stateVar.row0, 
					it.X.begin() + stateVar.row0 + stateVar.nRows);
			}

			// velCon has false for a velocity state if there is a 
			// discontinuity between the terminus of the segment and the 
			// terminal node
			std::vector<bool> velCon = it.pArcIn->getSegRefByIx_const(s).getVelCon();
			for(int i = 3; i < 6; i++){
				// Compute difference in velocity; if velCon[i-3] is true, 
				// then velocity should be continuous and any difference is 
				// numerical error, so set to zero by multiplying by not-true
				it.deltaVs[s*3+i-3] = !velCon[i-3]*(lastState[i] - state[i]);
			}
		}
	}
}//=====================================================

/**
 *	@brief Apply linear algebra to solve the update equation and obtain an 
 *	updated free-variable vector.
 *
 *	@details The update equation takes the form
 * 	\f[
 *		\vec{F}(\vec{X}) = D\vec{F}(\vec{X}) \left( \vec{X}_{n+1} - 
 *		\vec{X}_n \right)
 *	\f]
 *	We know \f$ \vec{X}_n \f$, \f$ \vec{F}(\vec{X}) \f$ and 
 *	\f$ D\vec{F}(\vec{X}) \f$, and we wish to solve for
 *	\f$ \vec{X}_{n+1} \f$. To do this, we need to invert or factor the 
 *	Jacobian matrix \f$ D\vec{F}(\vec{X}) \f$. Assuming the Jacobian is 
 *	non-singular, a solution can be obtained. Three possible scenarios can 
 *	occur: either the system is over-constrained (no solution), perfectly 
 *	constrained (square, one solution), or under constrained (infinitely many 
 *	solutions).
 *
 *	In the first case, we can use least squares to find the closest possible 
 *	solution. This is not implemented and the function will throw an error of 
 *	the system is over constrained. In the second case, the equation is solved
 *	via LU factorization using GSL's linear algebra functions. Finally, if the 
 *	system is under constrained, which is the most common case, we compute the 
 *	minimum-norm solution. We sequentially solve the following equations to
 *	arrive at the min-norm solution:
 *	\f{eqnarray*}{
 *		JJ^T \vec{w} &=& \vec{F}(\vec{X}) \\
 *		\vec{X}^* &=& J^T \vec{w}
 *	\f}
 *	where \f$ J \f$ is the Jacobian matrix, \f$ JJ^T \f$ is the associated 
 *	Gram Matrix, and \f$ \vec{X}^* \f$ is the min-norm solution. If the 
 *	Jacobian is non-singular, then the Gram Matrix will also be non-singular 
 *	and the system can be solved. Note that \f$ \vec{X}^* \f$ is the min-norm 
 *	solution for \f$ \vec{X}_{n+1} = \vec{X}_n \f$.
 *
 *	In all cases, errors will be thrown if the Jacobian is singular. This most 
 *	likely indicates that there has been a coding error in the corrector, 
 *	although singular Jacobians do occur when trajectories pass very near 
 *	primaries.
 *
 *	@param it reference to the MultShootData object associated with the 
 *	corrections process
 *	@param oldX constant reference to the current (previous) design vector
 *	@param FX constant reference to the current constraint vector
 *	@param newX reference to a vector in which to store the updated design 
 *	variable vector
 *	
 *	@throws LinAlgException if the problem is over constrained (i.e. Jacobian 
 *	has more rows than columns);
 *	
 *	@todo This can be updated to use a least-squares solution
 */
void MultShootEngine::solveUpdateEq(MultShootData& it, const Eigen::VectorXd& oldX, const Eigen::VectorXd& FX, Eigen::VectorXd& newX){
	if(it.totalCons == 1 && it.totalFree == 1){
		// If Jacobian is 1x1, skip all that linear algebra and just solve the 
		// equation. Update the design variable vector
		newX = oldX - FX/(it.DF_elements[0].value());
	}else{
		// Jacobian is not scalar; use linear algebra to solve update equation

		Eigen::VectorXd fullStep(it.totalFree, 1);	// Vector for solution
		
		// Construct the sparse Jacobian matrix
		SparseMatXCd J(it.totalCons, it.totalFree);
		J.setFromTriplets(it.DF_elements.begin(), it.DF_elements.end());
		
		if(it.totalCons == it.totalFree){	// J is square, use regular inverse

			/* Use LU decomposition to invert the Jacobian matrix and find a vec
			w. Multiplying J^T by w yields the minimum-norm solution x, where x 
			lies in the column-space of J^T, or in the orthogonal complement of
			the nullptrspace of J.
			Source: <http://www.math.usm.edu/lambers/mat419/lecture15.pdf>
			 */
			
			// Info about solving Sparse matrix systems with Eigen:
			// <https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html>

			// Solve the system Jw = b (In this case, w = fullStep)
			J.makeCompressed();
			factorizeJacobian(J, FX, fullStep, false);

			// For whatever reason, must cast to non-references?
			// MatrixXRd tempJ(J), tempFX(FX), tempX(oldX), tempStep(fullStep);
			// astrohelion::toCSV(tempJ, "DF_cpp.csv");
			// astrohelion::toCSV(tempFX, "FX_cpp.csv");
			// astrohelion::toCSV(tempX, "X_cpp.csv");
			// astrohelion::toCSV(tempStep, "dX_cpp.csv");
			// waitForUser();

		}else{
			if(it.totalCons < it.totalFree){	// Under-constrained

				/* Use LU decomposition to invert the Gramm matrix and find a 
				vector w. Multiplying J^T by w yields the minimum-norm solution
				x, where x lies in the column-space of J^T, or in the 
				orthogonal complement of the nullptrspace of J.
				Source: <http://www.math.usm.edu/lambers/mat419/lecture15.pdf>
				 */

				// Compute Gramm matrix
				SparseMatXCd JT = J.transpose();
				SparseMatXCd G = J*JT;		// G will always be symmetric

				// Solve the system Gw = b
				G.makeCompressed();
				Eigen::VectorXd w(G.cols(), 1);
				factorizeJacobian(G, FX, w, true);
				fullStep = JT*w;

				MatrixXRd mJ(J), mFX(FX), mStep(fullStep);
				toCSV(mJ, "DF_cpp.csv");
				toCSV(mFX, "FX_cpp.csv");
				toCSV(mStep, "dX_cpp.csv");
				waitForUser();
				// toCSV(oldX, "X_cpp.csv");

				// Alternative Method: SVD
				// NOTE: This takes approximately five times as much 
				// computation time
				// JacobiSVD<MatrixXd> svd(JT, Eigen::ComputeFullV | Eigen::ComputeFullU);
				// svd.setThreshold(tolF/100.f);

				// MatrixXd singVals = svd.singularValues();
				// // First singular value is the biggest one
				// double svdTol = (1e-12)*singVals(0);

				// MatrixXd Z = MatrixXd::Zero(JT.rows(), JT.cols());
				// for(unsigned int r = 0; r < singVals.rows(); r++){
				// 	Z(r,r) = std::abs(singVals(r)) > svdTol ? 1.0/singVals(r) : 0;
				// }

				// fullStep = svd.matrixU() * Z * svd.matrixV().transpose() * FX;

			}else{	// Over-constrained
				throw LinAlgException("MultShootEngine::solveUpdateEq: "
					"System is over constrained... No solution implemented");
			}			
		}

		if(bLineSearchAttenFactor){
			bool checkLocalMin = false;
			chooseStep_LineSearch(it, oldX, FX, J, fullStep, newX, 
				checkLocalMin);

			// This algorithm to check for the local min is lifted directly 
			// from the NumericalRecipes book, section 9.7 in the `newt` 
			// function example
			if(checkLocalMin){
				// f = half the magnitude of the constraint vector
				double f = 0.5*FX.transpose()*FX;
				double test = 0, temp = 0, denom = std::max(f, 0.5*FX.rows());
				Eigen::RowVectorXd grad_f = FX.transpose() * J;
				for(unsigned int i = 0; i < FX.rows(); i++){
					temp = std::abs(grad_f(i)) * std::max(newX(i), 1.0) / denom;
					if(temp > test)
						test = temp;
				}

				if(test < tolA){
					throw LinAlgException("MultShootEngine:solveUpdateEq: "
						"The line search has converged to a local minimum of "
						"f = (1/2) FX * FX. Try another initial guess for the "
						"design vector");
				}
			}
		}else{
			newX = oldX + fullStep;	// Update the design variable vector
		}
	}
}// End of solveUpdateEq() ============================

/**
 *  @brief Apply LU and/or QR factorization to solve the update equation
 *  @details LU factorization is attempted first; if this fails, QR 
 *  factorization is attempted as it can be more robust. The details of the 
 *  factorization process are explained in the Eigen documentation at
 *  <https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html>
 * 
 *  @param J the Jacobian (or Gramm) matrix that needs to be factorized
 *  to solve the update equation. J *must* be in compressed form!
 *  @param FX constant reference to the current constraint vector
 *  @param out the solution to the equation J*out = FX
 *  @param bIsSymmetric whether or not the matrix J is symmetric
 */
void MultShootEngine::factorizeJacobian(const SparseMatXCd &J, 
	const Eigen::VectorXd& FX, Eigen::VectorXd &out, bool bIsSymmetric){

	Eigen::ComputationInfo info;

	// If the LU solver fails, it is almost definitely because the
	// Jacobian is poorly scaled. QR performs better w/ poor scaling,
	// So use QR for the rest of the iterations
	if(bLUFailed){
		// use QR
		info = QR(J, FX, out);
	}else{
		// use LU
		info = LU(J, FX, out, bIsSymmetric);

		// If LU fails, use QR and remember that LU failed
		if(info != Eigen::Success){
			bLUFailed = true;
			printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, RED, "LU"
				" Factorization failed, trying QR factorization to solve "
				"update equation.\n");
			info = QR(J, FX, out);
		}
	}

	// Report any errors to the user
	if(info != Eigen::Success){
		checkDFSingularities(J);
		std::string err = eigenCompInfo2Str(info);
		char msg[256];
		sprintf(msg, "MultShootEngine::factorizeJacobian: "
			"Could not factorize Jacobian matrix.\nEigen error: %s\n",
			err.c_str());
		throw LinAlgException(msg);
	}
}//====================================================

/**
 *  @brief Apply QR factorization to solve the update equation
 *  @details The details of the factorization process
 *  are explained in the Eigen documentation at
 *  <https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html>
 * 
 *  @param J the Jacobian (or Gramm) matrix that needs to be factorized
 *  to solve the update equation. J *must* be in compressed form!
 *  @param FX constant reference to the current constraint vector
 *  @param out the solution to the equation J*out = FX
 */
Eigen::ComputationInfo MultShootEngine::QR(const SparseMatXCd &J, 
	const Eigen::VectorXd& FX, Eigen::VectorXd &out){

	// The structure of J is the same for every iteration, so only 
	// analyze pattern once
	if(!bInitQRSolver){
		qrSolver.analyzePattern(J);
		bInitQRSolver = true;
	}

	qrSolver.factorize(J);

	// If factorization was successful, solve the equation
	if(qrSolver.info() == Eigen::Success)
		out = qrSolver.solve(-FX);

	return qrSolver.info();
}//====================================================

/**
 *  @brief Apply LU factorization to solve the update equation
 *  @details The details of the factorization process
 *  are explained in the Eigen documentation at
 *  <https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html>
 * 
 *  @param J the Jacobian (or Gramm) matrix that needs to be factorized
 *  to solve the update equation. J *must* be in compressed form!
 *  @param FX constant reference to the current constraint vector
 *  @param out the solution to the equation J*out = FX
 *  @param bIsSymmetric whether or not the matrix J is symmetric
 */
Eigen::ComputationInfo MultShootEngine::LU(const SparseMatXCd &J, 
	const Eigen::VectorXd& FX, Eigen::VectorXd &out, bool bIsSymmetric){

	// The structure of J is the same for every iteration, so only 
	// analyze pattern once
	if(!bInitLUSolver){
		luSolver.analyzePattern(J);
		luSolver.isSymmetric(bIsSymmetric);
		bInitLUSolver = true;
	}

	luSolver.factorize(J);

	// If factorization was successful, solve the equation
	if(luSolver.info() == Eigen::Success)
		out = luSolver.solve(-FX);

	return luSolver.info();
}//====================================================

/**
 *  @brief Use a line search to scale the Newton step to ensure convergence
 *  @details Leveraging this method vastly improves the convergence behavior
 *  of the standard Newton-Raphson solver by scaling the default Newton step
 *  to ensure that each iteration is closer to the solution than the last. 
 *  
 *  This method uses the line search and backtracking algorithm detailed in 
 *  section 9.7, "Globally Convergent Methods for Nonlinear Systems of Equations"
 *  in the book, Numerical Recipes in C: The Art of Scientific Computing. The 
 *  implementation details are included in my MathSpec document.
 *  
 *  The algorithm uses information from the current and previous iteration to 
 *  scale the full Newton step and compute a new design vector, which is stored
 *  in the vector that newX references
 * 
 *  @param it reference to the current iteration data object. This stores the 
 *  current constraint vector elements, Jacobian matrix elements, and design 
 *  vector elements, among other variables.
 *  @param oldX reference to the previous iteration's design vector
 *  @param oldFX reference to the previous iteration's constraint vector
 *  @param J reference to the Jacobian matrix, in compressed form
 *  @param fullStep reference to the full Newton step computed in 
 *  solveUpdateEq()
 *  @param newX pointer to a vector to store the new design vector in. 
 *  @param checkLocalMin pointer to a boolean that is set to true if 
 *  X_old and X_new are too close to one another. In a minimization algorithm 
 *  this usually signals convergence and can be ignored. In a root-finding 
 *  algorithm, this should trigger additional checks.
 */
void MultShootEngine::chooseStep_LineSearch(MultShootData& it, const Eigen::VectorXd& oldX, const Eigen::VectorXd& oldFX,
	const SparseMatXCd &J, const Eigen::VectorXd& fullStep, Eigen::VectorXd& newX, bool& checkLocalMin){

	double f_old = 0.5*oldFX.norm();
	checkLocalMin = false;

	// Initialize full step
	Eigen::VectorXd fullStep_scaled = fullStep;
	
	// Scale if attempted step is too big (i.e., if there is some unbounded thing going on)
	if(fullStep_scaled.norm() > ls_maxStepSize){
		fullStep_scaled *= ls_maxStepSize/fullStep_scaled.norm();
	}

	// Compute the gradient of the constraint vector magnitude; this is the initial rate of decrease
	Eigen::RowVectorXd slopeProd = oldFX.transpose() * J;
	slopeProd = slopeProd * fullStep_scaled;
	assert(slopeProd.rows() == 1 && slopeProd.cols() == 1);
	double initROD = slopeProd(0);	// initial rate of decrease

	if(initROD >= 0.0){
		printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, RED, 
			"  Line search: initROD = %e >= 0... roundoff error in update\n",
			initROD);
		
		// Use the analytical initial rate of descent by bypassing the update equation
		initROD = -1 * oldFX.transpose() * oldFX;
		printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, BLUE, 
			"  Line Search: Using analytical solution for initROD = %e\n",
			initROD);
		// throw DivergeException(msg);
	}

	// Compute min step size
	double big = 0, temp = 0;
	for(unsigned int i = 0; i < fullStep_scaled.rows(); i++){
		temp = std::abs(fullStep_scaled(i))/std::max(std::abs(oldX(i)), 1.0);
		if(temp > big)
			big = temp;
	}
	double minLambda = tolX/big;	// smallest permissible attenuation factor
	double lambda = 1;				// attenuation factor begins at 1

	printVerbColor(verbosity >= Verbosity_tp::ALL_MSG, BLUE, 
		"  Line Search: Minimum attenuation factor = %.4e\n", minLambda);

	// Set up a simulation engine to propagate segments
	SimEngine sim;
	double simTol = tolF/1000 < 1e-15 ? 1e-15 : tolF/1000;
	sim.setAbsTol(simTol);
	sim.setRelTol(simTol);
	sim.setVarStepSize(false);
	sim.setNumSteps(2);
	sim.setMakeDefaultEvents(false);
	if(verbosity >= Verbosity_tp::DEBUG)
		sim.setVerbosity(static_cast<Verbosity_tp>(static_cast<int>(verbosity) - 1));
	else
		sim.setVerbosity(Verbosity_tp::NO_MSG);
	
	unsigned int maxCount = 10;		// Max number of line search iterations
	double lambda_next = 1;			// Storage for the next step size
	double lambda_prev = 1;			// Step size from the previous iteration of the line search
	double f = 1;					// Error term for the current iteration
	double f2 = 1;					// Error term from the previous iteration of the line search

	unsigned int coreDim = it.pArcIn->getSysData()->getDynamicsModel()->getCoreStateSize();

	for(unsigned int count = 0; count < maxCount; count++){
		// Compute next free variable vector give the current attenuation factor
		newX = oldX + lambda * fullStep_scaled;

		// Compute an updated constraint vector
		MultShootData tempData(it);
		tempData.FX.clear();
		tempData.FX.assign(tempData.totalCons, 0);
		tempData.DF_elements.clear();
		tempData.DF_elements.reserve(tempData.totalCons * coreDim);

		// Update free variable vector in the temporary data object
		assert(newX.rows() > 1 && newX.cols() == 1);
		tempData.X = std::vector<double>(newX.data(), newX.data() + newX.rows());
		
		// Do all those expensive function evaluations
		propSegsFromFreeVars(tempData, sim, verbosity);
		for(unsigned int c = 0; c < tempData.allCons.size(); c++){
			tempData.pArcIn->getSysData()->getDynamicsModel()->multShoot_applyConstraint(tempData, tempData.allCons[c], c);
		}

		// Compute magnitude of new FX vector
		temp = 0;
		for(unsigned int i = 0; i < tempData.FX.size(); i++){
			temp += tempData.FX[i]*tempData.FX[i];
		}
		// Scaled norm of the constraint vector (an error metric)
		f = 0.5*temp;

		if(lambda < minLambda){
			// Apparently, the "optimal" attenuation factor is approaching zero
			// Thus, the update delta X is zero, i.e., we've converged on X_old
			printVerbColor(verbosity >= Verbosity_tp::ALL_MSG, RED, 
				"  Line search decreased past minimum bound.\n");
			newX = oldX;
			checkLocalMin = true;
			return;
		}


		if(f <= f_old + ls_alpha*lambda*initROD){
			printVerbColor(verbosity >= Verbosity_tp::ALL_MSG, BLUE, 
				"  Line Search: attenuation factor = %.4e (%u its)\n", lambda,
				count);
			return;	// We're all set, so return!
		}else{	// Need to backtrack
			if(count == 0){
				// Approximate step-size as a quadratic, find attenuation factor that yields better guess
				lambda_next = -initROD/(2.0*(f - f_old - initROD));
			}else{
				// Approximate step-size as a cubic
				double term1 = f - f_old - lambda*initROD;
				double term2 = f2 - f_old - lambda_prev*initROD;
				double a = (term1/(lambda*lambda) - term2/(lambda_prev*lambda_prev))/(lambda - lambda_prev);
				double b = (-lambda_prev*term1/(lambda*lambda) + lambda*term2/(lambda_prev*lambda_prev))/(lambda - lambda_prev);

				if(a == 0){
					lambda_next = -initROD/(2.0*b);
				}else{
					double disc = b*b - 3*a*initROD;
					if(disc < 0){
						lambda_next = 0.5*lambda;
					}else if(b <= 0){
						lambda_next = (-b + sqrt(disc))/(3.0*a);
					}else{
						lambda_next = -initROD/(b + sqrt(disc));
					}

					// Ensure that attenuation factor decreases by at least a factor of 2
					if(lambda_next > 0.5*lambda){
						lambda_next = 0.5*lambda;
					}
				}
			}

			// Update storage for next iteration
			lambda_prev = lambda;
			f2 = f;

			// Ensure that attenuation factor decreases by no more than 90%
			lambda = std::max(lambda_next, 0.1*lambda);
		}
	}

	if(f <= f_old + ls_alpha*lambda*initROD)
		throw DivergeException("MultShootEngine::updateFreeVarVec: Line search for attenuation factor diverged.");

	printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, BLUE, "  Line Search: attenuation factor = %.4e (%d its)\n", lambda, maxCount);
}//====================================================

/**
 *  @brief Checks the Jacobian (DF) matrix for singularities, i.e., rows
 *  or columns that contain only zeros
 * 
 *  @param DF Jacobian matrix employed in the correction process
 */
void MultShootEngine::checkDFSingularities(MatrixXRd DF){
	if(verbosity == Verbosity_tp::NO_MSG)
		return;

	bool foundSingularity = false;
	if(verbosity >= Verbosity_tp::SOME_MSG){
		for(unsigned int r = 0; r < DF.rows(); r++){
			if(DF.row(r).norm() == 0){
				printErr("Singular Jacobian! Row %u is all zeros\n", r);
				foundSingularity = true;
			}
		}

		for(unsigned int c = 0; c < DF.cols(); c++){
			if(DF.col(c).norm() == 0){
				printErr("Singular Jacobian! Column %u is all zeros\n", c);
				foundSingularity = true;
			}
		}
	}

	if(!foundSingularity)
		printColor(BLUE, "No singularity found in Jacobian\n");
}//====================================================

/**
 *  @brief Print out the magnitude of each constraint.
 *  @details This can be useful when debugging to highlight which constraints are unsatisfied
 * 
 *  @param it reference to an MultShootData object associated with a corrections process
 */
void MultShootEngine::reportConMags(const MultShootData& it){
	// unsigned int conCount = 0;
	// for(unsigned int r = 0; r < (it.FX.size()); r++){
 //        if(r == 0 && it.totalCons > 0){
 //            printf("Node %d %s Constraint:\n", it.allCons[conCount].getID(), 
 //            	it.allCons[conCount].getTypeStr());
 //        }else if(conCount < it.allCons.size() && r >= 
 //        	static_cast<unsigned int>(it.conRows[conCount+1])){

 //            conCount++;
 //            printf("Node %d %s Constraint:\n", it.allCons[conCount].getID(), 
 //            	it.allCons[conCount].getTypeStr());
 //        }
 //        printf("  ||row %03u||: %.6e\n", r, std::abs(it.FX[r]));
 //    }
    unsigned int conCount = 0;
	for(unsigned int r = 0; r < it.FX.size(); r++){
        if(r == 0 && it.totalCons > 0){
            printf("Applies to %s %d: %s Constraint:\n", 
                Constraint::getAppTypeStr(it.allCons[conCount].getAppType()),
                it.allCons[conCount].getID(), it.allCons[conCount].getTypeStr());
        }else if(conCount+1 < it.allCons.size() && 
        	r >= static_cast<unsigned int>(it.conRows[conCount+1])){

            conCount++;
            printf("Applies to %s %d: %s Constraint:\n", 
                Constraint::getAppTypeStr(it.allCons[conCount].getAppType()),
                it.allCons[conCount].getID(), it.allCons[conCount].getTypeStr());
        }
        printColor(it.FX[r] > tolF ? RED : GREEN, "  row %03zu: %.6e \n", r, 
        	it.FX[r]);
    }

}//===============================================================

/**
 *	@brief clean up data so that engine can be used again (or deconstructed safely)
 */
void MultShootEngine::cleanEngine(){
	astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "Cleaning the engine...\n");
	
	bInitLUSolver = false;
	bInitQRSolver = false;
	bLUFailed = false;

	bIsClean = true;
}//====================================================

/**
 *  @brief Reset the multiple shooting engine to the default parameter values
 */
void MultShootEngine::reset(){
	if(!bIsClean)
		cleanEngine();

	tofTp = MSTOF_tp::VAR_FREE;
	maxIts = 20;
	tolF = 1e-12;
	tolX = 1e-14;
	tolA = 1e-12;
	bFindEvent = false;
	bIgnoreCrash = false;
	bIgnoreDiverge = false;
	bFullFinalProp = true;
	bSaveEachIt = false;
	ls_alpha = 1e-4;
	ls_maxStepSize = 100;
}//====================================================

/**
 *  @brief Check the DF matrix for the multiple shooting algorithm using finite differencing
 *  @details This function checks to make sure the Jacobian matrix (i.e. DF) is correct
 *  by computing each partial derivative numerically via forward differencing.
 * 
 *  @param pNodeset A arcset with some constraints
 *  @param verbosity Specify how verbose the output is
 *  @param writeToFile Whether or not to write .csv or .mat files with relevant information
 *  
 *  @return whether or not the Jacobian matrix is "correct". Issues that result in a return 
 *  value of FALSE are: The Jacobian is singular (row or column of all zeros), or the differences
 *  between the numerically computed Jacobian and analytical Jacobian are large.
 */
bool MultShootEngine::finiteDiff_checkMultShoot(const Arcset *pNodeset, Verbosity_tp verbosity, bool writeToFile){
    MultShootEngine engine;  // Create engine with default settings
    return finiteDiff_checkMultShoot(pNodeset, engine, verbosity, writeToFile);
}//====================================================

/**
 *  @brief Check the DF matrix for the multiple shooting algorithm using finite differencing
 *  @details This function checks to make sure the Jacobian matrix (i.e. DF) is correct
 *  by computing each partial derivative numerically via forward differencing.
 * 
 *  @param pNodeset A arcset with some constraints
 *  @param engine correction engine object configured with the appropriate settings (i.e.,
 *  equal arc time, etc.). Note that the maxIts, verbosity, and ignoreDiverge
 *  attributes of the engine will be overridden by this function.
 *  @param verbosity Specify how verbose the output is. If set to DEBUG, info about the
 *  errors will be printed regardless of their size.
 *  @param writeToFile Whether or not to write .csv or .mat files with relevant information
 *  
 *  @return whether or not the Jacobian matrix is "correct". Issues that result in a return 
 *  value of FALSE are: The Jacobian is singular (row or column of all zeros), or the differences
 *  between the numerically computed Jacobian and analytical Jacobian are large.
 */
bool MultShootEngine::finiteDiff_checkMultShoot(const Arcset *pNodeset, 
	MultShootEngine engine, Verbosity_tp verbosity, bool writeToFile){

    printVerb(verbosity >= Verbosity_tp::SOME_MSG, 
    	"Finite Diff: Checking DF matrix...\n");
    // Create multiple shooter that will only do 1 iteration
    MultShootEngine corrector(engine);
    corrector.setMaxIts(1);
    corrector.setVerbosity(verbosity < Verbosity_tp::DEBUG ? 
    	Verbosity_tp::NO_MSG : Verbosity_tp::ALL_MSG);
    corrector.setIgnoreDiverge(true);
    corrector.setIgnoreCrash(true);
    corrector.setFullFinalProp(false);

    // Run multiple shooter to get X, FX, and DF
    MultShootData it(pNodeset);
    corrector.multShoot(pNodeset, nullptr, &it);
    Eigen::VectorXd FX = Eigen::Map<Eigen::VectorXd>(&(it.FX[0]), it.totalCons, 1);
    SparseMatXCd sparseDF(it.totalCons, it.totalFree);
    sparseDF.setFromTriplets(it.DF_elements.begin(), it.DF_elements.end());
    MatrixXRd DF = MatrixXRd(sparseDF);
    // MatrixXRd DF = Eigen::Map<MatrixXRd>(&(it.DF[0]), it.totalCons, it.totalFree);
    MatrixXRd DFest = MatrixXRd::Zero(it.totalCons, it.totalFree);

    if(writeToFile){
        astrohelion::toCSV(DF, "FiniteDiff_DF.csv");
    }

    double pertSize = 1e-8;
    #pragma omp parallel for firstprivate(it, corrector) schedule(dynamic)
    for(int i = 0; i < it.totalFree; i++){
        std::vector<double> pertX = it.X0;      // Copy unperturbed state vetor
        pertX[i] += pertSize;                   // add perturbation
        it.X = pertX;                           // Copy into iteration data
        MultShootData pertIt = it;
        corrector.multShoot(&pertIt);     // Correct perturbred state
        Eigen::VectorXd FX_up = Eigen::Map<Eigen::VectorXd>(&(pertIt.FX[0]), 
        	it.totalCons, 1);

        // Do another process for opposite direction
        pertX = it.X0;
        pertX[i] -= pertSize;
        it.X = pertX;
        pertIt = it;
        corrector.multShoot(&pertIt);
        Eigen::VectorXd FX_down = Eigen::Map<Eigen::VectorXd>(&(pertIt.FX[0]), 
        	it.totalCons, 1);

        // An iteration for twice the perturbation up
        pertX = it.X0;
        pertX[i] += 2*pertSize;
        it.X = pertX;
        pertIt = it;
        corrector.multShoot(&pertIt);
        Eigen::VectorXd FX_2up = Eigen::Map<Eigen::VectorXd>(&(pertIt.FX[0]), 
        	it.totalCons, 1);

        // An iteration for twice the perturbation down
        pertX = it.X0;
        pertX[i] -= 2*pertSize;
        it.X = pertX;
        pertIt = it;
        corrector.multShoot(&pertIt);
        Eigen::VectorXd FX_2down = Eigen::Map<Eigen::VectorXd>(&(pertIt.FX[0]), 
        	it.totalCons, 1);


        // Compute central difference
        Eigen::VectorXd col = (-1*FX_2up + 8*FX_up - 8*FX_down + FX_2down)/
        	std::abs(12*pertSize);   // Five-point stencil
        // Eigen::VectorXd col = (FX_up - FX_down)/std::abs(2*pertSize);   // Central Difference
        // Eigen::VectorXd col = (FX_up - FX)/std::abs(pertSize);       // Forward difference
        DFest.block(0, i, it.totalCons, 1) = col;
    }


    MatrixXRd diff = DF - DFest;
    MatrixXRd DF_abs = DF.cwiseAbs();       // Get coefficient-wise absolute value
    MatrixXRd DFest_abs = DFest.cwiseAbs();

    diff = diff.cwiseAbs();	// Get coefficient-wise aboslute value
    MatrixXRd relDiff = MatrixXRd::Zero(diff.rows(), diff.cols());	// Relative differences (divide by magnitude of each value)

    // Divide each element by the magnitude of the DF element to get a relative difference magnitude
    for(int r = 0; r < diff.rows(); r++){
        for(int c = 0; c < diff.cols(); c++){
            // If one of the elements is zero, let the difference just be the difference; no division
            if(DF_abs(r,c) > 1e-13 && DFest_abs(r,c) > 1e-13)   // consider 1e-13 to be zero
                relDiff(r,c) = diff(r,c)/DF_abs(r,c);
        }
    }
    
    if(writeToFile){
        toCSV(DFest, "FiniteDiff_DFest.csv");
        toCSV(diff, "FiniteDiff_Diff.csv"); 
        toCSV(relDiff, "FiniteDiff_RelDiff.csv");
    }

    // Compute the largest coefficient in each row and column
    std::vector<unsigned int> rowMaxIndex(diff.rows(), 0), colMaxIndex(diff.cols(), 0);
    Eigen::VectorXd rowMax(diff.rows(), 1);
    Eigen::VectorXd colMax(diff.cols(), 1);
    for(unsigned int r = 0; r < diff.rows(); r++){
    	rowMax(r) = diff.row(r).maxCoeff(&rowMaxIndex[r]);
    }
	for(unsigned int c = 0; c < diff.cols(); c++){
		colMax(c) = diff.col(c).maxCoeff(&colMaxIndex[c]);
	}

    // Eigen::VectorXd rowMax = diff.rowwise().maxCoeff();
    // Eigen::RowVectorXd colMax = diff.colwise().maxCoeff();

    // Make a map that links the row number (column in the Jacobian matrix)
    // of the free variable to the MSVarMap_Key that represents the MSVarMap_Obj 
    // that contains information about the free variable
    std::map<int, MSVarMap_Key> freeVarMap_rowNumToKey;
    for(auto& obj : it.freeVarMap){
    	if(obj.second.row0 != -1){
	        for(int r = 0; r < obj.second.nRows; r++){
	        	freeVarMap_rowNumToKey[obj.second.row0 + r] = obj.first;
	        }
	    }
    }

    unsigned int rowMaxMaxIx = 0, colMaxMaxIx = 0;
    double rowMaxMax = rowMax.maxCoeff(&rowMaxMaxIx);
    double colMaxMax = colMax.maxCoeff(&colMaxMaxIx);
    int errScalar = 100000;

    bool goodDF = true;
    for(unsigned int r = 0; r < DF.rows(); r++){
    	if(DF.row(r).norm() == 0){
    		printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, BOLDRED,
    			"Singular Jacobian: row %u contains only zeros\n", r);

    		goodDF = goodDF && false;
    	}
    }

    for(unsigned int c = 0; c < DF.cols(); c++){
    	if(DF.col(c).norm() == 0){
    		printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, BOLDRED,
    			"Singular Jacobian: column %u contains only zeros\n", c);

    		goodDF = goodDF && false;
    	}
    }

    if(rowMaxMax < errScalar*pertSize && colMaxMax < errScalar*colMaxMax){
        if(verbosity >= Verbosity_tp::SOME_MSG)
        printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, BOLDGREEN,
        	"No significant errors! (Abs diff is small)\n");

        goodDF = goodDF && true;
    }else{

    	if(rowMaxMax >= errScalar*pertSize){
    		// row error is too big
    		if(relDiff(rowMaxMaxIx, rowMaxIndex[rowMaxMaxIx]) < errScalar*pertSize){
    			// The relative difference is small enough
    			printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, BOLDGREEN,
    				"No significant errors! (Abs diff is large, rel diff is small)\n");
    			goodDF = goodDF && true;
    		}else{
    			printColor(BOLDRED, "Significant errors (Abs diff and Rel diff)\n");
    			goodDF = goodDF && false;
    		}
    	}else{
    		// Column error is too big
    		if(relDiff(colMaxIndex[colMaxMaxIx], colMaxMaxIx) < errScalar*pertSize){
    			printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, BOLDGREEN,
    				"No significant errors! (Abs diff is large, rel diff is small)\n");
    			goodDF = goodDF && true;
    		}else{
    			printColor(BOLDRED, "Significant errors (Abs diff and Rel diff)\n");
    			goodDF = goodDF && false;
    		}
    	}
    }

    // Print out big table of information
    if((!goodDF && verbosity >= Verbosity_tp::SOME_MSG) || verbosity >= Verbosity_tp::DEBUG){
        printf("Maximum relative difference between computed DF and estimated DF\n");
        int conCount = 0;
        for(long r = 0; r < rowMax.size(); r++){
            if(r == 0 && it.totalCons > 0){
                printf("Applies to %s %d: %s Constraint:\n", 
                    Constraint::getAppTypeStr(it.allCons[conCount].getAppType()),
                    it.allCons[conCount].getID(), it.allCons[conCount].getTypeStr());
            }else if(conCount+1 < static_cast<int>(it.allCons.size()) && r >= it.conRows[conCount+1]){
                conCount++;
                printf("Applies to %s %d: %s Constraint:\n", 
                    Constraint::getAppTypeStr(it.allCons[conCount].getAppType()),
                    it.allCons[conCount].getID(), it.allCons[conCount].getTypeStr());
            }
            astrohelion::printColor(rowMax[r] > errScalar*pertSize || std::isnan(rowMax[r]) ? RED : GREEN,
                "  row %03zu: %.6e (Abs) %.6e (Rel)\n", r, rowMax[r], relDiff(r, rowMaxIndex[r]));
        }
        for(long c = 0; c < colMax.size(); c++){
            std::string parent = "unknown";
            std::string type = "unknown";
            MSVarMap_Key key;
            try{
                key = freeVarMap_rowNumToKey.at(c);
                MSVarMap_Obj obj = it.freeVarMap.at(key);
                parent = MSVarMap_Obj::parent2str(obj.parent);
                type = MSVarMap_Key::type2str(key.type);
            }catch(std::out_of_range &e){}

            astrohelion::printColor(colMax[c] > errScalar*pertSize || std::isnan(colMax[c]) ? RED : GREEN,
                "Col %03zu: %s (%d)-owned %s: %.6e (Abs) %.6e (Rel)\n", c, parent.c_str(), key.id,
                type.c_str(), colMax[c], relDiff(colMaxIndex[c], c));
        }
    }

    return goodDF;
}//====================================================

/**
 *  @brief Compute the total delta-V along a corrected arcset
 * 
 *  @param it reference to an MultShootData object associated with a corrections process
 *  @return the total delta-V, units consistent with the arcset's stored velocity states
 */
double MultShootEngine::getTotalDV(const MultShootData& it){
    double total = 0;
    for(unsigned int n = 0; n < it.deltaVs.size()/3; n++){
        total += sqrt(it.deltaVs[3*n + 0]*it.deltaVs[3*n + 0] +
            it.deltaVs[3*n + 1]*it.deltaVs[3*n + 1] + 
            it.deltaVs[3*n + 2]*it.deltaVs[3*n + 2]);
    }
    return total;
}//=====================================================

}// END of Astrohelion namespace
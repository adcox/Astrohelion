/**
 *	\file MultShootEngine.cpp
 *	\brief Engine object that applies differential corrections to nodesets
 *	
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
 */

/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
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
#include <Eigen/OrderingMethods>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
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
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *  \brief Default constructor
 */
MultShootEngine::MultShootEngine(){}

/**
 *	\brief Copy constructor - create this engine by copying the input engine
 *	\param e input correction engine
 */
MultShootEngine::MultShootEngine(const MultShootEngine &e){
	copyMe(e);
}//=======================================================

/**
 *	\brief Destructor
 */
MultShootEngine::~MultShootEngine(){}

/**
 *	\brief Copy all engine variables
 *	\param e an engine reference
 */
void MultShootEngine::copyMe(const MultShootEngine &e){
	Engine::copyBaseEngine(e);
	verbosity = e.verbosity;//
	maxIts = e.maxIts;//
	tol = e.tol;//
	bFindEvent = e.bFindEvent;//
	bIgnoreCrash = e.bIgnoreCrash;//
	bIgnoreDiverge = e.bIgnoreDiverge;
	bFullFinalProp = e.bFullFinalProp;
	tofTp = e.tofTp;
	bLineSearchStepSize = e.bLineSearchStepSize;
}//====================================================

//-----------------------------------------------------
//      Operator Functions
//-----------------------------------------------------

/**
 *	\brief Copy operator; make a copy of the input correction engine.
 *
 *	\param e
 *	\return this correction engine
 */
MultShootEngine& MultShootEngine::operator =(const MultShootEngine &e){
	copyMe(e);
	return *this;
}//====================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *  \brief Retrieve whether or not the engine is locating an event crossing
 *	\return whether or not the algorithm will optimize the process to find an event
 */
bool MultShootEngine::isFindingEvent() const { return bFindEvent; }

/**
 *  \brief Retreive whether or not the engine will use a full, variable-step
 *  propagation for the final propagation.
 *  \details By default, this setting is TRUE. For lower computation time, 
 *  set to false via setFullFinalProp().
 *  
 *  \return whether or not the engine will use a full, variable-step
 *  propagation for the final propagation
 */
bool MultShootEngine::doesFullFinalProp() const { return bFullFinalProp; }

/**
 *  \brief Retreive whether or not the engine conducts a line search to choose
 *  the Newton step size.
 *  \return whether or not the engine conducts a line search to choose
 *  the Newton step size.
 */
bool MultShootEngine::doesLineSearch() const { return bLineSearchStepSize; }

/**
 *  \brief Retrieve the maximum number of iterations to attempt
 *	\return the maximum number of iterations to attempt before giving up
 */
int MultShootEngine::getMaxIts() const { return maxIts; }

/**
 *  \brief Retrieve the enumerated type describing how time-of-flight values
 *  are encoded (if at all) in the free variable vector
 *  
 *  \return the time-of-flight type
 */
MSTOF_tp MultShootEngine::getTOFType() const{ return tofTp; }

/**
 *  \brief Retrieve the minimum error tolerance
 *	\return the minimum error tolerance (non-dimensional units); errors
 *	less than this value are considered negligible
 */
double MultShootEngine::getTol() const { return tol; }

/**
 *  \brief Set whether or not the engine conducts a line search to choose 
 *  the Newton step size.
 *  \details Although the Newton step direction is guaranteed to point toward a 
 *  decreasing constraint vector, the full step may be too large. By leveraging
 *  a line search, the step size is chosen such that the magnitude of the 
 *  constraint vector decreases. WARNING: This adds many evaluations of expensive
 *  functions and will greatly decrease the speed. The line search is only recommended
 *  for particularly stubborn correction processes.
 * 
 *  \param b 
 */
void MultShootEngine::setDoLineSearch(bool b){ bLineSearchStepSize = b; }

/**
 *  \brief Set whether or not the engine will use a full, variable-step
 *  propagation for the final propagation
 *  \details By default, this setting is TRUE. For lower computation time, 
 *  set to false.
 * 
 *  \param b whether or not the engine will use a full, variable-step
 *  propagation for the final propagation
 */
void MultShootEngine::setFullFinalProp(bool b){ bFullFinalProp = b; }

/**
 * \brief Tell the corrector to ignore crash events (or to not to).
 * \details By default, the corrector does monitor crashes and will 
 * run into issues if the trajectory being corrected passes through a 
 * primary.
 * 
 * \param b whether or not to ignore crashes (default is false)
 */
void MultShootEngine::setIgnoreCrash(bool b){ bIgnoreCrash = b; }

/**
 *  \brief Tell the corrector to ignore divergence and return the partially
 *  corrected iteration data instead of throwing an exception when divergence
 *  occurs.
 * 
 *  \param b Whether or not to ignore divergance
 */
void MultShootEngine::setIgnoreDiverge(bool b){ bIgnoreDiverge = b;}

/**
 *	\brief Set maximum iterations
 *	\param i the maximum number of iterations to attempt before giving up
 */
void MultShootEngine::setMaxIts(int i){ maxIts = i; }

/**
 *  \brief Set the attenuation scalar and the limiting tolerance
 * 
 *  \param scale Scale the multiple shooting step by this value (multiply)
 *  \param limit Do not scale step if corrector error is below this value
 */
void MultShootEngine::setAttenuation(double scale, double limit){
	attenuation = scale;
	attenuationLimitTol = limit;
}//====================================================

/**
 *  \brief Set the way times-of-flight are encoded (if at all) in the
 *  free variable vector
 * 
 *  \param tp Describes how times-of-flight are encoded in the free
 *  variable vector
 */
void MultShootEngine::setTOFType(MSTOF_tp tp){ tofTp = tp; }

/**
 *	\brief Set the error tolerance
 *	\param d errors below this value will be considered negligible
 */
void MultShootEngine::setTol(double d){
	tol = d;

	if(tol > 1)
		astrohelion::printWarn("MultShootEngine::setTol: tolerance is greater than 1... just FYI\n");
}//====================================================

/**
 *	\brief Set the findEven flag
 *	\param b whether or not the algorithm will be looking for an event
 */
void MultShootEngine::setFindEvent(bool b){ bFindEvent = b; }

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	\brief Correct a generic nodeset using multiple shooting
 *	\details This algorithm employs multiple shooting to correct a set of nodes
 *	subject to a set of constraints. The nodes and constraints are all stored in the 
 *	input nodeset object. 
 *	
 *	Numerical integration is used to generate state transition matrices between
 *	nodes. In order to optimize performance and accuracy, only two steps are computed
 *	by the simulation engine, and the step size is fixed to force the usage of the 
 *	Adams-Bashforth Adams-Moulton method.
 *	
 *	\param set pointer to the nodeset that needs to be corrected
 *	\param pNodesOut pointer to the nodeset object that will contain the results of
 *	the shooting process
 *	\return the iteration data object for this corrections process
 *	\throws DivergeException if the corrections process does not converge
 *	\throws Exception
 *	* if the input and output nodesets contain different system data objects
 *	* if the dynamic model associated with the input
 *	nodeset does not support one or more of the nodeset's constraints
 *	* if the input nodeset contains more than one delta-v constraint
 *	* if the input nodeset contains more than one TOF constraint
 */
MultShootData MultShootEngine::multShoot(const Arcset *set, Arcset *pNodesOut){
	if(pNodesOut != nullptr && *(set->getSysData()) != *(pNodesOut->getSysData()))
		throw Exception("MultShootEngine::multShoot: Input and Output nodesets must use the same system data object");

	if(!bIsClean)
		cleanEngine();

	bIsClean = false;

	// Create structure to store iteration data for easy sharing
	MultShootData it(set);
	it.nodesOut = pNodesOut;
	it.tofTp = tofTp;
	// it.bVarTime = bVarTime;	// Save in structure to pass easily to other functions
	// it.bEqualArcTime = bEqualArcTime;
	
	astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "Multiple Shooting Algorithm:\n");
	astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  it.numNodes = %d\n", it.numNodes);
	astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  sysType = %s\n", set->getSysData()->getTypeStr().c_str());

	// Get the model associated with the nodeset
	const DynamicsModel *pModel = set->getSysData()->getDynamicsModel();
	pModel->multShoot_initDesignVec(&it);

	// Create constraints that enforce continuity between nodes; this process
	// does account for velocity discontinuities specified in the nodeset
	it.allCons.clear();
	pModel->multShoot_createContCons(&it);

	// Add all node constraints
	for(unsigned int n = 0; n < set->getNumNodes(); n++){
		std::vector<Constraint> nodeCons = set->getNodeRefByIx_const(n).getConstraints();
		it.allCons.insert(it.allCons.end(), nodeCons.begin(), nodeCons.end());
	}

	// Add all segment constraints
	for(unsigned int s = 0; s < set->getNumSegs(); s++){
		std::vector<Constraint> segCons = set->getSegRefByIx_const(s).getConstraints();
		it.allCons.insert(it.allCons.end(), segCons.begin(), segCons.end());
	}

	// Add all arcset constraints
	std::vector<Constraint> arcCons = set->getArcConstraints();
	it.allCons.insert(it.allCons.end(), arcCons.begin(), arcCons.end());

	// Compute number of extra consraint functions to add
	it.numSlack = 0;

	// Initialize vector to keep track of which row each constraint begins on
	// Also add slack variables to the design variable vector
	it.conRows.assign(it.allCons.size(), -1);	// Fill with -1
	int conRow = 0;
	bool foundDVCon = false;
	bool foundTOFCon = false;
	for(unsigned int c = 0; c < it.allCons.size(); c++){
		int addToRows = 0;
		Constraint con = it.allCons[c];

		if(!pModel->supportsCon(con.getType()))
			throw Exception("MultShootEngine::multShoot: The dynamic model does not support one of the constraints!");

		switch(con.getType()){
			case Constraint_tp::CONT_PV:
			case Constraint_tp::CONT_EX:
			case Constraint_tp::CONT_CTRL:
			case Constraint_tp::SEG_CONT_PV:
			case Constraint_tp::STATE:
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
				it.X.push_back(pModel->multShoot_getSlackVarVal(&it, con));
				it.slackAssignCon.push_back(c);
				it.numSlack++;
				break;
			case Constraint_tp::SP_MAX_DIST:
				it.X.push_back(pModel->multShoot_getSlackVarVal(&it, con));
				it.slackAssignCon.push_back(c);
				it.numSlack++;
			case Constraint_tp::SP_DIST:
				addToRows = 1;
				break;
			case Constraint_tp::MAX_DIST:
			case Constraint_tp::MIN_DIST:
			case Constraint_tp::ENDSEG_MAX_DIST:
			case Constraint_tp::ENDSEG_MIN_DIST:
				it.X.push_back(pModel->multShoot_getSlackVarVal(&it, con));
				it.slackAssignCon.push_back(c);	// remember where this slack variable is hiding
				it.numSlack++;
				// do NOT break here, continue on to do stuff for Constraint_tp::DIST as well
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
						/* Add a slack variable to the end of design vector and keep track
						 * of which constraint it is assigned to; value of slack
						 * variable will be recomputed later
						 */
						it.X.push_back(pModel->multShoot_getSlackVarVal(&it, con));
						it.numSlack++;
						it.slackAssignCon.push_back(c);
					}
				}else{
					throw Exception("MultShootEngine::multShoot: You can only apply ONE delta-V constraint");
				}
				break;
			case Constraint_tp::APSE:
			case Constraint_tp::ENDSEG_APSE:
			case Constraint_tp::ENDSEG_JC:
			case Constraint_tp::JC:
			case Constraint_tp::PSEUDOARC:
			case Constraint_tp::SEG_CONT_EX:
				addToRows = 1;
				break;
			case Constraint_tp::TOF:
				if(to_underlying(tofTp) <= 0)
					astrohelion::printWarn("MultShootEngine::multShoot: Attempting to constraint TOF without variable time... won't work!\n");
				
				if(!foundTOFCon)
					addToRows = 1;
				else
					throw Exception("MultShootEngine::multShoot: You can only apply ONE TOF constraint");
				break;
			case Constraint_tp::EPOCH:
				if(to_underlying(tofTp) <= 0)
					printWarn("MultShootEngine::multShoot: Attempting to constrain Epoch without variable time... won't work!\n");

				addToRows = 1;
				break;
			case Constraint_tp::RM_STATE:
			case Constraint_tp::RM_EPOCH:
			case Constraint_tp::RM_CTRL:
				// These constraints are handled differently
				addToRows = 0;
				break;
			default:
				printWarn("MultShootEngine::multShoot: Unhandled constraint: %s\n", con.getTypeStr());
				break;
		}

		// Save the index of the first row for this constraint
		it.conRows[c] = conRow;
		conRow += addToRows;	// remember we've added rows
		it.totalCons += addToRows;
	}// END of loop through constraints

	// Determine the number of free/design variables based on the system type
	it.totalFree = it.X.size();

	// Save the initial free variable vector
	it.X0 = it.X;

	// Print debugging information
	astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  # Free: %d\n  # Constraints: %d\n", it.totalFree, it.totalCons);
	astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  -> # Slack Variables: %d\n", it.numSlack);

	astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "ALL CONSTRAINTS:\n\n");
	if(verbosity >= Verbosity_tp::ALL_MSG){
		for(unsigned int n = 0; n < it.allCons.size(); n++){
			it.allCons[n].print();
		}
	}
	
	// Run the multiple shooting process
	return multShoot(it);
}//==========================================================

/**
 *  \brief Run a multiple shooting algorithm given an MultShootData object
 * 
 *  \param it A completely formed MultShootData object that describes a 
 *  multiple shooting problem. These are created from Arcset and its
 *  derivative types by the other implementation of multShoot()
 *  \return A corrected MultShootData object
 *  @see multShoot(Arcset*)
 *  \throws DivergeException if the multiple shooting process does not converge
 */
MultShootData MultShootEngine::multShoot(MultShootData it){
	it.count = 0;

	// create a simulation engine
	SimEngine simEngine;
	simEngine.setVerbosity(static_cast<Verbosity_tp>(static_cast<int>(verbosity) - 1));
	
	// Set both tolerances of simulation engine to be three orders of magnitude less corrector
	double simTol = tol/1000 < 1e-15 ? 1e-15 : tol/1000;
	simEngine.setAbsTol(simTol);
	simEngine.setRelTol(simTol);

	// Only need info about the final state, no need to generate lots of intermediate data points
	// This forces the integrator to use the Adams-Bashforth Adams-Moulton integration method,
	// which is most similar to Matlab's ode113
	simEngine.setVarStepSize(false);
	simEngine.setNumSteps(2);

	if(bFindEvent || bIgnoreCrash){
		simEngine.setMakeDefaultEvents(false);	// don't use crash events when searching for an event
	}

	// Define values for use in corrections loop
	double err = 10*tol;
	unsigned int coreStateSize = it.nodesIn->getSysData()->getDynamicsModel()->getCoreStateSize();

	Eigen::VectorXd oldX(it.totalFree, 1), newX(it.totalFree, 1), FX(it.totalCons, 1);

	while( err > tol && it.count < maxIts){
		if(it.count > 0){
			// Solve for newX and copy into working vector X
			oldX = Eigen::Map<Eigen::VectorXd>(&(it.X[0]), it.totalFree, 1);
			solveUpdateEq(&it, &oldX, &FX, &newX);

			it.X.clear();
			it.X.insert(it.X.begin(), newX.data(), newX.data()+it.totalFree);
		}

		it.FX.clear();					// Clear vectors each iteration
		it.FX.assign(it.totalCons, 0);	// Size the vectors and fill with zeros

		it.DF_elements.clear();
		it.DF_elements.reserve(it.totalCons * coreStateSize);

		// Fill each trajectory object with a propagated arc
		propSegsFromFreeVars(&it, &simEngine);
		// waitForUser();

		// Loop through all constraints and compute the constraint values, partials, and
		// apply them to the FX and DF matrices
		for(unsigned int c = 0; c < it.allCons.size(); c++){
			it.nodesIn->getSysData()->getDynamicsModel()->multShoot_applyConstraint(&it, it.allCons[c], c);
			printVerb(verbosity >= Verbosity_tp::DEBUG, "* Applying %s constraint\n", it.allCons[c].getTypeStr());
		}

		// Check to see what the error is; if it's too high, update X and continue another iteration
		FX = Eigen::Map<Eigen::VectorXd>(&(it.FX[0]), it.totalCons, 1);
		double err_cons = FX.norm();

		if(verbosity >= Verbosity_tp::ALL_MSG)
			reportConMags(&it);

		// Compute error: difference between subsequent free variable vectors
		// Eigen::VectorXd diff = newX - oldX;
		// double err_it = diff.norm();

		// Choose the lower of the two?
		// err = err_cons < err_it ? err_cons : err_it;
		// std::string errType = err_cons < err_it ? "||F||" : "||X - X_old||";
		err = err_cons;
		std::string errType = "||F||";

		it.count++;
		astrohelion::printVerbColor((bFindEvent && verbosity >= Verbosity_tp::ALL_MSG) || (!bFindEvent && verbosity > Verbosity_tp::NO_MSG), YELLOW, "Iteration %02d: err = %.8e (%s)\n",
			it.count, err, errType.c_str());
	}// end of corrections loop

	if(err > tol && !bIgnoreDiverge){
		throw DivergeException();
	}

	if(it.nodesOut){
		try{
			// Save propagated data and free variable vector values to the output arcset
			if(bFullFinalProp){
				simEngine.setVarStepSize(true);
				propSegsFromFreeVars(&it, &simEngine);
			}
			it.nodesIn->getSysData()->getDynamicsModel()->multShoot_createOutput(&it);
		}catch(Exception &e){
			astrohelion::printErr("MultShootEngine::multShoot: Unable to create output nodeset\n  Err: %s\n", e.what());
			throw e;
		}
	}
	
	return it;
}//=====================================================

/**
 *  \brief Propagate all segments along the arcset
 *  \details Initial conditions and other parameters for each integrated
 *  arc are obtained from the free variable vector or the input nodeset
 * 
 *  \param pIt Pointer to the multiple shooting data structure
 *  \param pSim Pointer to a simulalation engine initialized for multiple
 *  shooting propagations
 */
void MultShootEngine::propSegsFromFreeVars(MultShootData *pIt, SimEngine *pSim){
	unsigned int coreStateSize = pIt->nodesIn->getSysData()->getDynamicsModel()->getCoreStateSize();

	pIt->deltaVs.clear();
	pIt->deltaVs.assign(3*pIt->numNodes, 0);

	// initialize a vector of trajectory objects to store each propagated segment
	pIt->propSegs.clear();
	pIt->nodesIn->getSysData()->getDynamicsModel()->multShoot_initIterData(pIt);

	for(unsigned int s = 0; s < pIt->nodesIn->getNumSegs(); s++){
		// printf("Retrieving ICs for segment (ix %02d):\n", s);
		// Get simulation conditions from design vector via dynamic model implementation
		double t0 = 0, tof = 0;
		std::vector<double> ic(coreStateSize, 0);

		ControlLaw *pLaw = pIt->nodesIn->getSegRefByIx_const(s).getCtrlLaw();
		std::vector<double> ctrl0;
		if(pLaw){
			ctrl0.assign(pLaw->getNumStates(), 0);
		}

		pIt->nodesIn->getSysData()->getDynamicsModel()->multShoot_getSimICs(pIt, pIt->nodesIn->getSegRefByIx_const(s).getID(),
			&(ic.front()), &(ctrl0.front()), &t0, &tof);

		pSim->setRevTime(tof < 0);
		
		printVerb(verbosity >= Verbosity_tp::DEBUG, "Simulating segment %d:\n  t0 = %.4f\n  tof = %.4f\n", s, t0, tof);

		try{
			pSim->runSim(ic, ctrl0, t0, tof, &(pIt->propSegs[s]), pLaw);
		}catch(DivergeException &e){
			printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, RED, "SimEngine integration diverged...\n");
		}catch(Exception &e){
			printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, RED, "SimEngine Error:\n%s\nEnding corrections.\n", e.what());
		}

		// if(verbosity >= Verbosity_tp::DEBUG){
		// 	pIt->propSegs[s].print();
		// }

		std::vector<double> lastState = pIt->propSegs[s].getStateByIx(-1);
		int termID = pIt->nodesIn->getSegRefByIx_const(s).getTerminus();
		if(termID != Linkable::INVALID_ID){

			// Get the state of the terminal node
			MSVarMap_Obj stateVar = pIt->getVarMap_obj(MSVar_tp::STATE, termID);
			std::vector<double> state;
			if(stateVar.row0 == -1)
				state = pIt->nodesIn->getState(termID);
			else
				state = std::vector<double>(pIt->X.begin() + stateVar.row0, pIt->X.begin() + stateVar.row0 + stateVar.nRows);

			// velCon has false for a velocity state if there is a discontinuity between
			// the terminus of the segment and the terminal node
			std::vector<bool> velCon = pIt->nodesIn->getSegRefByIx_const(s).getVelCon();
			for(int i = 3; i < 6; i++){
				// Compute difference in velocity; if velCon[i-3] is true, then velocity
				// should be continuous and any difference is numerical error, so set to
				// zero by multiplying by not-true
				pIt->deltaVs[s*3+i-3] = !velCon[i-3]*(lastState[i] - state[i]);
			}
		}
	}
}//=====================================================

/**
 *	\brief Apply linear algebra to solve the update equation and obtain an updated free-variable vector
 *
 *	The update equation takes the form
 * 	\f[
 *		\vec{F}(\vec{X}) = D\vec{F}(\vec{X}) \left( \vec{X}_{n+1} - \vec{X}_n \right)
 *	\f]
 *	We know \f$ \vec{X}_n \f$, \f$ \vec{F}(\vec{X}) \f$ and \f$ D\vec{F}(\vec{X}) \f$, and we wish to solve for
 *	\f$ \vec{X}_{n+1} \f$. To do this, we need to invert or factor the Jacobian matrix \f$ D\vec{F}(\vec{X}) \f$.
 *	Assuming the Jacobian is non-singular, a solution can be obtained. Three possible scenarios can occur: either
 *	the system is over-constrained (no solution), perfectly constrained (square, one solution), or under constrained
 *	(infinitely many solutions).
 *
 *	In the first case, we can use least squares to find the closest possible solution. This is not implemented and
 *	the function will throw an error of the system is over constrained. In the second case, the equation is solved
 *	via LU factorization using GSL's linear algebra functions. Finally, if the system is under constrained, which 
 *	is the most common case, we compute the minimum-norm solution. We sequentially solve the following equations to
 *	arrive at the min-norm solution:
 *	\f{eqnarray*}{
 *		JJ^T \vec{w} &=& \vec{F}(\vec{X}) \\
 *		\vec{X}^* &=& J^T \vec{w}
 *	\f}
 *	where \f$ J \f$ is the Jacobian matrix, \f$ JJ^T \f$ is the associated Gram Matrix, and \f$ \vec{X}^* \f$ is the
 *	min-norm solution. If the Jacobian is non-singular, then the Gram Matrix will also be non-singular and the system
 *	can be solved. Note that \f$ \vec{X}^* \f$ is the min-norm solution for \f$ \vec{X}_{n+1} = \vec{X}_n \f$.
 *
 *	In all cases, errors will be thrown if the Jacobian is singular. This most likely indicates that there has been
 *	a coding error in the corrector, although singular Jacobians do occur when trajectories pass very near primaries.
 *
 *	\param pIt pointer to the MultShootData object associated with the corrections process
 *	\param pOldX constant pointer to the current (previous) design vector
 *	\param pFX constant pointer to the current constraint vector
 *	\param pNewX pointer to a vector in which to store the updated design variable vector
 *	
 *	\throws Exception if the problem is over constrained (i.e. Jacobian has more rows than columns);
 *	This can be updated to use a least-squares solution (TODO)
 */
void MultShootEngine::solveUpdateEq(MultShootData* pIt, const Eigen::VectorXd* pOldX, const Eigen::VectorXd *pFX, Eigen::VectorXd *pNewX){
	Eigen::VectorXd fullStep(pIt->totalFree, 1);	// Create a vector to put the solution in

	if(pIt->totalCons == 1 && pIt->totalFree == 1){
		// If Jacobian is 1x1, skip all that linear algebra and just solve the equation
		fullStep = -(*pFX)/(pIt->DF_elements[0].value());
	}else{
		// Jacobian is not scalar; use linear algebra to solve update equation
		SparseMatXCd J(pIt->totalCons, pIt->totalFree);
		J.setFromTriplets(pIt->DF_elements.begin(), pIt->DF_elements.end());
		J.makeCompressed();
		
		if(pIt->totalCons == pIt->totalFree){	// J is square, use regular inverse

			/* Use LU decomposition to invert the Jacobian matrix and find a vector
			w. Multiplying J^T by w yields the minimum-norm solution x, where x 
			lies in the column-space of J^T, or in the orthogonal complement of
			the nullspace of J.
			Source: <http://www.math.usm.edu/lambers/mat419/lecture15.pdf>
			 */
			
			// Info about solving Sparse matrix systems with Eigen:
			// <https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html>

			// Solve the system Jw = b (In this case, w = fullStep)
			Eigen::SparseLU<SparseMatXCd, Eigen::COLAMDOrdering<int> > luSolver;
			luSolver.analyzePattern(J);
			luSolver.factorize(J);
			if(luSolver.info() != Eigen::Success){
				checkDFSingularities(J);
				throw LinAlgException("MultShootEngine::solveUpdateEq: Could not factorize Jacobian matrix");
			}

			fullStep = luSolver.solve(-(*pFX));
			if(luSolver.info() != Eigen::Success){
				checkDFSingularities(J);
				throw LinAlgException("MultShootEngine::solveUpdateEq: Could not solve update equation");
			}
		}else{
			if(pIt->totalCons < pIt->totalFree){	// Under-constrained

				/* Use LU decomposition to invert the Gramm matrix and find a vector
				w. Multiplying J^T by w yields the minimum-norm solution x, where x 
				lies in the column-space of J^T, or in the orthogonal complement of
				the nullspace of J.
				Source: <http://www.math.usm.edu/lambers/mat419/lecture15.pdf>
				 */

				// Compute Gramm matrix
				SparseMatXCd JT = J.transpose();
				SparseMatXCd G = J*JT;		// G will always be symmetric

				// astrohelion::toCSV(J, "DF_cpp.csv");
				// astrohelion::toCSV(FX, "FX_cpp.csv");
				// astrohelion::toCSV(oldX, "X_cpp.csv");

				// Solve the system Gw = b
				// LU-Factorization is faster and works most of the time
				Eigen::SparseLU<SparseMatXCd, Eigen::COLAMDOrdering<int> > luSolver;
				luSolver.compute(G);	
				if(luSolver.info() != Eigen::Success){
					// If LU factorization fails, try QR factorization - it tends to be more robust
					printVerb(verbosity >= Verbosity_tp::ALL_MSG, "LU Factorization failed, trying QR factorization to solve update equation.\n");
					Eigen::SparseQR<SparseMatXCd, Eigen::COLAMDOrdering<int> > qrSolver;
					qrSolver.compute(G);
					if(qrSolver.info() != Eigen::Success){
						checkDFSingularities(J);
						throw LinAlgException("MultShootEngine::solveUpdateEq: Could not factorize Gramm matrix.");
					}

					Eigen::VectorXd w = qrSolver.solve(-(*pFX));
					if(qrSolver.info() != Eigen::Success){
						checkDFSingularities(J);
						throw LinAlgException("MultShootEngine::solveUpdateEq: Could not solve update equation (Gramm)");
					}

					fullStep = JT*w;	// Compute optimal step from w
				}else{

					Eigen::VectorXd w = luSolver.solve(-(*pFX));
					if(luSolver.info() != Eigen::Success){
						checkDFSingularities(J);
						throw LinAlgException("MultShootEngine::solveUpdateEq: Could not solve update equation (Gramm)");
					}
					fullStep = JT*w;	// Compute optimal step from w
				}

				// Alternative Method: SVD
				// NOTE: This takes approximately five times as much computation time
				// JacobiSVD<MatrixXd> svd(JT, Eigen::ComputeFullV | Eigen::ComputeFullU);
				// svd.setThreshold(tol/100.f);

				// MatrixXd singVals = svd.singularValues();
				// double svdTol = (1e-12)*singVals(0);	// First singular value is the biggest one

				// MatrixXd Z = MatrixXd::Zero(JT.rows(), JT.cols());
				// for(unsigned int r = 0; r < singVals.rows(); r++){
				// 	Z(r,r) = std::abs(singVals(r)) > svdTol ? 1.0/singVals(r) : 0;
				// }

				// fullStep = svd.matrixU() * Z * svd.matrixV().transpose() * FX;

			}else{	// Over-constrained
				throw LinAlgException("MultShootEngine::solveUpdateEq: System is over constrained... No solution implemented");
			}
		}
	}

	if(bLineSearchStepSize){
		chooseStep_LineSearch(pIt, pOldX, pFX, &fullStep, pNewX);
	}else{
		double scale = pFX->norm() < attenuationLimitTol ? 1.0 : attenuation;
		*pNewX = *pOldX + scale*fullStep;	// newX = oldX + fullStep
	}
}// End of solveUpdateEq() ============================

void MultShootEngine::chooseStep_LineSearch(MultShootData* pIt, const Eigen::VectorXd* pOldX, const Eigen::VectorXd *pOldFX,
	const Eigen::VectorXd *pFullStep, Eigen::VectorXd *pNewX){

	// Fixed parameters
	double maxStepSize = 100;
	double alpha = 1e-4;

	double f_old = 0.5*pOldFX->norm();

	// Initialize full step
	Eigen::VectorXd fullStep = *pFullStep;
	
	// Scale if attempted step is too big (i.e., if there is some unbounded thing going on)
	if(fullStep.norm() > maxStepSize){
		fullStep *= maxStepSize/fullStep.norm();
	}

	SparseMatXCd J(pIt->totalCons, pIt->totalFree);
	J.setFromTriplets(pIt->DF_elements.begin(), pIt->DF_elements.end());
	J.makeCompressed();

	Eigen::VectorXd slopeProd = J.transpose() * (*pOldFX);
	slopeProd = slopeProd.transpose() * (*pFullStep);
	assert(slopeProd.rows() == 1 && slopeProd.cols() == 1);
	double slope = slopeProd(0);

	if(slope >= 0.0){
		throw DivergeException("MultShootEngine::updateFreeVarVec: Slope is positive... roundoff error in line search.\n");
	}

	// Compute min step size
	double big = 0, temp = 0;
	for(unsigned int i = 0; i < fullStep.rows(); i++){
		temp = std::abs(fullStep(i))/std::max(std::abs((*pOldX)(i)), 1.0);
		if(temp > big)
			big = temp;
	}
	double minStep = 1e-7/big;
	double step = 1;

	printVerb(verbosity >= Verbosity_tp::ALL_MSG, "  Line Search: Minimum Step Size = %.4e\n", minStep);

	// Set up a simulation engine to propagate segments
	SimEngine sim;
	double simTol = tol/1000 < 1e-15 ? 1e-15 : tol/1000;
	sim.setAbsTol(simTol);
	sim.setRelTol(simTol);
	sim.setVarStepSize(false);
	sim.setNumSteps(2);
	sim.setMakeDefaultEvents(false);
	sim.setVerbosity(verbosity);
	
	unsigned int maxCount = 10;		// Max number of line search iterations
	double tempStep = 1;			// Storage for the next step size
	double step2 = 1;				// Step size from the previous iteration of the line search
	double f = 1;					// Error term for the current iteration
	double f2 = 1;;					// Error term from the previous iteration of the line search

	unsigned int coreDim = pIt->nodesIn->getSysData()->getDynamicsModel()->getCoreStateSize();

	for(unsigned int count = 0; count < maxCount; count++){
		// Compute next free variable vector give the current step size
		*pNewX = *pOldX + step * fullStep;

		// Compute an updated constraint vector
		MultShootData tempData = *pIt;
		tempData.FX.clear();
		tempData.FX.assign(tempData.totalCons, 0);
		tempData.DF_elements.clear();
		tempData.DF_elements.reserve(tempData.totalCons * coreDim);

		// Update free variable vector in the temporary data object
		assert(pNewX->rows() > 1 && pNewX->cols() == 1);
		tempData.X = std::vector<double>(pNewX->data(), pNewX->data() + pNewX->rows());
		
		// Do all those expensive function evaluations
		propSegsFromFreeVars(&tempData, &sim);
		for(unsigned int c = 0; c < tempData.allCons.size(); c++){
			tempData.nodesIn->getSysData()->getDynamicsModel()->multShoot_applyConstraint(&tempData, tempData.allCons[c], c);
		}

		// Compute magnitude of new FX vector
		temp = 0;
		for(unsigned int i = 0; i < tempData.FX.size(); i++){
			temp += tempData.FX[i]*tempData.FX[i];
		}
		// Scaled norm of the constraint vector (an error metric)
		f = 0.5*sqrt(temp);

		if(step < minStep){
			printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, RED, "  Line search decreased past minimum bound.\n");
			return;
		}


		if(f <= f_old + alpha*step*slope){
			printVerb(verbosity >= Verbosity_tp::SOME_MSG, "  Line Search: Step Size = %.4e (%u its)\n", step, count);
			return;	// We're all set, so return!
		}else{	// Need to backtrack
			if(count == 0){
				// Approximate step-size as a quadratic, find step size that yields better guess
				tempStep = -slope/(2.0*(f - f_old - slope));
			}else{
				// Approximate step-size as a cubic
				double term1 = f - f_old - step*slope;
				double term2 = f2 - f_old - step2*slope;
				double a = (term1/(step*step) - term2/(step2*step2))/(step - step2);
				double b = (-step2*term1/(step*step) + step*term2/(step2*step2))/(step - step2);

				if(a == 0){
					tempStep = -slope/(2.0*b);
				}else{
					double disc = b*b - 3*a*slope;
					if(disc < 0){
						tempStep = 0.5*step;
					}else if(b <= 0){
						tempStep = (-b + sqrt(disc))/(3.0*a);
					}else{
						tempStep = -slope/(b + sqrt(disc));
					}

					// Ensure that step size decreases by at least a factor of 2
					if(tempStep > 0.5*step){
						tempStep = 0.5*step;
					}
				}
			}

			// Update storage for next iteration
			step2 = step;
			f2 = f;

			// Ensure that step size decreases by no more than 90%
			step = std::max(tempStep, 0.1*step);
		}
	}

	if(f <= f_old + alpha*step*slope)
		throw DivergeException("MultShootEngine::updateFreeVarVec: Line search for step size diverged.");

	printVerb(verbosity >= Verbosity_tp::SOME_MSG, "  Line Search: Step Size = %.4e (%d its)\n", step, maxCount);
}//====================================================

/**
 *  \brief Checks the Jacobian (DF) matrix for singularities, i.e., rows
 *  or columns that contain only zeros
 * 
 *  \param DF Jacobian matrix employed in the correction process
 */
void MultShootEngine::checkDFSingularities(MatrixXRd DF){
	if(verbosity >= Verbosity_tp::SOME_MSG){
		for(unsigned int r = 0; r < DF.rows(); r++){
			if(DF.row(r).norm() == 0){
				printErr("Singular Jacobian! Row %u is all zeros\n", r);
			}
		}

		for(unsigned int c = 0; c < DF.cols(); c++){
			if(DF.col(c).norm() == 0){
				printErr("Singular Jacobian! Column %u is all zeros\n", c);
			}
		}
	}
}//====================================================

/**
 *  \brief Print out the magnitude of each constraint.
 *  \details This can be useful when debugging to highlight which constraints are unsatisfied
 * 
 *  \param pIt pointer to an MultShootData object associated with a corrections process
 */
void MultShootEngine::reportConMags(const MultShootData *pIt){
	unsigned int conCount = 0;
	for(unsigned int r = 0; r < (pIt->FX.size()); r++){
        if(r == 0 && pIt->totalCons > 0){
            printf("Node %d %s Constraint:\n", pIt->allCons[conCount].getID(), pIt->allCons[conCount].getTypeStr());
        }else if(conCount < pIt->allCons.size() && r >= static_cast<unsigned int>(pIt->conRows[conCount+1])){
            conCount++;
            printf("Node %d %s Constraint:\n", pIt->allCons[conCount].getID(), pIt->allCons[conCount].getTypeStr());
        }
        printf("  ||row %03u||: %.6e\n", r, std::abs(pIt->FX[r]));
    }
}//===============================================================

/**
 *	\brief clean up data so that engine can be used again (or deconstructed safely)
 */
void MultShootEngine::cleanEngine(){
	astrohelion::printVerb(verbosity >= Verbosity_tp::ALL_MSG, "Cleaning the engine...\n");
	bIsClean = true;
}//====================================================

/**
 *  \brief Reset the multiple shooting engine to the default parameter values
 */
void MultShootEngine::reset(){
	if(!bIsClean)
		cleanEngine();

	// bVarTime = true;
	// bEqualArcTime = false;
	tofTp = MSTOF_tp::VAR_FREE;
	maxIts = 20;
	tol = 1e-12;
	bFindEvent = false;
	bIgnoreCrash = false;
	bIgnoreDiverge = false;
	bFullFinalProp = true;
	attenuation = 1;
	attenuationLimitTol = 1e-8;
}//====================================================

/**
 *  \brief Check the DF matrix for the multiple shooting algorithm using finite differencing
 *  \details This function checks to make sure the Jacobian matrix (i.e. DF) is correct
 *  by computing each partial derivative numerically via forward differencing.
 * 
 *  \param pNodeset A nodeset with some constraints
 *  \param verbosity Specify how verbose the output is
 *  \param writeToFile Whether or not to write .csv or .mat files with relevant information
 */
bool MultShootEngine::finiteDiff_checkMultShoot(const Arcset *pNodeset, Verbosity_tp verbosity, bool writeToFile){
    MultShootEngine engine;  // Create engine with default settings
    return finiteDiff_checkMultShoot(pNodeset, engine, verbosity, writeToFile);
}//====================================================

/**
 *  \brief Check the DF matrix for the multiple shooting algorithm using finite differencing
 *  \details This function checks to make sure the Jacobian matrix (i.e. DF) is correct
 *  by computing each partial derivative numerically via forward differencing.
 * 
 *  \param pNodeset A nodeset with some constraints
 *  \param engine correction engine object configured with the appropriate settings (i.e.,
 *  equal arc time, etc.). Note that the maxIts, verbosity, and ignoreDiverge
 *  attributes of the engine will be overridden by this function.
 *  \param verbosity Specify how verbose the output is
 *  \param writeToFile Whether or not to write .csv or .mat files with relevant information
 */
bool MultShootEngine::finiteDiff_checkMultShoot(const Arcset *pNodeset, MultShootEngine engine, Verbosity_tp verbosity, bool writeToFile){
    printVerb(verbosity >= Verbosity_tp::SOME_MSG, "Finite Diff: Checking DF matrix...\n");
    // Create multiple shooter that will only do 1 iteration
    MultShootEngine corrector(engine);
    corrector.setMaxIts(1);
    corrector.setVerbosity(verbosity < Verbosity_tp::DEBUG ? Verbosity_tp::NO_MSG : Verbosity_tp::ALL_MSG);
    corrector.setIgnoreDiverge(true);
    corrector.setFullFinalProp(false);

    // Run multiple shooter to get X, FX, and DF
    MultShootData it = corrector.multShoot(pNodeset, nullptr);
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
        MultShootData pertIt = corrector.multShoot(it);     // Correct perturbed state
        Eigen::VectorXd FX_up = Eigen::Map<Eigen::VectorXd>(&(pertIt.FX[0]), it.totalCons, 1);

        // Do another process for opposite direction
        pertX = it.X0;
        pertX[i] -= pertSize;
        it.X = pertX;
        pertIt = corrector.multShoot(it);
        Eigen::VectorXd FX_down = Eigen::Map<Eigen::VectorXd>(&(pertIt.FX[0]), it.totalCons, 1);

        // An iteration for twice the perturbation up
        pertX = it.X0;
        pertX[i] += 2*pertSize;
        it.X = pertX;
        pertIt = corrector.multShoot(it);
        Eigen::VectorXd FX_2up = Eigen::Map<Eigen::VectorXd>(&(pertIt.FX[0]), it.totalCons, 1);

        // An iteration for twice the perturbation down
        pertX = it.X0;
        pertX[i] -= 2*pertSize;
        it.X = pertX;
        pertIt = corrector.multShoot(it);
        Eigen::VectorXd FX_2down = Eigen::Map<Eigen::VectorXd>(&(pertIt.FX[0]), it.totalCons, 1);


        // Compute central difference
        Eigen::VectorXd col = (-1*FX_2up + 8*FX_up - 8*FX_down + FX_2down)/std::abs(12*pertSize);   // Five-point stencil
        // Eigen::VectorXd col = (FX_up - FX_down)/std::abs(2*pertSize);   // Central Difference
        // Eigen::VectorXd col = (FX_up - FX)/std::abs(pertSize);       // Forward difference
        DFest.block(0, i, it.totalCons, 1) = col;
    }


    MatrixXRd diff = DF - DFest;
    MatrixXRd DF_abs = DF.cwiseAbs();       // Get coefficient-wise absolute value
    MatrixXRd DFest_abs = DFest.cwiseAbs();

    diff = diff.cwiseAbs();                     	// Get coefficient-wise aboslute value
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
        astrohelion::toCSV(DFest, "FiniteDiff_DFest.csv");
        astrohelion::toCSV(diff, "FiniteDiff_Diff.csv"); 
        astrohelion::toCSV(relDiff, "FiniteDiff_RelDiff.csv");
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
 *  \brief Compute the total delta-V along a corrected nodeset
 * 
 *  \param pIt pointer to an MultShootData object associated with a corrections process
 *  \return the total delta-V, units consistent with the nodeset's stored velocity states
 */
double MultShootEngine::getTotalDV(const MultShootData *pIt){
    double total = 0;
    for(unsigned int n = 0; n < pIt->deltaVs.size()/3; n++){
        total += sqrt(pIt->deltaVs[3*n + 0]*pIt->deltaVs[3*n + 0] +
            pIt->deltaVs[3*n + 1]*pIt->deltaVs[3*n + 1] + 
            pIt->deltaVs[3*n + 2]*pIt->deltaVs[3*n + 2]);
    }
    return total;
}//=====================================================

}// END of Astrohelion namespace
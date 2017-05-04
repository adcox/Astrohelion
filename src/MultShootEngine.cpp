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

#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <vector>

using Eigen::MatrixXd;
using Eigen::JacobiSVD;

namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

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
	bVarTime = e.bVarTime;//
	bEqualArcTime = e.bEqualArcTime;//
	maxIts = e.maxIts;//
	tol = e.tol;//
	bFindEvent = e.bFindEvent;//
	bIgnoreCrash = e.bIgnoreCrash;//
	bIgnoreDiverge = e.bIgnoreDiverge;
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
 *  \brief Retrieve whether or not we are using variable time
 *	\return whether or not the corrector uses variable time (as opposed
 * 	to fixed time)
 */
bool MultShootEngine::usesVarTime() const { return bVarTime; }

/**
 *	\brief Retrieve whether or not we force all segments to have the same length
 *	(in time).
 *
 *	This setting only applies if variable time is turned on.
 *	\return whether or not each arc will be forced to have the same length in time
 */
bool MultShootEngine::usesEqualArcTime() const { return bEqualArcTime; }

/**
 *  \brief Retrieve whether or not we are located an event crossing
 *	\return whether or not the algorithm will optimize the process to find an event
 */
bool MultShootEngine::isFindingEvent() const { return bFindEvent; }

/**
 *  \brief Retrieve the maximum number of iterations to attempt
 *	\return the maximum number of iterations to attempt before giving up
 */
int MultShootEngine::getMaxIts() const { return maxIts; }

/**
 *  \brief Retrieve the minimum error tolerance
 *	\return the minimum error tolerance (non-dimensional units); errors
 *	less than this value are considered negligible
 */
double MultShootEngine::getTol() const { return tol; }

/**
 *	\brief Set bVarTime
 *	\param b whether or not the corrector should use variable time
 */
void MultShootEngine::setVarTime(bool b){
	bVarTime = b;
	// Turn off equal-time arcs too if bVarTime is false
	if(!bVarTime)
		bEqualArcTime = false;
}//==================================================

/**
 *	\brief Tell the corrector how to apply variable time
 *	\param b whether or not each arc will be forced to have the same duration
 */
void MultShootEngine::setEqualArcTime(bool b){
	if(!bVarTime && b){
		astrohelion::printErr("MultShootEngine::setequalArcTime: Cannot use equal-time arcs if variable time is disabled; please turn bVarTime ON first\n");
		bEqualArcTime = false;
	}else{
		bEqualArcTime = b;
	}
}//==================================================

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
	it.bVarTime = bVarTime;	// Save in structure to pass easily to other functions
	it.bEqualArcTime = bEqualArcTime;
	
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
		std::vector<Constraint> nodeCons = set->getNodeByIx(n).getConstraints();
		it.allCons.insert(it.allCons.end(), nodeCons.begin(), nodeCons.end());
	}

	// Add all segment constraints
	for(unsigned int s = 0; s < set->getNumSegs(); s++){
		std::vector<Constraint> segCons = set->getSegByIx(s).getConstraints();
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
			case Constraint_tp::SEG_CONT_PV:
			case Constraint_tp::SEG_CONT_EX:
			case Constraint_tp::STATE:
			case Constraint_tp::MATCH_CUST:
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
				it.X.push_back(pModel->multShoot_getSlackVarVal(&it, con));
				it.slackAssignCon.push_back(c);	// remember where this slack variable is hiding
				it.numSlack++;
				// do NOT break here, continue on to do stuff for Constraint_tp::DIST as well
			case Constraint_tp::DIST:
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
			case Constraint_tp::JC:
			case Constraint_tp::PSEUDOARC:
				addToRows = 1;
				break;
			case Constraint_tp::TOF:
				if(!bVarTime)
					astrohelion::printWarn("MultShootEngine::multShoot: Attempting to constraint TOF without variable time... won't work!\n");
				
				if(!foundTOFCon)
					addToRows = 1;
				else
					throw Exception("MultShootEngine::multShoot: You can only apply ONE TOF constraint");
				break;
			case Constraint_tp::EPOCH:
				if(!bVarTime)
					printWarn("MultShootEngine::multShoot: Attempting to constrain Epoch without variable time... won't work!\n");

				addToRows = 1;
				break;
			default: break;
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

	Arcset arc(it.nodesIn->getSysData());
	while( err > tol && it.count < maxIts){
		it.FX.clear();					// Clear vectors each iteration
		it.DF.clear();
		it.deltaVs.clear();
		it.propSegs.clear();
		it.FX.assign(it.totalCons, 0);	// Size the vectors and fill with zeros
		it.DF.assign(it.totalCons*it.totalFree, 0);
		it.deltaVs.assign(3*it.numNodes, 0);

		// initialize a vector of trajectory objects to store each propagated segment
		it.nodesIn->getSysData()->getDynamicsModel()->multShoot_initIterData(&it);

		for(unsigned int s = 0; s < it.nodesIn->getNumSegs(); s++){
			// printf("Retrieving ICs for segment (ix %02d):\n", s);
			// Get simulation conditions from design vector via dynamic model implementation
			double t0 = 0, tof = 0;
			std::vector<double> ic(coreStateSize, 0);
			it.nodesIn->getSysData()->getDynamicsModel()->multShoot_getSimICs(&it, it.nodesIn->getSegByIx(s).getID(),
				&(ic[0]), &t0, &tof);

			simEngine.setRevTime(tof < 0);
			
			printVerb(verbosity >= Verbosity_tp::DEBUG, "Simulating segment %d:\n  t0 = %.4f\n  tof = %.4f\n", s, t0, tof);
			simEngine.setCtrlLaw(it.nodesIn->getSegByIx(s).getCtrlLaw());

			try{
				simEngine.runSim(ic, t0, tof, &(it.propSegs[s]));
			}catch(DivergeException &e){
				printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, RED, "SimEngine integration diverged...\n");
			}catch(Exception &e){
				printVerbColor(verbosity >= Verbosity_tp::SOME_MSG, RED, "SimEngine Error:\n%s\nEnding corrections.\n", e.what());
			}

			// if(verbosity >= Verbosity_tp::DEBUG){
			// 	it.propSegs[s].print();
			// }
		}

		// waitForUser();
		
		// Compute Delta-Vs between node segments
		for(unsigned int s = 0; s < it.nodesIn->getNumSegs(); s++){
			std::vector<double> lastState = it.propSegs[s].getStateByIx(-1);
			int termID = it.nodesIn->getSegByIx(s).getTerminus();
			if(termID != Linkable::INVALID_ID){
				int termNodeIx = it.nodesIn->getNodeIx(termID);
				// velCon has false for a velocity state if there is a discontinuity between
				// the terminus of the segment and the terminal node
				std::vector<bool> velCon = it.nodesIn->getSegByIx(s).getVelCon();
				for(int i = 3; i < 6; i++){
					// Compute difference in velocity; if velCon[i-3] is true, then velocity
					// should be continuous and any difference is numerical error, so set to
					// zero by multiplying by not-true
					it.deltaVs[s*3+i-3] = !velCon[i-3]*(lastState[i] - it.X[6*termNodeIx+i]);
				}
			}
		}

		// Loop through all constraints and compute the constraint values, partials, and
		// apply them to the FX and DF matrices
		for(unsigned int c = 0; c < it.allCons.size(); c++)
			it.nodesIn->getSysData()->getDynamicsModel()->multShoot_applyConstraint(&it, it.allCons[c], c);

		// Solve for newX and copy into working vector X
		Eigen::VectorXd oldX = Eigen::Map<Eigen::VectorXd>(&(it.X[0]), it.totalFree, 1);
		Eigen::VectorXd newX = solveUpdateEq(&it);

		it.X.clear();
		it.X.insert(it.X.begin(), newX.data(), newX.data()+it.totalFree);

		// Compute error; norm of constraint vector
		Eigen::VectorXd FX = Eigen::Map<Eigen::VectorXd>(&(it.FX[0]), it.totalCons, 1);
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
			it.nodesIn->getSysData()->getDynamicsModel()->multShoot_createOutput(&it);
		}catch(Exception &e){
			astrohelion::printErr("MultShootEngine::multShoot: Unable to create output nodeset\n  Err: %s\n", e.what());
			throw e;
		}
	}
	
	return it;
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
 *
 *	\return the updated free variable vector \f$ \vec{X}_{n+1} \f$
 *	\throws Exception if the problem is over constrained (i.e. Jacobian has more rows than columns);
 *	This can be updated to use a least-squares solution (TODO)
 */
Eigen::VectorXd MultShootEngine::solveUpdateEq(MultShootData* pIt){
	// Create matrices for X, Jacobian matrix DF, and constraint vector FX
	Eigen::VectorXd oldX = Eigen::Map<Eigen::VectorXd>(&(pIt->X[0]), pIt->totalFree, 1);
	MatrixXRd J = Eigen::Map<MatrixXRd>(&(pIt->DF[0]), pIt->totalCons, pIt->totalFree);
	Eigen::VectorXd FX = Eigen::Map<Eigen::VectorXd>(&(pIt->FX[0]), pIt->totalCons, 1);

	// change sign for matrix multiplication
	FX *= -1;

	// Create a vector to put the solution in
	Eigen::VectorXd X_diff(pIt->totalFree, 1);
	if(pIt->totalCons == pIt->totalFree){	// J is square, use regular inverse

		/* Use LU decomposition to invert the Jacobian matrix and find a vector
		w. Multiplying J^T by w yields the minimum-norm solution x, where x 
		lies in the column-space of J^T, or in the orthogonal complement of
		the nullspace of J.
		Source: <http://www.math.usm.edu/lambers/mat419/lecture15.pdf>
		 */
		
		// Solve the system Jw = b (In this case, w = X_diff)
		Eigen::FullPivLU<MatrixXRd> lu(J);
		// lu.setThreshold(1e-20);

		// if(!lu.isInvertible())
		// 	throw LinAlgException("MultShootEngine::solveUpdateEq: Jacobian is singular; cannot solve");

		X_diff = lu.solve(FX);
	}else{
		if(pIt->totalCons < pIt->totalFree){	// Under-constrained
			// Compute Gramm matrix
			MatrixXRd JT = J.transpose();
			MatrixXRd G = J*JT;

			// astrohelion::toCSV(J, "DF_cpp.csv");
			// astrohelion::toCSV(FX, "FX_cpp.csv");
			// astrohelion::toCSV(oldX, "X_cpp.csv");

			// waitForUser();
			/* Use LU decomposition to invert the Gramm matrix and find a vector
			w. Multiplying J^T by w yields the minimum-norm solution x, where x 
			lies in the column-space of J^T, or in the orthogonal complement of
			the nullspace of J.
			Source: <http://www.math.usm.edu/lambers/mat419/lecture15.pdf>
			 */
			
			// Solve the system Gw = b
			Eigen::FullPivLU<MatrixXRd> lu(G);
			// lu.setThreshold(1e-20);

			// if(!lu.isInvertible()){
			// 	astrohelion::printErr("Gramm Matrix rank = %ld / %ld\n", lu.rank(), G.rows());
			// 	throw LinAlgException("MultShootEngine::solveUpdateEq: Gramm matrix is singular; cannot solve");
			// }

			Eigen::VectorXd w = lu.solve(FX);
			
			// Compute optimal x from w
			X_diff = JT*w;



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

			// X_diff = svd.matrixU() * Z * svd.matrixV().transpose() * FX;

		}else{	// Over-constrained
			throw LinAlgException("System is over constrained... No solution implemented");
		}
	}

	double scale = FX.norm() < attenuationLimitTol ? 1.0 : attenuation;
	return oldX + scale*X_diff;	// newX = oldX + X_diff
}// End of solveUpdateEq() =====================================

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

void MultShootEngine::reset(){
	if(!bIsClean)
		cleanEngine();

	bVarTime = true;
	bEqualArcTime = false;
	maxIts = 20;
	tol = 1e-12;
	bFindEvent = false;
	bIgnoreCrash = false;
	bIgnoreDiverge = false;
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
    printVerb(verbosity >= Verbosity_tp::SOME_MSG, "Finite Diff: Checking DF matrix... ");
    // Create multiple shooter that will only do 1 iteration
    MultShootEngine corrector(engine);
    corrector.setMaxIts(1);
    corrector.setVerbosity(Verbosity_tp::NO_MSG);
    corrector.setIgnoreDiverge(true);

    // Run multiple shooter to get X, FX, and DF
    MultShootData it = corrector.multShoot(pNodeset, nullptr);
    Eigen::VectorXd FX = Eigen::Map<Eigen::VectorXd>(&(it.FX[0]), it.totalCons, 1);
    MatrixXRd DF = Eigen::Map<MatrixXRd>(&(it.DF[0]), it.totalCons, it.totalFree);
    MatrixXRd DFest = MatrixXRd::Zero(it.totalCons, it.totalFree);

    double pertSize = 1e-8;
    #pragma omp parallel for firstprivate(it, corrector)
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

    diff = diff.cwiseAbs();                     // Get coefficient-wise aboslute value

    // // Divide each element by the magnitude of the DF element to get a relative difference magnitude
    // for(int r = 0; r < diff.rows(); r++){
    //     for(int c = 0; c < diff.cols(); c++){
    //         // If one of the elements is zero, let the difference just be the difference; no division
    //         if(DF_abs(r,c) > 1e-13 && DFest_abs(r,c) > 1e-13)   // consider 1e-13 to be zero
    //             diff(r,c) = diff(r,c)/DF_abs(r,c);

    //         // if(r == 50 && c == 98){
    //         //     printf("DF(50, 98) = %.4e\n", DF(r,c));
    //         //     printf("DFest(50, 98) = %.4e\n", DFest(r,c));
    //         // }
    //     }
    // }
    
    if(writeToFile){
        astrohelion::toCSV(DF, "FiniteDiff_DF.csv");
        astrohelion::toCSV(DFest, "FiniteDiff_DFest.csv");
        astrohelion::toCSV(diff, "FiniteDiff_Diff.csv"); 
    }
    // 

    Eigen::VectorXd rowMax = diff.rowwise().maxCoeff();
    Eigen::RowVectorXd colMax = diff.colwise().maxCoeff();

    // Make a map that links the row number (column in the Jacobian matrix)
    // of the free variable to the MSVarMap_Key that represents the MSVarMap_Obj 
    // that contains information about the free variable
    std::map<int, MSVarMap_Key> freeVarMap_rowNumToKey;
    for(auto& obj : it.freeVarMap){
        for(int r = 0; r < obj.second.nRows; r++){
            freeVarMap_rowNumToKey[obj.second.row0 + r] = obj.first;
        }
    }

    double rowMaxMax = rowMax.maxCoeff();
    double colMaxMax = colMax.maxCoeff();
    int errScalar = 10000;

    if(rowMaxMax < errScalar*pertSize && colMaxMax < errScalar*colMaxMax){
        if(verbosity >= Verbosity_tp::SOME_MSG)
            astrohelion::printColor(BOLDGREEN, "No significant errors!\n");

        return true;
    }else{
        if(verbosity >= Verbosity_tp::SOME_MSG){
            astrohelion::printColor(BOLDRED, "Significant errors!\n");
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
                    "  row %03zu: %.6e\n", r, rowMax[r]);
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
                    "Col %03zu: %s (%d)-owned %s: %.6e\n", c, parent.c_str(), key.id, type.c_str(), colMax[c]);
            }
        }

        return false;
    }
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
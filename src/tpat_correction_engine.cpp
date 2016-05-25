/**
 *	@file tpat_correction_engine.cpp
 *	@brief Engine object that applies differential corrections to nodesets
 */

/*
 *	Trajectory Propagation and Analysis Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
 *
 *  TPAT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  TPAT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with TPAT.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#include "tpat_correction_engine.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_calculations.hpp"
#include "tpat_exceptions.hpp" 
#include "tpat_multShoot_data.hpp"
#include "tpat_node.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sim_engine.hpp"
#include "tpat_utilities.hpp"

#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <vector>

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Copy constructor - create this engine by copying the input engine
 *	@param e input correction engine
 */
TPAT_Correction_Engine::TPAT_Correction_Engine(const TPAT_Correction_Engine &e){
	copyEngine(e);
}//=======================================================

/**
 *	@brief Destructor
 */
TPAT_Correction_Engine::~TPAT_Correction_Engine(){}

/**
 *	@brief Copy all engine variables
 *	@param e an engine reference
 */
void TPAT_Correction_Engine::copyEngine(const TPAT_Correction_Engine &e){
	verbose = e.verbose;//
	varTime = e.varTime;//
	equalArcTime = e.equalArcTime;//
	maxIts = e.maxIts;//
	tol = e.tol;//
	findEvent = e.findEvent;//
	ignoreCrash = e.ignoreCrash;//
	ignoreDiverge = e.ignoreDiverge;
	scaleVars = e.scaleVars;
	isClean = e.isClean;
}//====================================================

//-----------------------------------------------------
//      Operator Functions
//-----------------------------------------------------

/**
 *	@brief Copy operator; make a copy of the input correction engine.
 *
 *	@param e
 *	@return this correction engine
 */
TPAT_Correction_Engine& TPAT_Correction_Engine::operator =(const TPAT_Correction_Engine &e){
	copyEngine(e);
	return *this;
}//====================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *  @brief Retrieve whether or not we are using variable time
 *	@return whether or not the corrector uses variable time (as opposed
 * 	to fixed time)
 */
bool TPAT_Correction_Engine::usesVarTime() const { return varTime; }

/**
 *	@brief Retrieve whether or not we force all segments to have the same length
 *	(in time).
 *
 *	This setting only applies if variable time is turned on.
 *	@return whether or not each arc will be forced to have the same length in time
 */
bool TPAT_Correction_Engine::usesEqualArcTime() const { return equalArcTime; }

/**
 *  @brief Retrieve whether or not the multiple shooting algorithm uses variable scaling
 *  @return whether or not the multiple shooting algorithm uses variable scaling
 */
bool TPAT_Correction_Engine::usesScaledVars() const { return scaleVars; }

/**
 *  @brief Retrieve the verbosity setting
 *	@return whether or not the corrector will be verbose
 */
TPAT_Verbosity_Tp TPAT_Correction_Engine::isVerbose() const { return verbose; }

/**
 *  @brief Retrieve whether or not we are located an event crossing
 *	@return whether or not the algorithm will optimize the process to find an event
 */
bool TPAT_Correction_Engine::isFindingEvent() const { return findEvent; }

/**
 *  @brief Retrieve the maximum number of iterations to attempt
 *	@return the maximum number of iterations to attempt before giving up
 */
int TPAT_Correction_Engine::getMaxIts() const { return maxIts; }

/**
 *  @brief Retrieve the minimum error tolerance
 *	@return the minimum error tolerance (non-dimensional units); errors
 *	less than this value are considered negligible
 */
double TPAT_Correction_Engine::getTol() const { return tol; }

/**
 *	@brief Set varTime
 *	@param b whether or not the corrector should use variable time
 */
void TPAT_Correction_Engine::setVarTime(bool b){
	varTime = b;
	// Turn off equal-time arcs too if varTime is false
	if(!varTime)
		equalArcTime = false;
}//==================================================

/**
 *	@brief Tell the corrector how to apply variable time
 *	@param b whether or not each arc will be forced to have the same duration
 */
void TPAT_Correction_Engine::setEqualArcTime(bool b){
	if(!varTime && b){
		printErr("TPAT_Correction_Engine::setequalArcTime: Cannot use equal-time arcs if variable time is disabled; please turn varTime ON first\n");
		equalArcTime = false;
	}else{
		equalArcTime = b;
	}
}//==================================================

/**
 * @brief Tell the corrector to ignore crash events (or to not to).
 * @details By default, the corrector does monitor crashes and will 
 * run into issues if the trajectory being corrected passes through a 
 * primary.
 * 
 * @param b whether or not to ignore crashes (default is false)
 */
void TPAT_Correction_Engine::setIgnoreCrash(bool b){ ignoreCrash = b; }

/**
 *  @brief Tell the corrector to ignore divergence and return the partially
 *  corrected iteration data instead of throwing an exception when divergence
 *  occurs.
 * 
 *  @param b Whether or not to ignore divergance
 */
void TPAT_Correction_Engine::setIgnoreDiverge(bool b){ ignoreDiverge = b;}

/**
 *	@brief Set verbosity
 *	@param b whether or not the corrector should be verbose in its outputs
 */
void TPAT_Correction_Engine::setVerbose(TPAT_Verbosity_Tp b){ verbose = b; }

/**
 *	@brief Set maximum iterations
 *	@param i the maximum number of iterations to attempt before giving up
 */
void TPAT_Correction_Engine::setMaxIts(int i){ maxIts = i; }

/**
 *  @brief Set the scaleVar flag
 * 
 *  @param b whether or not the multiple shooting algorithm should use variable scaling
 */
void TPAT_Correction_Engine::setScaleVars(bool b){ scaleVars = b; }

/**
 *	@brief Set the error tolerance
 *	@param d errors below this value will be considered negligible
 */
void TPAT_Correction_Engine::setTol(double d){
	tol = d;

	if(tol > 1)
		printWarn("TPAT_Correction_Engine::setTol: tolerance is greater than 1... just FYI\n");
}//====================================================

/**
 *	@brief Set the findEven flag
 *	@param b whether or not the algorithm will be looking for an event
 */
void TPAT_Correction_Engine::setFindEvent(bool b){ findEvent = b; }

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Correct a generic nodeset using multiple shooting
 *	@details This algorithm employs multiple shooting to correct a set of nodes
 *	subject to a set of constraints. The nodes and constraints are all stored in the 
 *	input nodeset object. 
 *	
 *	Numerical integration is used to generate state transition matrices between
 *	nodes. In order to optimize performance and accuracy, only two steps are computed
 *	by the simulation engine, and the step size is fixed to force the usage of the 
 *	Adams-Bashforth Adams-Moulton method.
 *	
 *	@param set pointer to the nodeset that needs to be corrected
 *	@param nodesOut pointer to the nodeset object that will contain the results of
 *	the shooting process
 *	@return the iteration data object for this corrections process
 *	@throws TPAT_Diverge if the corrections process does not converge
 *	@throws TPAT_Exception
 *	* if the input and output nodesets contain different system data objects
 *	* if the dynamic model associated with the input
 *	nodeset does not support one or more of the nodeset's constraints
 *	* if the input nodeset contains more than one delta-v constraint
 *	* if the input nodeset contains more than one TOF constraint
 */
TPAT_MultShoot_Data TPAT_Correction_Engine::multShoot(const TPAT_Nodeset *set, TPAT_Nodeset *nodesOut){
	if(nodesOut != NULL && *(set->getSysData()) != *(nodesOut->getSysData()))
		throw TPAT_Exception("TPAT_Correction_Engine::multShoot: Input and Output nodesets must use the same system data object");

	if(!isClean)
		cleanEngine();

	isClean = false;

	// Create structure to store iteration data for easy sharing
	TPAT_MultShoot_Data it(set);
	it.varTime = varTime;	// Save in structure to pass easily to other functions
	it.equalArcTime = equalArcTime;
	
	printVerb(verbose == TPAT_Verbosity_Tp::ALL_MSG, "Multiple Shooting Algorithm:\n");
	printVerb(verbose == TPAT_Verbosity_Tp::ALL_MSG, "  it.numNodes = %d\n", it.numNodes);
	printVerb(verbose == TPAT_Verbosity_Tp::ALL_MSG, "  sysType = %s\n", set->getSysData()->getTypeStr().c_str());

	// Get the model associated with the nodeset
	const TPAT_Model *model = set->getSysData()->getModel();
	model->multShoot_initDesignVec(&it, set);

	// Set up scaling
	it.freeVarScale.assign(4, 1);	// Assign all variable scalings to be one -> NOTE: ADD MORE ENTRIES IF YOU NEED MORE!!
	if(scaleVars)
		model->multShoot_scaleDesignVec(&it, set);

	// Create constraints that enforce continuity between nodes; this process
	// does account for velocity discontinuities specified in the nodeset
	it.allCons.clear();
	model->multShoot_createContCons(&it, set);

	// Add all node constraints
	for(int n = 0; n < set->getNumNodes(); n++){
		std::vector<TPAT_Constraint> nodeCons = set->getNodeByIx(n).getConstraints();
		it.allCons.insert(it.allCons.end(), nodeCons.begin(), nodeCons.end());
	}

	// Add all segment constraints
	for(int s = 0; s < set->getNumSegs(); s++){
		std::vector<TPAT_Constraint> segCons = set->getSegByIx(s).getConstraints();
		it.allCons.insert(it.allCons.end(), segCons.begin(), segCons.end());
	}

	// Add all arcset constraints
	std::vector<TPAT_Constraint> arcCons = set->getArcConstraints();
	it.allCons.insert(it.allCons.end(), arcCons.begin(), arcCons.end());

	// Compute number of extra consraint functions to add
	it.numSlack = 0;

	// Initialize vector to keep track of which row each constraint begins on
	// Also add slack variables to the design variable vector
	it.conRows.assign(it.allCons.size(), -1);	// Fill with -1
	int conRow = 0;
	bool foundDVCon = false;
	bool foundTOFCon = false;
	for(size_t c = 0; c < it.allCons.size(); c++){
		int addToRows = 0;
		TPAT_Constraint con = it.allCons[c];

		if(!model->supportsCon(con.getType()))
			throw TPAT_Exception("TPAT_Correction_Engine::multShoot: The dynamic model does not support one of the constraints!");

		switch(con.getType()){
			case TPAT_Constraint_Tp::CONT_PV:
			case TPAT_Constraint_Tp::CONT_EX:
			case TPAT_Constraint_Tp::SEG_CONT_PV:
			case TPAT_Constraint_Tp::SEG_CONT_EX:
			case TPAT_Constraint_Tp::STATE:
			case TPAT_Constraint_Tp::MATCH_CUST:
				addToRows = con.countConstrainedStates();
				break;
			case TPAT_Constraint_Tp::MATCH_ALL:
				addToRows = 6;
				break;
			case TPAT_Constraint_Tp::SP:
				addToRows = 3;
				break;
			case TPAT_Constraint_Tp::SP_RANGE:
				addToRows = 1;
				it.X.push_back(model->multShoot_getSlackVarVal(&it, con));
				it.slackAssignCon.push_back(c);
				it.numSlack++;
				break;
			case TPAT_Constraint_Tp::SP_MAX_DIST:
				it.X.push_back(model->multShoot_getSlackVarVal(&it, con));
				it.slackAssignCon.push_back(c);
				it.numSlack++;
			case TPAT_Constraint_Tp::SP_DIST:
				addToRows = 1;
				break;
			case TPAT_Constraint_Tp::MAX_DIST:
			case TPAT_Constraint_Tp::MIN_DIST:
				it.X.push_back(model->multShoot_getSlackVarVal(&it, con));
				it.slackAssignCon.push_back(c);	// remember where this slack variable is hiding
				it.numSlack++;
				// do NOT break here, continue on to do stuff for TPAT_Constraint_Tp::DIST as well
			case TPAT_Constraint_Tp::DIST:
				addToRows = 1;
				break;
			case TPAT_Constraint_Tp::MAX_DELTA_V:
			case TPAT_Constraint_Tp::DELTA_V:
				if(!foundDVCon){
					addToRows = 1;
					foundDVCon = true;

					if(con.getType() == TPAT_Constraint_Tp::MAX_DELTA_V){
						/* Add a slack variable to the end of design vector and keep track
						 * of which constraint it is assigned to; value of slack
						 * variable will be recomputed later
						 */
						it.X.push_back(model->multShoot_getSlackVarVal(&it, con));
						it.numSlack++;
						it.slackAssignCon.push_back(c);
					}
				}else{
					throw TPAT_Exception("TPAT_Correction_Engine::multShoot: You can only apply ONE delta-V constraint");
				}
				break;
			case TPAT_Constraint_Tp::JC:
				addToRows = 1;
				break;
			case TPAT_Constraint_Tp::TOF:
				if(!varTime)
					printWarn("TPAT_Correction_Engine::multShoot: Attempting to constraint TOF without variable time... won't work!");
				
				if(!foundTOFCon)
					addToRows = 1;
				else
					throw TPAT_Exception("TPAT_Correction_Engine::multShoot: You can only apply ONE TOF constraint");
				break;
			case TPAT_Constraint_Tp::APSE:
				addToRows = 1;
				break;
			case TPAT_Constraint_Tp::PSEUDOARC:
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

	printVerb(verbose == TPAT_Verbosity_Tp::ALL_MSG, "  # Free: %d\n  # Constraints: %d\n", it.totalFree, it.totalCons);
	printVerb(verbose == TPAT_Verbosity_Tp::ALL_MSG, "  -> # Slack Variables: %d\n", it.numSlack);

	printVerb(verbose == TPAT_Verbosity_Tp::ALL_MSG, "ALL CONSTRAINTS:\n\n");
	if(verbose == TPAT_Verbosity_Tp::ALL_MSG){
		for(size_t n = 0; n < it.allCons.size(); n++){
			it.allCons[n].print();
		}
	}
	
	// Run the multiple shooting process
	return multShoot(it, nodesOut);
}//==========================================================

/**
 *  @brief Run a multiple shooting algorithm given an TPAT_MultShoot_Data object
 * 
 *  @param it A completely formed TPAT_MultShoot_Data object that describes a 
 *  multiple shooting problem. These are created from TPAT_Nodeset and its
 *  derivative types by the other implementation of multShoot()
 *  @param nodesOut pointer to a nodeset object that will contain the results 
 *  of the shooting process
 *  @return A corrected TPAT_MultShoot_Data object
 *  @see multShoot(TPAT_Nodeset*)
 *  @throws TPAT_Diverge if the multiple shooting process does not converge
 */
TPAT_MultShoot_Data TPAT_Correction_Engine::multShoot(TPAT_MultShoot_Data it, TPAT_Nodeset *nodesOut){
	it.count = 0;

	// create a simulation engine
	TPAT_Sim_Engine simEngine;
	simEngine.setVerbose(verbose);
	
	// Set both tolerances of simulation engine to be three orders of magnitude less corrector
	double simTol = tol/1000 < 1e-15 ? 1e-15 : tol/1000;
	simEngine.setAbsTol(simTol);
	simEngine.setRelTol(simTol);

	// Only need info about the final state, no need to generate lots of intermediate data points
	// This forces the integrator to use the Adams-Bashforth Adams-Moulton integration method,
	// which is most similar to Matlab's ode113
	simEngine.setVarStepSize(false);
	simEngine.setNumSteps(2);

	if(findEvent || ignoreCrash){
		simEngine.clearEvents();	// don't use crash events when searching for an event
	}

	// Define values for use in corrections loop
	double err = 10*tol;
	while( err > tol && it.count < maxIts){
		it.FX.clear();					// Clear vectors each iteration
		it.DF.clear();
		it.deltaVs.clear();
		it.propSegs.clear();
		it.FX.assign(it.totalCons, 0);	// Size the vectors and fill with zeros
		it.DF.assign(it.totalCons*it.totalFree, 0);
		it.deltaVs.assign(3*it.numNodes, 0);

		// initialize a vector of trajectory objects to store each propagated segment
		it.sysData->getModel()->multShoot_initIterData(&it);

		for(int s = 0; s < it.nodeset->getNumSegs(); s++){
			// printf("Retrieving ICs for segment (ix %02d):\n", s);
			// Get simulation conditions from design vector via dynamic model implementation
			double t0 = 0, tof = 0;
			double ic[] = {0,0,0,0,0,0};
			it.sysData->getModel()->multShoot_getSimICs(&it, it.nodeset, it.nodeset->getSegByIx(s).getID(),
				ic, &t0, &tof);

			simEngine.setRevTime(tof < 0);
			simEngine.runSim(ic, t0, tof, &(it.propSegs[s]));
		}

		// waitForUser();
		
		// Compute Delta-Vs between node segments
		for(int s = 0; s < it.nodeset->getNumSegs(); s++){
			std::vector<double> lastState = it.propSegs[s].getStateByIx(-1);
			int termID = it.nodeset->getSegByIx(s).getTerminus();
			if(termID != TPAT_Linkable::INVALID_ID){
				int termNodeIx = it.nodeset->getNodeIx(termID);
				// velCon has false for a velocity state if there is a discontinuity between
				// the terminus of the segment and the terminal node
				std::vector<bool> velCon = it.nodeset->getSegByIx(s).getVelCon();
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
		for(size_t c = 0; c < it.allCons.size(); c++)
			it.sysData->getModel()->multShoot_applyConstraint(&it, it.allCons[c], c);

		// Solve for newX and copy into working vector X
		Eigen::VectorXd oldX = Eigen::Map<Eigen::VectorXd>(&(it.X[0]), it.totalFree, 1);
		Eigen::VectorXd newX = solveUpdateEq(&it);

		it.X.clear();
		it.X.insert(it.X.begin(), newX.data(), newX.data()+it.totalFree);

		// Compute error; norm of constraint vector
		Eigen::VectorXd FX = Eigen::Map<Eigen::VectorXd>(&(it.FX[0]), it.totalCons, 1);
		double err_cons = FX.norm();

		if(verbose == TPAT_Verbosity_Tp::ALL_MSG)
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
		printVerbColor((findEvent && verbose == TPAT_Verbosity_Tp::ALL_MSG) || (!findEvent && verbose > TPAT_Verbosity_Tp::NO_MSG), YELLOW, "Iteration %02d: err = %.4e (%s)\n",
			it.count, err, errType.c_str());
	}// end of corrections loop

	if(err > tol && !ignoreDiverge){
		throw TPAT_Diverge();
	}

	if(nodesOut != NULL){
		try{
			it.sysData->getModel()->multShoot_createOutput(&it, it.nodeset, findEvent, nodesOut);
		}catch(TPAT_Exception &e){
			printErr("TPAT_Correction_Engine::multShoot: Unable to create output nodeset\n  Err: %s\n", e.what());
			throw e;
		}
	}
	
	return it;
}//=====================================================

/**
 *	@brief Apply linear algebra to solve the update equation and obtain an updated free-variable vector
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
 *	@param it the TPAT_MultShoot_Data object associated with the corrections process
 *
 *	@return the updated free variable vector \f$ \vec{X}_{n+1} \f$
 *	@throws TPAT_Exception if the problem is over constrained (i.e. Jacobian has more rows than columns);
 *	This can be updated to use a least-squares solution (TODO)
 */
Eigen::VectorXd TPAT_Correction_Engine::solveUpdateEq(TPAT_MultShoot_Data* it){
	// Create matrices for X, Jacobian matrix DF, and constraint vector FX
	Eigen::VectorXd oldX = Eigen::Map<Eigen::VectorXd>(&(it->X[0]), it->totalFree, 1);
	MatrixXRd J = Eigen::Map<MatrixXRd>(&(it->DF[0]), it->totalCons, it->totalFree);
	Eigen::VectorXd FX = Eigen::Map<Eigen::VectorXd>(&(it->FX[0]), it->totalCons, 1);

	// change sign for matrix multiplication
	FX *= -1;

	// Create a vector to put the solution in
	Eigen::VectorXd X_diff(it->totalFree, 1);
	if(it->totalCons == it->totalFree){	// J is square, use regular inverse

		/* Use LU decomposition to invert the Gramm matrix and find a vector
		w. Multiplying J^T by w yields the minimum-norm solution x, where x 
		lies in the column-space of J^T, or in the orthogonal complement of
		the nullspace of J.
		Source: <http://www.math.usm.edu/lambers/mat419/lecture15.pdf>
		 */
		
		// Solve the system Jw = b (In this case, w = X_diff)
		Eigen::FullPivLU<MatrixXRd> lu(J);
		// lu.setThreshold(1e-20);

		// if(!lu.isInvertible())
		// 	throw TPAT_LinAlg_Err("TPAT_Correction_Engine::solveUpdateEq: Jacobian is singular; cannot solve");

		X_diff = lu.solve(FX);
	}else{
		if(it->totalCons < it->totalFree){	// Under-constrained
			// Compute Gramm matrix
			MatrixXRd JT = J.transpose();
			MatrixXRd G = J*JT;

			// toCSV(J, "DF_cpp.csv");
			// toCSV(FX, "FX_cpp.csv");
			// toCSV(oldX, "X_cpp.csv");

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
			// 	printErr("Gramm Matrix rank = %ld / %ld\n", lu.rank(), G.rows());
			// 	throw TPAT_LinAlg_Err("TPAT_Correction_Engine::solveUpdateEq: Gramm matrix is singular; cannot solve");
			// }

			Eigen::VectorXd w = lu.solve(FX);
			
			// Compute optimal x from w
			X_diff = JT*w;
		}else{	// Over-constrained
			throw TPAT_LinAlg_Err("System is over constrained... No solution implemented");
		}
	}

	return oldX + X_diff;	// newX = oldX + X_diff
}// End of solveUpdateEq() =====================================

/**
 *  @brief Print out the magnitude of each constraint.
 *  @details This can be useful when debugging to highlight which constraints are unsatisfied
 * 
 *  @param it pointer to an TPAT_MultShoot_Data object associated with a corrections process
 */
void TPAT_Correction_Engine::reportConMags(const TPAT_MultShoot_Data *it){
	int conCount = 0;
	for(long r = 0; r < (int)(it->FX.size()); r++){
        if(r == 0 && it->totalCons > 0){
            printf("Node %d %s Constraint:\n", it->allCons[conCount].getID(), it->allCons[conCount].getTypeStr());
        }else if(conCount < (int)(it->allCons.size()) && r >= it->conRows[conCount+1]){
            conCount++;
            printf("Node %d %s Constraint:\n", it->allCons[conCount].getID(), it->allCons[conCount].getTypeStr());
        }
        printf("  ||row %03zu||: %.6e\n", r, std::abs(it->FX[r]));
    }
}//===============================================================

/**
 *	@brief clean up data so that engine can be used again (or deconstructed safely)
 */
void TPAT_Correction_Engine::cleanEngine(){
	printVerb(verbose == TPAT_Verbosity_Tp::ALL_MSG, "Cleaning the engine...\n");
	isClean = true;
}//====================================================





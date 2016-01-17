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
#include "tpat.hpp"
 
#include "tpat_correction_engine.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_calculations.hpp"
#include "tpat_exceptions.hpp" 
#include "tpat_nodeset_bcr4bp.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"
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
tpat_correction_engine::tpat_correction_engine(const tpat_correction_engine &e){
	copyEngine(e);
}//=======================================================

/**
 *	@brief Destructor
 */
tpat_correction_engine::~tpat_correction_engine(){
	cleanEngine();
}//=================================================

/**
 *	@brief Copy all engine variables
 *	@param e an engine reference
 */
void tpat_correction_engine::copyEngine(const tpat_correction_engine &e){
	verbose = e.verbose;
	varTime = e.varTime;
	equalArcTime = e.equalArcTime;
	maxIts = e.maxIts;
	tol = e.tol;
	createdNodesetOut = e.createdNodesetOut;
	findEvent = e.findEvent;
	ignoreCrash = e.ignoreCrash;
	
	if(createdNodesetOut){
		tpat_sys_data::system_t type = e.nodeset_out->getSysData()->getType();
		switch(type){
			case tpat_sys_data::CR3BP_SYS:
				nodeset_out = new tpat_nodeset_cr3bp (* static_cast<tpat_nodeset_cr3bp *>(e.nodeset_out));
				break;
			case tpat_sys_data::BCR4BPR_SYS:
				nodeset_out = new tpat_nodeset_bcr4bp (*static_cast<tpat_nodeset_bcr4bp *>(e.nodeset_out));
				break;
			default: nodeset_out = 0; break;
		}
	}else{
		nodeset_out = 0;
	}
}//=================================================

//-----------------------------------------------------
//      Operator Functions
//-----------------------------------------------------

/**
 *	@brief Copy operator; make a copy of the input correction engine. The dynamically allocated
 *	<tt>nodeset_out</tt> is copied, if possible, or set to 0 (it's a pointer) otherwise.
 *
 *	@param e
 *	@return this correction engine
 */
tpat_correction_engine& tpat_correction_engine::operator =(const tpat_correction_engine &e){
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
bool tpat_correction_engine::usesVarTime() const { return varTime; }

/**
 *	@brief Retrieve whether or not we force all segments to have the same length
 *	(in time).
 *
 *	This setting only applies if variable time is turned on.
 *	@return whether or not each arc will be forced to have the same length in time
 */
bool tpat_correction_engine::usesEqualArcTime() const { return equalArcTime; }

/**
 *  @brief Retrieve whether or not the multiple shooting algorithm uses variable scaling
 *  @return whether or not the multiple shooting algorithm uses variable scaling
 */
bool tpat_correction_engine::usesScaledVars() const { return scaleVars; }

/**
 *  @brief Retrieve the verbosity setting
 *	@return whether or not the corrector will be verbose
 */
verbosity_t tpat_correction_engine::isVerbose() const { return verbose; }

/**
 *  @brief Retrieve whether or not we are located an event crossing
 *	@return whether or not the algorithm will optimize the process to find an event
 */
bool tpat_correction_engine::isFindingEvent() const { return findEvent; }

/**
 *  @brief Retrieve the maximum number of iterations to attempt
 *	@return the maximum number of iterations to attempt before giving up
 */
int tpat_correction_engine::getMaxIts() const { return maxIts; }

/**
 *  @brief Retrieve the minimum error tolerance
 *	@return the minimum error tolerance (non-dimensional units); errors
 *	less than this value are considered negligible
 */
double tpat_correction_engine::getTol() const { return tol; }

/**
 *	@brief Retrieve the output CR3BP nodeset (after corrections). 
 *
 *	Note that this method will throw an
 *	error if the corrections process has not been run or failed to produce an output.
 *
 *	@return a CR3BP nodeset object with the corrected trajectory data stored inside
 */
tpat_nodeset_cr3bp tpat_correction_engine::getCR3BP_Output(){
	if(createdNodesetOut){
		if(nodeset_out->getSysData()->getType() == tpat_sys_data::CR3BP_SYS){
			// Create a copy of the nodeset, return it
			tpat_nodeset_cr3bp temp( *(static_cast<tpat_nodeset_cr3bp *>(nodeset_out)) );
			return temp;
		}else{
			throw tpat_exception("tpat_correction_engine::getCR3BP_Output: Wrong system type");
		}
	}else{
		throw tpat_exception("tpat_correction_engine::getCR3BP_Output: Output nodeset has not been created, cannot return CR3BP output");
	}
}//=========================================

/**
 *	@brief Retrieve the output BCR4BPR nodeset (after corrections). 
 *
 *	Note that this method will throw an
 *	error if the corrections process has not been run or failed to produce an output.
 *
 *	@return a BCR4BPR nodeset object with the corrected trajectory data stored inside
 */
tpat_nodeset_bcr4bp tpat_correction_engine::getBCR4BPR_Output(){
	if(createdNodesetOut){
		if(nodeset_out->getSysData()->getType() == tpat_sys_data::BCR4BPR_SYS){
			// Create a copy of the nodeset, return it
			tpat_nodeset_bcr4bp temp( *(static_cast<tpat_nodeset_bcr4bp *>(nodeset_out)) );
			return temp;
		}else{
			throw tpat_exception("tpat_correction_engine::getBCR4BPR_Output: Wrong system type");
		}
	}else{
		throw tpat_exception("tpat_correction_engine::getBCR4BPR_Output: Output nodeset has not been created, cannot return CR3BP output");
	}
}//==================================================

/**
 *	@brief Set varTime
 *	@param b whether or not the corrector should use variable time
 */
void tpat_correction_engine::setVarTime(bool b){
	varTime = b;
	// Turn off equal-time arcs too if varTime is false
	if(!varTime)
		equalArcTime = false;
}//==================================================

/**
 *	@brief Tell the corrector how to apply variable time
 *	@param b whether or not each arc will be forced to have the same duration
 */
void tpat_correction_engine::setEqualArcTime(bool b){
	if(!varTime && b){
		printErr("tpat_correction_engine::setequalArcTime: Cannot use equal-time arcs if variable time is disabled; please turn varTime ON first\n");
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
void tpat_correction_engine::setIgnoreCrash(bool b){ ignoreCrash = b; }

/**
 *  @brief Tell the corrector to ignore divergence and return the partially
 *  corrected iteration data instead of throwing an exception when divergence
 *  occurs.
 * 
 *  @param b Whether or not to ignore divergance
 */
void tpat_correction_engine::setIgnoreDiverge(bool b){ ignoreDiverge = b;}

/**
 *	@brief Set verbosity
 *	@param b whether or not the corrector should be verbose in its outputs
 */
void tpat_correction_engine::setVerbose(verbosity_t b){ verbose = b; }

/**
 *	@brief Set maximum iterations
 *	@param i the maximum number of iterations to attempt before giving up
 */
void tpat_correction_engine::setMaxIts(int i){ maxIts = i; }

/**
 *  @brief Set the scaleVar flag
 * 
 *  @param b whether or not the multiple shooting algorithm should use variable scaling
 */
void tpat_correction_engine::setScaleVars(bool b){ scaleVars = b; }

/**
 *	@brief Set the error tolerance
 *	@param d errors below this value will be considered negligible
 */
void tpat_correction_engine::setTol(double d){
	tol = d;

	if(tol > 1)
		printWarn("tpat_correction_engine::setTol: tolerance is greater than 1... just FYI\n");
}

/**
 *	@brief Set the findEven flag
 *	@param b whether or not the algorithm will be looking for an event
 */
void tpat_correction_engine::setFindEvent(bool b){ findEvent = b; }

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
 *	@param set a pointer to a nodeset
 *	@return the iteration data object for this corrections process
 */
iterationData tpat_correction_engine::multShoot(tpat_nodeset *set){
	if(!isClean)
		cleanEngine();

	isClean = false;

	// Make sure all constraints have the propper node numbers
	set->updateCons();

	// Create structure to store iteration data for easy sharing
	iterationData it;
	it.varTime = varTime;	// Save in structure to pass easily to other functions
	it.equalArcTime = equalArcTime;
	it.sysData = set->getSysData();
	it.nodeset = set;

	// Save original nodes for later access (particularly when variable time is off)
	for(int n = 0; n < set->getNumNodes(); n++){
		it.origNodes.push_back(set->getNode(n));
	}

	// Get some basic data from the input nodeset
	it.numNodes = set->getNumNodes();
	
	printVerb(verbose == ALL_MSG, "Multiple Shooting Algorithm:\n");
	printVerb(verbose == ALL_MSG, "  it.numNodes = %d\n", it.numNodes);
	printVerb(verbose == ALL_MSG, "  sysType = %s\n", set->getSysData()->getTypeStr().c_str());

	// Get the model associated with the nodeset
	const tpat_model *model = set->getSysData()->getModel();
	model->multShoot_initDesignVec(&it, set);

	// Set up scaling
	it.freeVarScale.assign(4, 1);	// Assign all variable scalings to be one -> NOTE: ADD MORE ENTRIES IF YOU NEED MORE!!
	if(scaleVars)
		model->multShoot_scaleDesignVec(&it, set);

	// Create constraints that enforce continuity between nodes; this process
	// does account for velocity discontinuities specified in the nodeset
	it.allCons.clear();
	model->multShoot_createContCons(&it, set);

	// Add all extra constraints from the nodeset to the total constraint vector
	for(int n = 0; n < set->getNumNodes(); n++){
		std::vector<tpat_constraint> nodeCons = set->getNodeCons(n);
		for(size_t c = 0; c < nodeCons.size(); c++){
			it.allCons.push_back(nodeCons[c]);
		}
	}

	// Compute number of extra consraint functions to add
	it.numSlack = 0;

	// Initialize vector to keep track of which row each constraint begins on
	// Also add slack variables to the 
	it.conRows.assign(it.allCons.size(), NAN);	// Fill with NAN
	int conRow = 0;
	bool foundDVCon = false;
	bool foundTOFCon = false;
	for(size_t c = 0; c < it.allCons.size(); c++){
		int addToRows = 0;
		tpat_constraint con = it.allCons[c];

		if(!model->supportsCon(con.getType()))
			throw tpat_exception("tpat_correction_engine::correct: The dynamic model does not support one of the constraints!");

		switch(con.getType()){
			case tpat_constraint::CONT_PV:
			case tpat_constraint::CONT_EX:
			case tpat_constraint::STATE:
			case tpat_constraint::MATCH_CUST:
				addToRows = con.countConstrainedStates();
				break;
			case tpat_constraint::MATCH_ALL:
				addToRows = 6;
				break;
			case tpat_constraint::SP:
				addToRows = 3;
				break;
			case tpat_constraint::SP_RANGE:
				addToRows = 1;
				it.X.push_back(model->multShoot_getSlackVarVal(&it, con));
				it.slackAssignCon.push_back(c);
				it.numSlack++;
				break;
			case tpat_constraint::SP_MAX_DIST:
				it.X.push_back(model->multShoot_getSlackVarVal(&it, con));
				it.slackAssignCon.push_back(c);
				it.numSlack++;
			case tpat_constraint::SP_DIST:
				addToRows = 1;
				break;
			case tpat_constraint::MAX_DIST:
			case tpat_constraint::MIN_DIST:
				it.X.push_back(model->multShoot_getSlackVarVal(&it, con));
				it.slackAssignCon.push_back(c);	// remember where this slack variable is hiding
				it.numSlack++;
				// do NOT break here, continue on to do stuff for DIST as well
			case tpat_constraint::DIST:
				addToRows = 1;
				break;
			case tpat_constraint::MAX_DELTA_V:
			case tpat_constraint::DELTA_V:
				if(!foundDVCon){
					addToRows = 1;
					foundDVCon = true;

					if(con.getType() == tpat_constraint::MAX_DELTA_V){
						/* Add a slack variable to the end of design vector and keep track
						 * of which constraint it is assigned to; value of slack
						 * variable will be recomputed later
						 */
						it.X.push_back(model->multShoot_getSlackVarVal(&it, con));
						it.numSlack++;
						it.slackAssignCon.push_back(c);
					}
				}else{
					throw tpat_exception("You can only apply ONE delta-V constraint");
				}
				break;
			case tpat_constraint::JC:
				addToRows = 1;
				break;
			case tpat_constraint::TOF:
				if(!varTime)
					printWarn("Attempting to constraint TOF without variable time... won't work!");
				
				if(!foundTOFCon)
					addToRows = 1;
				else
					throw tpat_exception("You can only apply ONE TOF constraint");
				break;
			case tpat_constraint::APSE:
				addToRows = 1;
				break;
			case tpat_constraint::PSEUDOARC:
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

	printVerb(verbose == ALL_MSG, "  # Free: %d\n  # Constraints: %d\n", it.totalFree, it.totalCons);
	printVerb(verbose == ALL_MSG, "  -> # Slack Variables: %d\n", it.numSlack);

	printVerb(verbose == ALL_MSG, "ALL CONSTRAINTS:\n\n");
	if(verbose == ALL_MSG){
		for(size_t n = 0; n < it.allCons.size(); n++){
			it.allCons[n].print();
		}
	}
	
	// Run the multiple shooting process
	return multShoot(it);
}//==========================================================

/**
 *  @brief Run a multiple shooting algorithm given an iterationData object
 * 
 *  @param it A completely formed iterationData object that describes a 
 *  multiple shooting problem. These are created from tpat_nodeset and its
 *  derivative types by the other implementation of multShoot()
 *  @return A corrected iterationData object
 *  @see multShoot(tpat_nodeset*)
 */
iterationData tpat_correction_engine::multShoot(iterationData it){
	it.count = 0;

	// create a simulation engine
	tpat_simulation_engine simEngine(it.sysData);
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
	double err = 1e10;
	while( err > tol && it.count < maxIts){
		it.FX.clear();					// Clear vectors each iteration
		it.DF.clear();
		it.deltaVs.clear();
		it.allSegs.clear();
		it.FX.assign(it.totalCons, 0);	// Size the vectors and fill with zeros
		it.DF.assign(it.totalCons*it.totalFree, 0);
		it.deltaVs.assign(3*it.numNodes, 0);

		for(int n = 0; n < it.numNodes-1; n++){
			// Get simulation conditions from design vector via dynamic model implementation
			double t0 = 0;
			double tof = 0;
			double ic[] = {0,0,0,0,0,0};
			it.sysData->getModel()->multShoot_getSimICs(&it, it.nodeset, n, ic, &t0, &tof);

			simEngine.setRevTime(tof < 0);
			simEngine.runSim(ic, t0, tof);
			it.allSegs.push_back(simEngine.getTraj());
		}// end of loop through nodes

		// Compute Delta-Vs between node segments
		for(int n = 0; n < it.numNodes - 1; n++){
			std::vector<double> lastState = it.allSegs[n].getState(-1);
			// velCon has false for a velocity state if there is a discontinuity between v_n,f and v_n+1
			std::vector<bool> velCon = it.nodeset->getNode(n+1).getVelCon();
			for(int s = 3; s < 6; s++){
				// Compute difference in velocity; if velCon[s-3] is true, then velocity
				// should be continuous and any difference is numerical error, so set to
				// zero by multiplying by not-true
				it.deltaVs[n*3+s-3] = !velCon[s-3]*(lastState[s] - it.X[6*(n+1)+s]);
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

		if(verbose == ALL_MSG)
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
		printVerbColor((findEvent && verbose == ALL_MSG) || (!findEvent && verbose > NO_MSG), YELLOW, "Iteration %02d: err = %.4e (%s)\n",
			it.count, err, errType.c_str());
	}// end of corrections loop

	if(err > tol && !ignoreDiverge){
		throw tpat_diverge();
	}

	nodeset_out = it.sysData->getModel()->multShoot_createOutput(&it, it.nodeset, findEvent);
	createdNodesetOut = true;

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
 *	@param it the iterationData object associated with the corrections process
 *
 *	@return the updated free variable vector \f$ \vec{X}_{n+1} \f$
 */
Eigen::VectorXd tpat_correction_engine::solveUpdateEq(iterationData* it){
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
		// 	throw tpat_linalg_err("tpat_correction_engine::solveUpdateEq: Jacobian is singular; cannot solve");

		X_diff = lu.solve(FX);
	}else{
		if(it->totalCons < it->totalFree){	// Under-constrained
			// Compute Gramm matrix
			MatrixXRd JT = J.transpose();
			MatrixXRd G = J*JT;

			toCSV(J, "DF_cpp.csv");
			toCSV(FX, "FX_cpp.csv");
			toCSV(oldX, "X_cpp.csv");

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
			// 	throw tpat_linalg_err("tpat_correction_engine::solveUpdateEq: Gramm matrix is singular; cannot solve");
			// }

			Eigen::VectorXd w = lu.solve(FX);
			
			// Compute optimal x from w
			X_diff = JT*w;
		}else{	// Over-constrained
			throw tpat_linalg_err("System is over constrained... No solution implemented");
		}
	}

	return oldX + X_diff;	// newX = oldX + X_diff
}// End of solveUpdateEq() =====================================

void tpat_correction_engine::reportConMags(iterationData *it){
	int conCount = 0;
	for(long r = 0; r < (int)(it->FX.size()); r++){
        if(r == 0 && it->totalCons > 0){
            printf("Node %d %s Constraint:\n", it->allCons[conCount].getNode(), it->allCons[conCount].getTypeStr());
        }else if(conCount < (int)(it->allCons.size()) && r >= it->conRows[conCount+1]){
            conCount++;
            printf("Node %d %s Constraint:\n", it->allCons[conCount].getNode(), it->allCons[conCount].getTypeStr());
        }
        printf("  ||row %03zu||: %.6e\n", r, std::abs(it->FX[r]));
    }
}//===============================================================

/**
 *	@brief clean up data so that engine can be used again (or deconstructed safely)
 */
void tpat_correction_engine::cleanEngine(){
	printVerb(verbose == ALL_MSG, "Cleaning the engine...\n");
	if(createdNodesetOut){
		delete nodeset_out;
	}

	createdNodesetOut = false;
	isClean = true;
}//===================================
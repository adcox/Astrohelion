/**
 *	@file tpat_correction_engine.cpp
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
#include "tpat_nodeset_bcr4bpr.hpp"
#include "tpat_traj_bcr4bpr.hpp"
#include "tpat_calculations.hpp"
#include "tpat_constraint.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_matrix.hpp"
#include "tpat_nodeset.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data.hpp"
#include "tpat_traj.hpp"
#include "tpat_utilities.hpp"

#include <algorithm>
#include <cmath>
#include <gsl/gsl_linalg.h>
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
	maxIts = e.maxIts;
	tol = e.tol;
	receivedNodesetIn = e.receivedNodesetIn;
	createdNodesetOut = e.createdNodesetOut;
	findEvent = e.findEvent;
	nodeset_in = e.nodeset_in;		//POINTER, COPYING ADDRESS - passed in, so should point to same parent object

	if(createdNodesetOut){
		tpat_sys_data::system_t type = e.nodeset_out->getSysData()->getType();
		switch(type){
			case tpat_sys_data::CR3BP_SYS:
				nodeset_out = new tpat_nodeset_cr3bp (* static_cast<tpat_nodeset_cr3bp *>(e.nodeset_out));
				break;
			case tpat_sys_data::BCR4BPR_SYS:
				nodeset_out = new tpat_nodeset_bcr4bpr (*static_cast<tpat_nodeset_bcr4bpr *>(e.nodeset_out));
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
 *  @brief Retrieve the verbosity setting
 *	@return whether or not the corrector will be verbose
 */
bool tpat_correction_engine::isVerbose() const { return verbose; }

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
	if(createdNodesetOut && receivedNodesetIn){
		if(nodeset_in->getSysData()->getType() == tpat_sys_data::CR3BP_SYS){
			// Create a copy of the nodeset, return it
			tpat_nodeset_cr3bp temp( *(static_cast<tpat_nodeset_cr3bp *>(nodeset_out)) );
			return temp;
		}else{
			throw tpat_exception("Wrong system type");
		}
	}else{
		throw tpat_exception("Output nodeset has not been created, cannot return CR3BP output");
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
tpat_nodeset_bcr4bpr tpat_correction_engine::getBCR4BPR_Output(){
	if(createdNodesetOut && receivedNodesetIn){
		if(nodeset_in->getSysData()->getType() == tpat_sys_data::BCR4BPR_SYS){
			// Create a copy of the nodeset, return it
			tpat_nodeset_bcr4bpr temp( *(static_cast<tpat_nodeset_bcr4bpr *>(nodeset_out)) );
			return temp;
		}else{
			throw tpat_exception("Wrong system type");
		}
	}else{
		throw tpat_exception("Output nodeset has not been created, cannot return CR3BP output");
	}
}//=========================================

/**
 *	@brief Set varTime
 *	@param b whether or not the corrector should use variable time
 */
void tpat_correction_engine::setVarTime(bool b){ varTime = b; }

/**
 *	@brief Set verbosity
 *	@param b whether or not the corrector should be verbose in its outputs
 */
void tpat_correction_engine::setVerbose(bool b){ verbose = b; }

/**
 *	@brief Set maximum iterations
 *	@param i the maximum number of iterations to attempt before giving up
 */
void tpat_correction_engine::setMaxIts(int i){ maxIts = i; }

/**
 *	@brief Set the error tolerance
 *	@param d errors below this value will be considered negligible
 */
void tpat_correction_engine::setTol(double d){ tol = d; }

/**
 *	@brief Set the findEven flag
 *	@param b whether or not the algorithm will be looking for an event
 */
void tpat_correction_engine::setFindEvent(bool b){ findEvent = b; }

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Correct a CR3BP nodeset
 *	@param set a pointer to a nodeset
 *	@throws a <tt>tpat_diverge</tt> exception if the corrector cannot converge on a solution
 */
void tpat_correction_engine::correct_cr3bp(tpat_nodeset_cr3bp* set){
	if(!isClean)
		cleanEngine();

	nodeset_in = set;
	receivedNodesetIn = true;
	correct(nodeset_in);
	isClean = false;
}//========================================

/**
 *	@brief Correct a BCR4BP nodeset
 *	@param set a pointer to a nodeset
 *	@throws a <tt>tpat_diverge</tt> exception if the corrector cannot converge on a solution
 */
void tpat_correction_engine::correct_bcr4bpr(tpat_nodeset_bcr4bpr* set){
	if(!isClean)
		cleanEngine();

	nodeset_in = set;
	receivedNodesetIn = true;
	correct(nodeset_in);
	isClean = false;
}//======================================

/**
 *	@brief Correct a generic nodeset; equipped to handle any type
 *	@param set a pointer to a nodeset
 */
void tpat_correction_engine::correct(tpat_nodeset *set){
	// Create structure to store iteration data for easy sharing
	iterationData it;
	it.varTime = varTime;	// Save in structure to pass easily to other functions
	it.sysData = set->getSysData();

	// Save original nodes for later access (particularly when variable time is off)
	for(int n = 0; n < set->getNumNodes(); n++){
		it.origNodes.push_back(set->getNode(n));
	}

	// Get some basic data from the input nodeset
	it.numNodes = set->getNumNodes();
	tpat_sys_data::system_t sysType = set->getSysData()->getType();
	
	printVerb(verbose, "Corrector:\n");
	printVerb(verbose, "  it.numNodes = %d\n", it.numNodes);
	printVerb(verbose, "  sysType = %s\n", set->getSysData()->getTypeStr().c_str());

	// Create specific nodeset objects using static_cast
	tpat_nodeset_bcr4bpr *bcSet;
	if(sysType == tpat_sys_data::BCR4BPR_SYS)
		bcSet = static_cast<tpat_nodeset_bcr4bpr *>(set);

	// Dummy code for now; this will be replaced with better code later
	tpat_model *model = set->getSysData()->getModel();
	model->corrector_initDesignVec(&it, set);

	// Create constraints that enforce continuity between nodes; this process
	// does account for velocity discontinuities specified in the nodeset
	it.allCons.clear();
	model->corrector_createContCons(&it, set);

	// Add all extra constraints from the nodeset to the total constraint vector
	for(int i = 0; i < set->getNumCons(); i++){
		it.allCons.push_back(set->getConstraint(i));
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
			case tpat_constraint::MATCH_ALL:
			case tpat_constraint::MATCH_CUST:
				addToRows = con.countConstrainedStates();
				break;
			case tpat_constraint::SP:
				addToRows = 3;
				break;
			case tpat_constraint::MAX_DIST:
			case tpat_constraint::MIN_DIST:
				it.X.push_back(1e-4);			// initial value for slack variable
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
						/* Add a slack variable to the design vector and keep track
						 * of which constraint it is assigned to; value of slack
						 * variable will be recomputed later
						 */
						it.X.push_back(1e-4);
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
			default: break;
		}

		// Save the index of the first row for this constraint
		it.conRows[c] = conRow;
		conRow += addToRows;	// remember we've added rows
		it.totalCons += addToRows;
	}// END of loop through constraints

	// Determine the number of free/design variables based on the system type
	it.totalFree = it.X.size();

	printVerb(verbose, "  # Free: %d\n  # Constraints: %d\n", it.totalFree, it.totalCons);
	printVerb(verbose, "  -> # Slack Variables: %d\n", it.numSlack);

	printVerb(verbose, "ALL CONSTRAINTS:\n\n");
	if(verbose){
		for(size_t n = 0; n < it.allCons.size(); n++){
			it.allCons[n].print();
		}
	}
	// create a simulation engine
	tpat_sys_data *sysData = set->getSysData();
	tpat_simulation_engine simEngine(sysData);
	simEngine.setVerbose(verbose);
	
	// Set both tolerances of simulation engine to be same as corrector
	simEngine.setAbsTol(tol);
	simEngine.setRelTol(tol);

	// Only need info about the final state, no need to generate lots of intermediate data points
	simEngine.setVarStepSize(false);
	simEngine.setNumSteps(2);

	if(findEvent){
		simEngine.clearEvents();	// don't use crash events when searching for an event
	}

	// Define values for use in corrections loop
	double err = 1000000000;
	while( err > tol && it.count < maxIts){
		it.FX.clear();					// Clear vectors each iteration
		it.DF.clear();
		it.deltaVs.clear();
		it.allSegs.clear();
		it.FX.assign(it.totalCons, 0);	// Size the vectors and fill with zeros
		it.DF.assign(it.totalCons*it.totalFree, 0);
		it.deltaVs.assign(3*it.numNodes, 0);

		for(int n = 0; n < it.numNodes-1; n++){
			double t0 = 0;
			double tof = 0;
			double ic[] = {0,0,0,0,0,0};
			model->corrector_getSimICs(&it, set, n, ic, &t0, &tof);

			// Copy IC for node n from free variable vector
			// double *ic = &(it.X[6*n]);
			simEngine.setRevTime(tof < 0);
			simEngine.runSim(ic, t0, tof);
			it.allSegs.push_back(simEngine.getTraj());
		}// end of loop through nodes

		// Compute Delta-Vs between node segments
		for(int n = 0; n < set->getNumNodes()-1; n++){
			std::vector<double> lastState = it.allSegs[n].getState(-1);
			for(int s = 3; s < 6; s++){
				it.deltaVs[n*3+s-3] = lastState[s] - it.X[6*(n+1)+s];
			}
		}

		// Loop through all constraints and compute the constraint values, partials, and
		// apply them to the FX and DF matrices
		for(size_t c = 0; c < it.allCons.size(); c++)
			model->corrector_applyConstraint(&it, it.allCons[c], c);

		// Solve for newX and copy into working vector X
		tpat_matrix newX = solveUpdateEq(&it);
		it.X.clear();
		it.X.insert(it.X.begin(), newX.getDataPtr(), newX.getDataPtr()+it.totalFree);

		// Compute error; norm of constraint vector
		tpat_matrix FX_mat(it.totalCons, 1, it.FX);
		err = norm(FX_mat);

		it.count++;
		printVerbColor((!findEvent && !verbose) || verbose, YELLOW, "Iteration %02d: ||F|| = %.4e\n", it.count, err);
	}// end of corrections loop

	if(err > tol){
		throw tpat_diverge();
	}

	nodeset_out = model->corrector_createOutput(&it, nodeset_in, findEvent);
	createdNodesetOut = true;
}//==========================================================

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
tpat_matrix tpat_correction_engine::solveUpdateEq(iterationData* it){
	// Create matrices for X, Jacobian matrix DF, and constraint vector FX
	tpat_matrix oldX(it->totalFree, 1, it->X);
	tpat_matrix J(it->totalCons, it->totalFree, it->DF);
	tpat_matrix FX_mat(it->totalCons, 1, it->FX);
	
	// change sign for matrix multiplication
	FX_mat*=-1;

	// Create vector out of constraint vector FX
	gsl_vector_view b = gsl_vector_view_array(FX_mat.getDataPtr(), FX_mat.getRows());
	
	// Allocate memory for intermediate vector w
	gsl_vector *w;
	gsl_permutation *perm;
	tpat_matrix X_diff(it->totalFree, 1);
	int permSign;	// store sign (even/odd) of permutation matrix
	int status;		// status for GSL functions
	if(it->totalCons == it->totalFree){	// J is square, use regular inverse
		// J.toCSV("J.csv");
		// FX_mat.toCSV("FX.csv");
		// oldX.toCSV("X.csv");

		// Solve the system Jw = b
		w = gsl_vector_alloc(it->totalFree);
		perm = gsl_permutation_alloc(J.getRows());
		status = gsl_linalg_LU_decomp(J.getGSLMat(), perm, &permSign);
		if(status){
			printErr("tpat_correction_engine::solveUpdateEq: GSL ERR: %s\n", gsl_strerror(status));
			throw tpat_linalg_err("Unable to decompose J into L and U");
		}
		status = gsl_linalg_LU_solve(J.getGSLMat(), perm, &(b.vector), w);
		if(status){
			printErr("tpat_correction_engine::solveUpdateEq: GSL ERR: %s\n", gsl_strerror(status));
			throw tpat_linalg_err("Unable to invert J, likely singular");
		}
		// w, in this case, is X_diff
		X_diff = tpat_matrix(w, false);
	}else{
		if(it->totalCons < it->totalFree){	// Under-constrained
			// J.toCSV("J.csv");
			// FX_mat.toCSV("FX.csv");
			// oldX.toCSV("X.csv");

			// Compute Gramm matrix
			tpat_matrix G = J*trans(J);

			/* Use LU decomposition to invert the Gramm matrix and find a vector
			w. Multiplying J^T by w yields the minimum-norm solution x, where x 
			lies in the column-space of J^T, or in the orthogonal complement of
			the nullspace of J.
			Source: <http://www.math.usm.edu/lambers/mat419/lecture15.pdf>
			 */
			// Solve the system Gw = b
			w = gsl_vector_alloc(it->totalCons);
			perm = gsl_permutation_alloc(G.getRows());
			status = gsl_linalg_LU_decomp(G.getGSLMat(), perm, &permSign);
			if(status){
				printErr("tpat_correction_engine::solveUpdateEq: GSL ERR: %s\n", gsl_strerror(status));
				throw tpat_linalg_err("Unable to decompose J into L and U");
			}
			status = gsl_linalg_LU_solve(G.getGSLMat(), perm, &(b.vector), w);
			if(status){
				printErr("tpat_correction_engine::solveUpdateEq: GSL ERR: %s\n", gsl_strerror(status));
				throw tpat_linalg_err("Unable to invert G = JJ', likely singular");
			}

			// Compute the optimal x from w
			tpat_matrix W(w, false);	// create column vector
			X_diff = trans(J)*W;	//X_diff = X_new - X_old
		}else{	// Over-constrained
			// dummy allocations to avoid errors when cleaning up
			perm = gsl_permutation_alloc(J.getRows());
			w = gsl_vector_alloc(J.getRows());
			throw tpat_linalg_err("System is over constrained... No solution implemented");
		}
	}

	tpat_matrix newX = X_diff + oldX;

	// Free up memory used to invert G or J
	gsl_permutation_free(perm);
	gsl_vector_free(w);

	return newX;
}// End of solveUpdateEq() =====================================

/**
 *	@brief clean up data so that engine can be used again (or deconstructed safely)
 */
void tpat_correction_engine::cleanEngine(){
	printVerb(verbose, "Cleaning the engine...\n");
	if(createdNodesetOut){
		delete nodeset_out;
	}

	nodeset_in = 0;
	receivedNodesetIn = false;
	createdNodesetOut = false;
	isClean = true;
}//===================================
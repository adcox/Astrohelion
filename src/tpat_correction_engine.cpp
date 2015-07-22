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

/**
 *	@brief a custom data class to encapsulate data used in each iteration
 *	of the corrections process.
 *
 *	This data object can be passed to other functions, allowing us to break 
 *	the master corrections loop into smaller functions without requiring an
 *	obscene amount of arguments to be passed in.
 */
class iterationData{
	public:
		std::vector<double> X;				//!< Free-Variable Vector
		std::vector<double> FX;				//!< Constraint Function Vector
		std::vector<double> DF;				//!< Jacobian Matrix
		std::vector<double> deltaVs;		//!< nx3 vector of non-dim delta-Vs
		std::vector<double> lastState;		//!< Final state on most recent integrated arc
		std::vector<double> lastAccel;		//!< Final acceleration on most recent integrated arc
		std::vector<double> last_dqdT;		//!< Final dqdT on most recent integrated arc
		std::vector<double> primPos;		//!< Store the positions of the primaries
		std::vector<double> primVel;		//!< Store the velocities of the primaries
		std::vector<int> velConNodes;		//!< Indices of nodes that are continuous in velocity
		std::vector<int> slackAssignCon;	//!< Indices of constraints, index of entry corresponds to a slack variable
		std::vector<int> conRows;			//!< Each entry holds the row # for the constraint; i.e. 0th element holds row # for 0th constraint
		tpat_traj newSeg;				//!< Most recent integrated arc
		tpat_traj_bcr4bpr bcNewSeg;			//!< Most recent integrated arc, cast to BCR4BPR
		tpat_matrix lastSTM = tpat_matrix(6,6);	//!< Final STM on most recent integrated arc

		int numNodes = 0;			//!< Number of nodes in the entire nodeset
		int count = 0;				//!< Count of number of iterations through corrections process

		int posVelCons = 0;			//!< # position and velocity constraints
		int timeCons = 0;			//!< # time constraints
		int numSlack = 0;			//!< # slack variables
		int totalCons = 0;			//!< Total # constraints -> # rows of DF
		int totalFree = 0;			//!< Total # free var. -> # cols of DF

		int pvConCount = 0;			//!< Rolling count of position and velocity constraints applied

		bool upToDatePrimPos = false;	//!< Whether or not primPos has been updated for this node/iteration
		bool upToDatePrimVel = false;	//!< Whether or not primVel has been updated for this node/iteration
};

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
 *	@return whether or not the corrector uses variable time (as opposed
 * 	to fixed time)
 */
bool tpat_correction_engine::usesVarTime() const { return varTime; }

/**
 *	@return whether or not the corrector will be verbose
 */
bool tpat_correction_engine::isVerbose() const { return verbose; }

/**
 *	@return whether or not the algorithm will optimize the process to find an event
 */
bool tpat_correction_engine::isFindingEvent() const { return findEvent; }

/**
 *	@return the maximum number of iterations to attempt before giving up
 */
int tpat_correction_engine::getMaxIts() const { return maxIts; }

/**
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

	// Create the initial state vector
	it.X.clear();

	// Copy in the state vectors
	it.X.insert(it.X.begin(), set->getNodes()->begin(), set->getNodes()->end());

	// Append the TOF and (if applicable) node epochs
	if(varTime){
		it.X.insert(it.X.end(), set->getTOFs()->begin(), set->getTOFs()->end());
		if(sysType == tpat_sys_data::BCR4BPR_SYS){
			it.X.insert(it.X.end(), bcSet->getEpochs()->begin(), bcSet->getEpochs()->end());
		}
	}

	// Get the indices of the nodes that are continuous in velocity
	it.velConNodes = set->getVelConNodes();

	printVerb(verbose, "  Velocity-Continuous Nodes: %d\n", ((int)it.velConNodes.size()));

	// Compute number of position and velocity continuity constraints
	it.posVelCons = 3*(it.numNodes - 1) + 3*it.velConNodes.size();

	// Determine number of time constraints
	it.timeCons = 0;
	if(sysType == tpat_sys_data::BCR4BPR_SYS && varTime){
		it.timeCons = it.numNodes - 1;
	}

	// Compute number of extra consraint functions to add
	it.numSlack = 0;

	it.conRows.assign(set->getNumCons(), NAN);	// Fill with NAN
	int conRow = it.posVelCons + it.timeCons;

	bool foundDVCon = false;
	bool foundTOFCon = false;
	
	for(int c = 0; c < set->getNumCons(); c++){
		tpat_constraint con = set->getConstraint(c);

		switch(con.getType()){
			case tpat_constraint::STATE:
				it.conRows.at(c) = conRow;
				conRow += con.countConstrainedStates();
				break;
			case tpat_constraint::MATCH_ALL:
				it.conRows.at(c) = conRow;
				conRow += 6;
				break;
			case tpat_constraint::MATCH_CUST:
				it.conRows.at(c) = conRow;
				conRow += con.countConstrainedStates();
				break;
			case tpat_constraint::SP:
				it.conRows.at(c) = conRow;
				conRow += 3;
				break;
			case tpat_constraint::MAX_DIST:
			case tpat_constraint::MIN_DIST:
				it.X.push_back(1e-4);
				it.slackAssignCon.push_back(c);
				it.numSlack++;
				// do NOT break here, continue on to do stuff for DIST as well
			case tpat_constraint::DIST:
				it.conRows.at(c) = conRow;
				conRow++;
				break;
			case tpat_constraint::MAX_DELTA_V:
			case tpat_constraint::DELTA_V:
				if(!foundDVCon){
					if(((int)it.velConNodes.size()) == it.numNodes-1){
						printErr("No velocity discontinuities are allowed, but a delta-V requirement is present.\n");
						printErr("The gramm matrix will likely be singular...\n");
					}

					it.conRows.at(c) = conRow;
					conRow++;
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
				if(sysType == tpat_sys_data::CR3BP_SYS){
					it.conRows.at(c) = conRow;
					conRow++;
				}else{
					throw tpat_exception("Cannot apply JC constraint to systems other than CR3BP");
				}
				break;
			case tpat_constraint::TOF:
				if(!varTime){
					printWarn("Attempting to constraint TOF without variable time... won't work!");
				}
				if(!foundTOFCon){
					it.conRows.at(c) = conRow;
					conRow++;
				}else
					throw tpat_exception("You can only apply ONE TOF constraint");
				break;
			default: break;
		}

		if(con.getNode() < 0 || con.getNode() > it.numNodes){
			printErr("Constraint #%d applies to a non-existsant node!\n", c);
			return;
		}
	}// end of loop through constraints

	// Define values for use in corrections loop
	double err = 1000000000;

	// Determine the number of free/design variables based on the system type
	it.totalFree = 0;
	if(sysType == tpat_sys_data::CR3BP_SYS){
		it.totalFree = varTime ? (7*it.numNodes - 1 + it.numSlack) : (6*it.numNodes + it.numSlack);
	}else{
		it.totalFree = varTime ? (8*it.numNodes - 1 + it.numSlack) : (6*it.numNodes + it.numSlack);
	}

	it.totalCons = conRow;

	printVerb(verbose, "  Pos + Vel Constraints: %d\n  Time Constraints: %d\n",
		it.posVelCons, it.timeCons);
	printVerb(verbose, "  Extra Constraints: %d\n  # Free: %d\n  # Constraints: %d\n",
		it.totalCons-it.posVelCons-it.timeCons, it.totalFree, it.totalCons);
	printVerb(verbose, "  -> # Slack Variables: %d\n", it.numSlack);

	// create a simulation engine
	tpat_sys_data *sysData = set->getSysData();
	tpat_simulation_engine simEngine(*sysData);
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

	while( err > tol && it.count < maxIts){
		it.FX.clear();					// Clear vectors each iteration
		it.DF.clear();
		it.deltaVs.clear();
		it.FX.assign(it.totalCons, 0);	// Size the vectors and fill with zeros
		it.DF.assign(it.totalCons*it.totalFree, 0);
		it.deltaVs.assign(3*it.numNodes, 0);

		it.pvConCount = 0;		// reset to zero every iteration

		for(int n = 0; n < it.numNodes; n++){
			it.upToDatePrimPos = false;
			it.upToDatePrimVel = false;
			double tof = 0;
			double t0 = 0;

			// For all but the final node, integrate for the specified time
			if(n < it.numNodes-1){
				// Clear old variables to avoid any potential issues
				it.lastState.clear();
				it.lastAccel.clear();
				it.last_dqdT.clear();
				it.lastSTM = tpat_matrix::I(6);
				
				// Determine TOF and t0
				tof = varTime ? it.X[6*it.numNodes+n] : set->getTOF(n);

				if(sysType == tpat_sys_data::CR3BP_SYS)
					t0 = 0;
				else
					t0 = varTime ? it.X[7*it.numNodes-1+n] : bcSet->getEpoch(n);
			
				// Copy IC for node n from free variable vector
				double *ic = &(it.X[6*n]);
				simEngine.setRevTime(tof < 0);
				simEngine.runSim(ic, t0, tof);
				it.newSeg = simEngine.getTraj();

				// Determine the final state and final STM on the integrated segment
				it.lastState = it.newSeg.getState(-1);
				it.lastAccel = it.newSeg.getAccel(-1);
				it.lastSTM = it.newSeg.getSTM(-1);

				// Extract specific trajectory type for access to system-specific data
				if(sysType == tpat_sys_data::BCR4BPR_SYS){
					it.bcNewSeg = simEngine.getBCR4BPR_Traj();
					it.last_dqdT = it.bcNewSeg.get_dqdT(-1);
				}

				createPosVelCons(&it, sysType, n);

			}// End of if(n < it.numNodes-1)
			
			// Add extra constraints and form the partials for those constraints
			for(int c = 0; c < set->getNumCons(); c++){
				tpat_constraint con = set->getConstraint(c);
				int row0 = it.conRows[c];

				if(con.getNode() == n){
					std::vector<double> conData = con.getData();

					switch(con.getType()){
						case tpat_constraint::STATE:
							targetState(&it, con, n, row0);
							break;
						case tpat_constraint::MATCH_ALL:
							targetMatchAll(&it, con, n, row0);
							break;
						case tpat_constraint::MATCH_CUST:
							targetMatchCust(&it, con, n, row0);
							break;
						case tpat_constraint::MAX_DIST:
						case tpat_constraint::MIN_DIST:
						case tpat_constraint::DIST:
							updatePrimPos(&it, sysData, t0);
							updatePrimVel(&it, sysData, t0);
							targetDist(&it, con, sysData, n, c);
							break;
						case tpat_constraint::SP:
							if(sysType == tpat_sys_data::BCR4BPR_SYS){
								updatePrimPos(&it, sysData, t0);
								updatePrimVel(&it, sysData, t0);
								targetSP(&it, static_cast<tpat_sys_data_bcr4bpr *>(sysData), n, row0);
							}else{
								throw tpat_exception("Cannot apply SP constraint to systems other than BCR4BPR");
							}
							break;
						case tpat_constraint::JC:
							targetJC(&it, con, sysData, n, row0);
							break;
						case tpat_constraint::DELTA_V:		// handled outside this loop
						case tpat_constraint::MAX_DELTA_V:	// handled outside this loop
						case tpat_constraint::TOF:			// handled outside this loop
							break;
						default:
							throw tpat_exception("Unrecognized consraint type");
					}// End switch/case
				}// End if(con.getNode() == n)

				/* If this node allows velocity discontinuity AND this constraint is the
				delta-V constraint; */
				if(n < it.numNodes-1 && 
					std::find(it.velConNodes.begin(), it.velConNodes.end(), n+1) == it.velConNodes.end() &&
					(con.getType() == tpat_constraint::DELTA_V || 
						con.getType() == tpat_constraint::MAX_DELTA_V)){
					
					updateDeltaVCon(&it, sysData, n, row0);
				}
			}// End of loop through constraints
		}// end of loop through nodes

		// Once all nodes have been integrated, finish computing node-independent constraints
		for(int c = 0; c < set->getNumCons(); c++){

			tpat_constraint con = set->getConstraint(c);
			int row0 = it.conRows[c];
			switch(con.getType()){
				case tpat_constraint::DELTA_V:
					it.FX[row0] -= con.getData()[0];
					break;
				case tpat_constraint::MAX_DELTA_V:
				{
					// figure out which of the slack variables correspond to this constraint
					std::vector<int>::iterator slackIx = std::find(it.slackAssignCon.begin(),
						it.slackAssignCon.end(), c);

					// which column of the DF matrix the slack variable is in
					int slackCol = it.totalFree - it.numSlack + (slackIx - it.slackAssignCon.begin());

					/* FIRST ITERATION ONLY (before the targeter computes a value for beta)
					Set the slack variable such that the constraint will evaluate to zero if 
					the actual deltaV is less than the required deltaV */
					if(it.count == 0){
						it.X[slackCol] = sqrt(std::abs(con.getData()[0] - it.FX[row0]));
					}
					it.FX[row0] -= con.getData()[0] - it.X[slackCol]*it.X[slackCol];
					it.DF[it.totalFree*row0 + slackCol] = 2*it.X[slackCol];
					break;
				}
				case tpat_constraint::TOF:
				{
					// Sum all TOF for total, set partials w.r.t. integration times equal to one
					for(int i = 0; i < it.numNodes-1; i++){
						it.FX[row0] += it.X[6*it.numNodes+i];
						it.DF[it.totalFree*row0 + 6*it.numNodes+i] = 1;
					}
					
					// subtract the desired TOF from the constraint to finish its computation
					it.FX[row0] -= con.getData()[0];
					break;
				}	
				default: break;
			}// End switch(conType)
		}// END loop through constraints

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

	createOutput(&it);
}//==========================================================

/**
 *	@brief Compute position and velocity constraint values and partial derivatives
 *
 *	This function computes and stores the default position continuity constraints as well
 *	as velocity constraints for all nodes marked continuous in velocity. The delta-Vs
 *	between arc segments and node states are recorded, and the partial derivatives of each
 *	node with respect to other nodes, integration time, and epoch time are all computed
 *	and placed in the appropriate spots in the Jacobian matrix.
 *
 *	@param it contains current data for the corrections process
 *	@param sysType the type of system we're working in
 *	@param n the index of the current node
 *
 */
void tpat_correction_engine::createPosVelCons(iterationData* it, tpat_sys_data::system_t sysType, int n){
	// Add the default position constraint for this node
	it->FX[it->pvConCount+0] = it->lastState[0] - it->X[6*(n+1)+0];
	it->FX[it->pvConCount+1] = it->lastState[1] - it->X[6*(n+1)+1];
	it->FX[it->pvConCount+2] = it->lastState[2] - it->X[6*(n+1)+2];

	// Record delta-V info
	it->deltaVs[n*3+0] = it->lastState[3] - it->X[6*(n+1)+3];
	it->deltaVs[n*3+1] = it->lastState[4] - it->X[6*(n+1)+4];
	it->deltaVs[n*3+2] = it->lastState[5] - it->X[6*(n+1)+5];

	/* Add velocity constraint if applicable; we just integrated
	from node n to where (hopefully) node n+1 is. If node n+1 has
	a velocity constraint, we match the velocity at the end of newSeg
	to that node's velocity states */
	int maxRowCol = 3;	// if no vel constraint, we only need to loop through position values below
	if(std::find(it->velConNodes.begin(), it->velConNodes.end(), n+1) != it->velConNodes.end()){	// then node n+1 has a velocity constraint
		it->FX[it->pvConCount+3] = it->deltaVs[n*3+0];
		it->FX[it->pvConCount+4] = it->deltaVs[n*3+1];
		it->FX[it->pvConCount+5] = it->deltaVs[n*3+2];
		maxRowCol = 6;	// velocity constraint, use all 6 states in constraints below
	}

	// Create the partials for the default constraints. Rows correspond
	// to constraints, columns to free variables
	for(int r = 0; r < maxRowCol; r++){
		for(int c = 0; c < 6; c++){
			// put STM elements into DF matrix
			it->DF[it->totalFree*(it->pvConCount+r) + 6*n+c] = it->lastSTM.at(r,c);
			// Negative identity matrix
			if(r == c)
				it->DF[it->totalFree*(it->pvConCount+r) + 6*(n+1)+c] = -1;
		}

		// Columns of DF based on time constraints
		if(varTime){
			// Column of state derivatives: [vel; accel]
			if(r < 3)
				it->DF[it->totalFree*(it->pvConCount+r) + 6*it->numNodes+n] = it->lastState[r+3];
			else
				it->DF[it->totalFree*(it->pvConCount+r) + 6*it->numNodes+n] = it->lastAccel[r-3];

			if(sysType == tpat_sys_data::BCR4BPR_SYS){
				// Epoch dependencies
				it->DF[it->totalFree*(it->pvConCount+r) + 7*it->numNodes-1+n] = it->last_dqdT[r];
			}
		}
	}

	it->pvConCount += maxRowCol;	// update count of applied constraints
	
	/* Add time-continuity constraints if applicable; we need to match
	the epoch time of node n+1 to the sum of node n's epoch and TOF */
	if(varTime && sysType == tpat_sys_data::BCR4BPR_SYS){
		it->FX[it->posVelCons+n] = it->X[7*it->numNodes+n] - (it->X[7*it->numNodes-1+n] +
			it->X[6*it->numNodes+n]);
		it->DF[it->totalFree*(it->posVelCons+n) + 6*it->numNodes+n] = -1;
		it->DF[it->totalFree*(it->posVelCons+n) + 7*it->numNodes-1+n] = -1;
		it->DF[it->totalFree*(it->posVelCons+n) + 7*it->numNodes-1+n+1] = 1;
	}
}// End of createPosVelCons() =====================

/**
 *	@brief Compute partials and constraint functions for nodes constrained with <tt>STATE</tt>
 *
 *	@param it a pointer to the class containing all the data relevant to the corrections process
 *	@param con a copy of the constraint object
 *	@param n the index of the node that has been constrained
 *	@param row0 the index of the row this constraint begins at
 */
void tpat_correction_engine::targetState(iterationData* it, tpat_constraint con, int n, int row0){
	std::vector<double> conData = con.getData();

	// Allow user to constrain all 7 states
	
	int count = 0; 	// Count # rows since some may be skipped (NAN)
	for(int s = 0; s < ((int)con.getData().size()); s++){
		if(!isnan(conData[s])){
			if(s < 6){
				it->FX[row0+count] = it->X[6*n+s] - conData[s];
				it->DF[it->totalFree*(row0 + count) + 6*n + s] = 1;
				count++;
			}else{
				if(s == 6){
					it->FX[row0+count] = it->X[7*it->numNodes-1+n] - conData[s];
					it->DF[it->totalFree*(row0 + count) + 7*it->numNodes-1+n] = 1;
					count++;
				}else{
					throw tpat_exception("State constraints must have <= 7 elements");
				}
			}
		}
	}
}//=================================================

/**
 *	@brief Compute partials and constraint functions for nodes constrained with <tt>MATCH_ALL</tt>
 *
 *	@param it a pointer to the class containing all the data relevant to the corrections process
 *	@param con a copy of the constraint object
 *	@param n the index of the node that has been constrained
 *	@param row0 the index of the row this constraint begins at
 */
void tpat_correction_engine::targetMatchAll(iterationData* it, tpat_constraint con, int n, int row0){
	// Only allow matching 6 states, not epoch time (state 7)
	int cn = con.getData()[0];
	for(int row = 0; row < 6; row++){
		// Constrain the states of THIS node to be equal to the node 
		// with index stored in conData[0]
		it->FX[row0+row] = it->X[6*n+row] - it->X[6*cn+row];

		// Partial of this constraint wrt THIS node = I
		it->DF[it->totalFree*(row0 + row) + 6*n + row] = 1;

		// Partial of this constraint wrt other node = -I
		it->DF[it->totalFree*(row0 + row) + 6*cn+row] = -1;
	}
}//=============================================

/**
 *	@brief Compute partials and constraint functions for nodes constrained with <tt>MATCH_CUST</tt>
 *
 *	@param it a pointer to the class containing all the data relevant to the corrections process
 *	@param con a copy of the constraint object
 *	@param n the index of the node that has been constrained
 *	@param row0 the index of the row this constraint begins at
 */
void tpat_correction_engine::targetMatchCust(iterationData* it, tpat_constraint con, int n, int row0){
	std::vector<double> conData = con.getData();
	int count = 0;
	// Only allow matching 6 states, not epoch time (state 7)
	for(int s = 0; s < 6; s++){
		if(!isnan(conData[s])){
			int cn = conData[0];
			it->FX[row0 + count] = it->X[6*n+s] - it->X[6*cn+s];

			// partial of this constraint wrt THIS node = 1
			it->DF[it->totalFree*(row0 + count) + 6*n+s] = 1;

			// partial of this constraint wrt other node = -1
			it->DF[it->totalFree*(row0 + count) + 6*cn+s] = -1;

			count++;
		}
	}
}//===============================================

void tpat_correction_engine::targetJC(iterationData* it, tpat_constraint con, tpat_sys_data *sysData, int n, int row0){
	std::vector<double> conData = con.getData();
	tpat_sys_data_cr3bp *crSys = static_cast<tpat_sys_data_cr3bp *> (sysData);

	// Compute the value of Jacobi at this node
	double mu = crSys->getMu();
	double *nodeState = &(it->X[6*n]);
	double nodeJC = cr3bp_getJacobi(nodeState, mu);
	
	// temp variables to make equations more readable; compute partials w.r.t. node state
	double x = nodeState[0];
	double y = nodeState[1];
	double z = nodeState[2];
	double vx = nodeState[3];
	double vy = nodeState[4];
	double vz = nodeState[5];

	double d = sqrt((x + mu)*(x + mu) + y*y + z*z);
	double r = sqrt((x + mu - 1)*(x + mu - 1) + y*y + z*z);

	it->FX[row0] = nodeJC - conData[0];
	it->DF[it->totalFree*row0 + 6*n + 0] = -2*(x + mu)*(1 - mu)/pow(d,3) - 2*(x + mu - 1)*mu/pow(r,3) + 2*x;	//dFdx
	it->DF[it->totalFree*row0 + 6*n + 1] = -2*y*(1 - mu)/pow(d,3) - 2*y*mu/pow(r,3) - 2*y;						//dFdy
	it->DF[it->totalFree*row0 + 6*n + 2] = -2*z*(1 - mu)/pow(d,3) - 2*z*mu/pow(r,3);							//dFdz
	it->DF[it->totalFree*row0 + 6*n + 3] = -2*vx;	//dFdx_dot
	it->DF[it->totalFree*row0 + 6*n + 4] = -2*vy;	//dFdy_dot
	it->DF[it->totalFree*row0 + 6*n + 5] = -2*vz;	//dFdz_dot
}//=============================================

/**
 *	@brief Compute partials and constraint functions for nodes constrained with <tt>DIST</tt>, 
 *	<tt>MIN_DIST</tt>, or <tt>MAX_DIST</tt>
 *
 *	This function requires updated primary positions and velocities, so call updatePrimPos() and
 *	updatePrimVel() before calling this function.
 *
 *	@param it a pointer to the class containing all the data relevant to the corrections process
 *	@param con a copy of the constraint object
 *	@param sysData a pointer to the system data object for this corrections process
 *	@param n the index of the node that has been constrained
 *	@param c the index of this constraint in the constraint vector object
 */
void tpat_correction_engine::targetDist(iterationData* it, tpat_constraint con,
		tpat_sys_data* sysData, int n, int c){

	std::vector<double> conData = con.getData();
	int Pix = (int)(conData[0]);	// index of primary
	int row0 = it->conRows[c];

	// Get distance between node and primary in x, y, and z-coordinates
	double dx = it->X[6*n+0] - it->primPos[Pix*3+0];
	double dy = it->X[6*n+1] - it->primPos[Pix*3+1];
	double dz = it->X[6*n+2] - it->primPos[Pix*3+2];

	double h = sqrt(dx*dx + dy*dy + dz*dz); 	// true distance

	// Compute difference between desired distance and true distance
	it->FX[row0] = h - conData[1];

	// Partials with respect to node position states
	it->DF[it->totalFree*row0 + 6*n + 0] = dx/h;
	it->DF[it->totalFree*row0 + 6*n + 1] = dy/h;
	it->DF[it->totalFree*row0 + 6*n + 2] = dz/h;

	// Extra stuff for inequality constraints
	if(con.getType() == tpat_constraint::MIN_DIST || 
		con.getType() == tpat_constraint::MAX_DIST ){
		// figure out which of the slack variables correspond to this constraint
		std::vector<int>::iterator slackIx = std::find(it->slackAssignCon.begin(), 
			it->slackAssignCon.end(), c);

		// which column of the DF matrix the slack variable is in
		int slackCol = it->totalFree - it->numSlack + (slackIx - it->slackAssignCon.begin());
		int sign = con.getType() == tpat_constraint::MAX_DIST ? 1 : -1;

		// Subtract squared slack variable from constraint
		it->FX[row0] += sign*it->X[slackCol]*it->X[slackCol];

		// Partial with respect to slack variable
		it->DF[it->totalFree*row0 + slackCol] = sign*2*it->X[slackCol];

		// Epoch dependencies from primary positions
		if(sysData->getType() == tpat_sys_data::BCR4BPR_SYS){
			double dhdr_data[] = {-dx/h, -dy/h, -dz/h};

			tpat_matrix dhdr(1, 3, dhdr_data);
			tpat_matrix vel(3, 1, &(it->primVel[Pix*3]) );
			it->DF[it->totalFree*row0 + 7*it->numNodes - 1 + n];
		}
	}
}// End of targetDist() =========================================

/**
 *	@brief Compute partials and constraint values for nodes constrained with <tt>SP</tt>
 *
 *	This function computes three constraint values and three rows of partials for the Jacobian.
 *	Each row/function corresponds to one position state. The FX and DF matrices are updated
 *	in place by editing their values stored in <tt>it</tt>
 *
 *	This function requires updated primary positions and velocities, so call updatePrimPos() and
 *	updatePrimVel() before calling this function.
 *
 *	Note that this constraint can only be applied to the BCR4BPR as it hasn't been formulated
 *	in other systems
 *
 *	@param it the iterationData object holding the current data for the corrections process
 *	@param sysData BCR4BPR system data
 *	@param n the index of the node that is being constriant to colocate with the SP
 *	@param row0 the index of the first row for this constraint
 */
void tpat_correction_engine::targetSP(iterationData* it, tpat_sys_data_bcr4bpr* sysData, int n, int row0){

	tpat_sys_data_bcr4bpr bcSysData = *sysData;

	// Get primary positions at the specified epoch time
	tpat_matrix primPos(3, 3, &(it->primPos[0]));

	double *X = &(it->X[0]);
	tpat_matrix r(3, 1, X+6*n);		// position vector

	// Create relative position vectors between s/c and primaries
    tpat_matrix r_p1 = r - trans(primPos.getRow(0));
    tpat_matrix r_p2 = r - trans(primPos.getRow(1));
    tpat_matrix r_p3 = r - trans(primPos.getRow(2));

    double d1 = norm(r_p1);
    double d2 = norm(r_p2);
    double d3 = norm(r_p3);

	double k = bcSysData.getK();
	double mu = bcSysData.getMu();
	double nu = bcSysData.getNu();

	// Evaluate three constraint function values 
	tpat_matrix conEval = -1*(1/k - mu)*r_p1/pow(d1, 3) - (mu - nu)*r_p2/pow(d2,3) - nu*r_p3/pow(d3, 3);

	// Parials w.r.t. node position r
	double dFdq_data[9] = {0};
	dFdq_data[0] = k*k - (1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1.at(0),2)/pow(d1,5)) -
            (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2.at(0),2)/pow(d2,5)) - nu*(1/pow(d3,3) - 
            	3*pow(r_p3.at(0),2)/pow(d3,5));		//dxdx
    dFdq_data[1] = (1/k - mu)*3*r_p1.at(0)*r_p1.at(1)/pow(d1,5) + 
    		(mu - nu)*3*r_p2.at(0)*r_p2.at(1)/pow(d2,5) +
            nu*3*r_p3.at(0)*r_p3.at(1)/pow(d3,5);	//dxdy
    dFdq_data[2] = (1/k - mu)*3*r_p1.at(0)*r_p1.at(2)/pow(d1,5) +
    		(mu - nu)*3*r_p2.at(0)*r_p2.at(2)/pow(d2,5) +
            nu*3*r_p3.at(0)*r_p3.at(2)/pow(d3,5);	//dxdz
    dFdq_data[3] = dFdq_data[1];	// dydx = dxdy
    dFdq_data[4] = k*k - (1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1.at(1),2)/pow(d1,5)) -
            (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2.at(1),2)/pow(d2,5)) - 
            nu*(1/pow(d3,3) - 3*pow(r_p3.at(1),2)/pow(d3,5));	//dydy
    dFdq_data[5] = (1/k - mu)*3*r_p1.at(1)*r_p1.at(2)/pow(d1,5) +
    		(mu - nu)*3*r_p2.at(1)*r_p2.at(2)/pow(d2,5) +
            nu*3*r_p3.at(1)*r_p3.at(2)/pow(d3,5);	//dydz
    dFdq_data[6] = dFdq_data[2];	//dzdx = dxdz
    dFdq_data[7] = dFdq_data[5];	//dzdy = dydz
    dFdq_data[8] = -(1/k - mu)*(1/pow(d1,3) - 3*pow(r_p1.at(2),2)/pow(d1,5)) -
            (mu-nu)*(1/pow(d2,3) - 3*pow(r_p2.at(2),2)/pow(d2,5)) - nu*(1/pow(d3,3) - 
            3*pow(r_p3.at(2),2)/pow(d3,5));	//dzdz

    tpat_matrix dFdq(3, 3, dFdq_data);

    // Get primary velocities at the specified epoch time
    tpat_matrix primVel(3, 3, &(it->primVel[0]));

    // Compute partials of state w.r.t. primary positions; dont' compute partials
    // for P1 because its velocity is zero in the rotating frame
    double dfdr2_data[18] = {0};   double dfdr3_data[18] = {0};

    dfdr2_data[9] = -1/pow(d2,3) + 3*pow(r_p2.at(0),2)/pow(d2,5);        //dxdx2
    dfdr2_data[10] = 3*r_p2.at(0)*r_p2.at(1)/pow(d2,5);                  //dxdy2
    dfdr2_data[11] = 3*r_p2.at(0)*r_p2.at(2)/pow(d2,5);                  //dxdz2
    dfdr2_data[13] = -1/pow(d2,3) + 3*pow(r_p2.at(0),2)/pow(d2,5);       //dydy2
    dfdr2_data[14] = 3*r_p2.at(1)*r_p2.at(2)/pow(d2,5);                  //dydz2
    dfdr2_data[17] = -1/pow(d2,3) + 3*pow(r_p2.at(0),2)/pow(d2,5);       //dzdz2

    dfdr2_data[12] = dfdr2_data[10];      // Fill in symmetric matrix
    dfdr2_data[15] = dfdr2_data[11];
    dfdr2_data[16] = dfdr2_data[14];

    dfdr3_data[9] = -1/pow(d3,3) + 3*pow(r_p3.at(0),2)/pow(d3,5);        //dxdx3
    dfdr3_data[10] = 3*r_p3.at(0)*r_p3.at(1)/pow(d3,5);                  //dxdy3
    dfdr3_data[11] = 3*r_p3.at(0)*r_p3.at(2)/pow(d3,5);                  //dxdz3
    dfdr3_data[13] = -1/pow(d3,3) + 3*pow(r_p3.at(0),2)/pow(d3,5);       //dydy3
    dfdr3_data[14] = 3*r_p3.at(1)*r_p3.at(2)/pow(d3,5);                  //dydz3
    dfdr3_data[17] = -1/pow(d3,3) + 3*pow(r_p3.at(0),2)/pow(d3,5);       //dzdz3

    dfdr3_data[12] = dfdr3_data[10];      // Fill in symmetric matrix
    dfdr3_data[15] = dfdr3_data[11];
    dfdr3_data[16] = dfdr3_data[14];

    tpat_matrix dFdr2(6,3, dfdr2_data);
    tpat_matrix dFdr3(6,3, dfdr3_data);

    // scale matrices by constants
    dFdr2 *= -1*(mu - nu);
    dFdr3 *= -1*nu;

    // Compute partials of constraint function w.r.t. epoch time
    tpat_matrix dFdT = dFdr2*trans(primVel.getRow(1)) + dFdr3*trans(primVel.getRow(2));

    // Copy data into the correct vectors/matrices
    double* conEvalPtr = conEval.getDataPtr();
    double* dFdq_ptr = dFdq.getDataPtr();
    double* dFdT_ptr = dFdT.getDataPtr();

    double *FX = &(it->FX[0]);
    double *DF = &(it->DF[0]);

    std::copy(conEvalPtr, conEvalPtr+3, FX+row0);
    std::copy(dFdq_ptr, dFdq_ptr+3, DF + it->totalFree*row0 + 6*n);
    std::copy(dFdq_ptr+3, dFdq_ptr+6, DF + it->totalFree*(row0+1) + 6*n);
    std::copy(dFdq_ptr+6, dFdq_ptr+9, DF + it->totalFree*(row0+2) + 6*n);

    if(varTime){
    	std::copy(dFdT_ptr, dFdT_ptr+1, DF + it->totalFree*row0 + 7*it->numNodes-1+n);
    	std::copy(dFdT_ptr+1, dFdT_ptr+2, DF + it->totalFree*(row0+1) + 7*it->numNodes-1+n);
    	std::copy(dFdT_ptr+2, dFdT_ptr+3, DF + it->totalFree*(row0+2) + 7*it->numNodes-1+n);
    }
}// End of SP Targeting ==============================

/**
 *	@brief Update partials and constraint functions for nodes constrained with <tt>DELTA_V</tt>
 *	or <tt>MIN_DELTA_V</tt>
 *
 *	Because the delta-V constraint applies to the entire trajectory, the constraint function values
 *	and partial derivatives must be computed and upated for each node along the trajectory. You can 
 *	call this function for every node, or for every node with velocity discontinuity; the result
 *	should be the same.
 *
 *	@param it a pointer to the class containing all the data relevant to the corrections process
 *	@param sysData a pointer to the system data object for this corrections process
 *	@param n the index of the node that has been constrained
 *	@param row0 the index of the first row for this constraint
 */
void tpat_correction_engine::updateDeltaVCon(iterationData* it, tpat_sys_data* sysData, int n, int row0){
	// Compute deltaV magnitude, add to constraint function
	double dvMag = sqrt(it->deltaVs[n*3]*it->deltaVs[n*3] +
		it->deltaVs[n*3+1]*it->deltaVs[n*3+1] + 
		it->deltaVs[n*3+2]*it->deltaVs[n*3+2]);
	
	// Add the desired delta-V
	it->FX[row0] += dvMag;

	// Compute parial w.r.t. node n+1 (where velocity is discontinuous)
	double dFdq_ndf_data[] = {0, 0, 0, -1*it->deltaVs[n*3]/dvMag, 
		-1*it->deltaVs[n*3+1]/dvMag, -1*it->deltaVs[n*3+2]/dvMag};
	tpat_matrix dFdq_n2(1, 6, dFdq_ndf_data);

	// Partial w.r.t. integrated path (newSeg) from node n
	tpat_matrix dFdq_nf = -1*dFdq_n2*it->lastSTM;

	// Compute partial w.r.t. integration time n
	tpat_matrix state_dot(6, 1, &(it->lastState[3]));
	tpat_matrix dFdt_n = -1*dFdq_n2 * state_dot;

	for(int i = 0; i < 6; i++){
		it->DF[it->totalFree*row0 + 6*(n+1) + i] = dFdq_n2.at(0, i);
		it->DF[it->totalFree*row0 + 6*n + i] = dFdq_nf.at(0, i);
	}

	it->DF[it->totalFree*row0 + 6*it->numNodes+n] = dFdt_n.at(0,0);

	if(sysData->getType() == tpat_sys_data::BCR4BPR_SYS){
		// Compute partial w.r.t. epoch time n
		tpat_matrix dqdT(6, 1, it->last_dqdT);
		tpat_matrix dFdT_n = -1*dFdq_n2 * dqdT;
		it->DF[it->totalFree*row0 + 7*it->numNodes-1+n] = dFdT_n.at(0,0);
	}
}//==============================================

/**
 *	@brief Update the primary positions stored in the iterationData object
 *
 *	The iterationData object stores the positions of all system primaries, but these positions
 *	must be updated for each node (for non-autonomous systems) because nodes have different epoch
 *	times. Before doing any math, a check is performed to see if the primary position data is up
 * 	to date. If it is, no additional computations will be performed and the function will return.
 *
 *	@param it a pointer to the class containing all the data relevant to the corrections process
 *	@param sysData a pointer to the system data object for this corrections process
 *	@param t the time the current node occurs at
 */
void tpat_correction_engine::updatePrimPos(iterationData* it, tpat_sys_data* sysData, double t){
	if(!it->upToDatePrimPos){
		it->primPos.clear();

		// Compute primary positions
		switch(sysData->getType()){
			case tpat_sys_data::CR3BP_SYS:
			{
				it->primPos.assign(6,0);
				tpat_sys_data_cr3bp *crSysData = static_cast<tpat_sys_data_cr3bp *>(sysData);
				it->primPos[0] = -1*crSysData->getMu();
				it->primPos[3] = 1 - crSysData->getMu();
				break;
			}
			case tpat_sys_data::BCR4BPR_SYS:
			{
				it->primPos.assign(9,0);
				tpat_sys_data_bcr4bpr *bcSysData = static_cast<tpat_sys_data_bcr4bpr *>(sysData);
				bcr4bpr_getPrimaryPos(t, *(bcSysData), &(it->primPos[0]));
				break;
			}
			default:
				throw tpat_exception("Cannot compute primary position for system type");
		}
		it->upToDatePrimPos = true;
	}
}//================================================

/**
 *	@brief Update the primary velocities stored in the iterationData object
 *
 *	The iterationData object stores the velocities of all system primaries, but these velocities
 *	must be updated for each node (for non-autonomous systems) because nodes have different epoch
 *	times. Before doing any math, a check is performed to see if the primary velocity data is up
 * 	to date. If it is, no additional computations will be performed and the function will return.
 *
 *	@param it a pointer to the class containing all the data relevant to the corrections process
 *	@param sysData a pointer to the system data object for this corrections process
 *	@param t the time the current node occurs at
 */
void tpat_correction_engine::updatePrimVel(iterationData* it, tpat_sys_data* sysData, double t){
	if(!it->upToDatePrimVel){
		it->primVel.clear();

		switch(sysData->getType()){
			case tpat_sys_data::CR3BP_SYS:
				it->primVel.assign(6,0);
				break;
			case tpat_sys_data::BCR4BPR_SYS:
			{
				it->primVel.assign(9,0);
				tpat_sys_data_bcr4bpr *bcSysData = static_cast<tpat_sys_data_bcr4bpr *>(sysData);
				bcr4bpr_getPrimaryVel(t, *(bcSysData), &(it->primVel[0]));
				break;
			}
			default:
				throw tpat_exception("Cannot compute primary velocity for system type");
		}
		it->upToDatePrimVel = true;
	}
}//===================================================

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
		// waitForUser();

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
			if(it->count == 0){
				J.toCSV("J.csv");
				FX_mat.toCSV("FX.csv");
				oldX.toCSV("X.csv");
			}

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
 *	@brief Take the final, corrected free variable vector <tt>X</tt> and the simulation engine and create
 * 	an output nodeset. 
 *
 *	If <tt>findEvent</tt> is set to true, the
 *	output nodeset will contain extra information for the simulation engine to use. Rather than
 *	returning only the position and velocity states, the output nodeset will contain the STM 
 *	and dqdT (if non-autonomous) values for the final node; this information will be appended to 
 *	the end of the regular n x 6 vector of nodes.
 *
 *	@param it an iteration data object containing all info from the corrections process
 */
void tpat_correction_engine::createOutput(iterationData *it){

	std::vector<double> X = it->X;

	// get objects for the system type
	tpat_sys_data::system_t sysType = nodeset_in->getSysData()->getType();

	switch(sysType){
		case tpat_sys_data::CR3BP_SYS:
		{
			// Cast input nodeset to its specific type
			tpat_nodeset_cr3bp *nodeInCast = static_cast<tpat_nodeset_cr3bp *>(nodeset_in);
			// Copy that nodeset into the output guy
			nodeset_out = new tpat_nodeset_cr3bp(*nodeInCast);
			createdNodesetOut = true;
			break;
		}
		case tpat_sys_data::BCR4BPR_SYS:
		{
			// Cast input nodeset to its specific type
			tpat_nodeset_bcr4bpr *nodeInCast = static_cast<tpat_nodeset_bcr4bpr *>(nodeset_in);
			// Copy that nodeset into the output guy
			nodeset_out = new tpat_nodeset_bcr4bpr(*nodeInCast);
			createdNodesetOut = true;
			break;
		}
		default: 
			printErr("System type not supported\n");
			return;
	}

	nodeset_out->setNodeDistro(tpat_nodeset::DISTRO_NONE);
	int numNodes = nodeset_in->getNumNodes();

	// Get pointers to relevant data vectors
	std::vector<double> *nodes = nodeset_out->getNodes();
	std::vector<double> *tofs = nodeset_out->getTOFs();

	// Clear nodes and TOF and repopulate using data from output free variable vector
	nodes->clear();
	tofs->clear();
	nodes->insert(nodes->begin(), X.begin(), X.begin() + numNodes*6);
	tofs->insert(tofs->begin(), X.begin() + numNodes*6, X.begin() + numNodes*7 - 1);

	/* To avoid re-integrating in the simulation engine, we will return the entire 42 or 48-length
	state for the last node. We do this by appending the STM elements and dqdT elements to the
	end of the node array. This output nodeset should have two "nodes": the first 6 elements
	are the first node, the final 42 or 48 elements are the second node with STM and dqdT 
	information*/
	if(findEvent){
		// Append the 36 STM elements to the node vector
		nodes->insert(nodes->end(), it->lastSTM.getDataPtr(), it->lastSTM.getDataPtr()+36);

		// Append dqdT as well
		if(sysType == tpat_sys_data::BCR4BPR_SYS){
			nodes->insert(nodes->end(), it->last_dqdT.begin(), it->last_dqdT.end());
		}
	}

	// Clear the epoch vector and copy in the final, corrected epochs (if applicable)
	if(sysType == tpat_sys_data::BCR4BPR_SYS){
		tpat_nodeset_bcr4bpr *nodeOutCast = static_cast<tpat_nodeset_bcr4bpr *>(nodeset_out);
		std::vector<double> *epochs = nodeOutCast->getEpochs();
		epochs->clear();
		epochs->insert(epochs->begin(), X.begin() + 7*numNodes-1, X.begin() + 8*numNodes-1);
	}
}//=================================================================================

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
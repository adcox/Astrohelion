/**
 *	@file adtk_correction_engine.cpp
 */

/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (ADTK).
 *
 *  ADTK is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ADTK is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ATDK.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "adtk_correction_engine.hpp"

#include "adtk_ascii_output.hpp"
#include "adtk_calculations.hpp"
#include "adtk_constraint.hpp"
#include "adtk_cr3bp_nodeset.hpp"
#include "adtk_cr3bp_traj.hpp"
#include "adtk_bcr4bpr_nodeset.hpp"
#include "adtk_bcr4bpr_traj.hpp"
#include "adtk_matrix.hpp"
#include "adtk_nodeset.hpp"
#include "adtk_simulation_engine.hpp"
#include "adtk_sys_data.hpp"
#include "adtk_trajectory.hpp"
#include "adtk_utilities.hpp"

#include <cmath>
#include <gsl/gsl_linalg.h>
#include <iostream>
#include <vector>

using namespace std;

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Copy constructor - create this engine by copying the input engine
 *	@param e input correction engine
 */
adtk_correction_engine::adtk_correction_engine(const adtk_correction_engine &e){
	copyEngine(e);
}//=======================================================

/**
 *	@brief Destructor
 */
adtk_correction_engine::~adtk_correction_engine(){
	if(createdNodesetOut){
		delete nodeset_out;
	}
}//=================================================

void adtk_correction_engine::copyEngine(const adtk_correction_engine &e){
	verbose = e.verbose;
	varTime = e.varTime;
	maxIts = e.maxIts;
	tol = e.tol;
	receivedNodesetIn = e.receivedNodesetIn;
	createdNodesetOut = e.createdNodesetOut;
	findEvent = e.findEvent;
	nodeset_in = e.nodeset_in;		//POINTER, COPYING ADDRESS - passed in, so should point to same parent object

	if(createdNodesetOut){
		adtk_sys_data::system_t type = e.nodeset_out->getSysData()->getType();
		switch(type){
			case adtk_sys_data::CR3BP_SYS:
				nodeset_out = new adtk_cr3bp_nodeset (* static_cast<adtk_cr3bp_nodeset *>(e.nodeset_out));
				break;
			case adtk_sys_data::BCR4BPR_SYS:
				nodeset_out = new adtk_bcr4bpr_nodeset (*static_cast<adtk_bcr4bpr_nodeset *>(e.nodeset_out));
				break;
			default: nodeset_out = 0;
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
adtk_correction_engine& adtk_correction_engine::operator =(const adtk_correction_engine &e){
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
bool adtk_correction_engine::usesVarTime() const { return varTime; }

/**
 *	@return whether or not the corrector will be verbose
 */
bool adtk_correction_engine::isVerbose() const { return verbose; }

/**
 *	@return whether or not the algorithm will optimize the process to find an event
 */
bool adtk_correction_engine::isFindingEvent() const { return findEvent; }

/**
 *	@return the maximum number of iterations to attempt before giving up
 */
int adtk_correction_engine::getMaxIts() const { return maxIts; }

/**
 *	@return the minimum error tolerance (non-dimensional units); errors
 *	less than this value are considered negligible
 */
double adtk_correction_engine::getTol() const { return tol; }

/**
 *	@brief Retrieve the output CR3BP nodeset (after corrections). 
 *
 *	Note that this method will throw an
 *	error if the corrections process has not been run or failed to produce an output.
 *
 *	@return a CR3BP nodeset object with the corrected trajectory data stored inside
 */
adtk_cr3bp_nodeset adtk_correction_engine::getCR3BPOutput(){
	if(createdNodesetOut && receivedNodesetIn){
		if(nodeset_in->getSysData()->getType() == adtk_sys_data::CR3BP_SYS){
			// Create a copy of the nodeset, return it
			adtk_cr3bp_nodeset temp( *(static_cast<adtk_cr3bp_nodeset *>(nodeset_out)) );
			return temp;
		}else{
			printErr("Wrong system type!\n");
			throw;
		}
	}else{
		printErr("Output nodeset has not been created, cannot return CR3BP output\n");
		throw;
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
adtk_bcr4bpr_nodeset adtk_correction_engine::getBCR4BPROutput(){
	if(createdNodesetOut && receivedNodesetIn){
		if(nodeset_in->getSysData()->getType() == adtk_sys_data::BCR4BPR_SYS){
			// Create a copy of the nodeset, return it
			adtk_bcr4bpr_nodeset temp( *(static_cast<adtk_bcr4bpr_nodeset *>(nodeset_out)) );
			return temp;
		}else{
			printErr("Wrong system type!\n");
			throw;
		}
	}else{
		printErr("Output nodeset has not been created, cannot return CR3BP output\n");
		throw;
	}
}//=========================================

/**
 *	@brief Set varTime
 *	@param b whether or not the corrector should use variable time
 */
void adtk_correction_engine::setVarTime(bool b){ varTime = b; }

/**
 *	@brief Set verbosity
 *	@param b whether or not the corrector should be verbose in its outputs
 */
void adtk_correction_engine::setVerbose(bool b){ verbose = b; }

/**
 *	@brief Set maximum iterations
 *	@param i the maximum number of iterations to attempt before giving up
 */
void adtk_correction_engine::setMaxIts(int i){ maxIts = i; }

/**
 *	@brief Set the error tolerance
 *	@param d errors below this value will be considered negligible
 */
void adtk_correction_engine::setTol(double d){ tol = d; }

/**
 *	@brief Set the findEven flag
 *	@param b whether or not the algorithm will be looking for an event
 */
void adtk_correction_engine::setFindEvent(bool b){ findEvent = b; }

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Correct a CR3BP nodeset
 *	@param set a pointer to a nodeset
 */
void adtk_correction_engine::correct_cr3bp(adtk_cr3bp_nodeset* set){
	nodeset_in = set;
	receivedNodesetIn = true;
	correct(nodeset_in);
}

/**
 *	@brief Correct a BCR4BP nodeset
 *	@param set a pointer to a nodeset
 */
void adtk_correction_engine::correct_bcr4bpr(adtk_bcr4bpr_nodeset* set){
	nodeset_in = set;
	receivedNodesetIn = true;
	correct(nodeset_in);
}

/**
 *	@brief Correct a generic nodeset; equipped to handle any type
 *	@param set a pointer to a nodeset
 */
void adtk_correction_engine::correct(adtk_nodeset *set){

	// Get some basic data from the input nodeset
	int numNodes = set->getNumNodes();
	adtk_sys_data::system_t sysType = set->getSysData()->getType();
	
	printVerb(verbose, "Corrector:\n");
	printVerb(verbose, "  numNodes = %d\n", numNodes);
	printVerb(verbose, "  sysType = %s\n", set->getSysData()->getTypeStr().c_str());

	// Create specific nodeset objects using static_cast
	adtk_bcr4bpr_nodeset *bcSet;
	if(sysType == adtk_sys_data::BCR4BPR_SYS)
		bcSet = static_cast<adtk_bcr4bpr_nodeset *>(set);

	// Create the initial state vector
	vector<double> X0;

	// Copy in the state vectors
	X0.insert(X0.begin(), set->getNodes()->begin(), set->getNodes()->end());
	
	// Append the TOF and (if applicable) node epochs
	if(varTime){
		X0.insert(X0.end(), set->getTOFs()->begin(), set->getTOFs()->end());

		if(sysType == adtk_sys_data::BCR4BPR_SYS){
			X0.insert(X0.end(), bcSet->getEpochs()->begin(), bcSet->getEpochs()->end());
		}
	}

	// Get the indices of the nodes that are continuous in velocity
	vector<int> velConNodes = set->getVelConNodes();

	printVerb(verbose, "  Velocity-Continuous Nodes: %d\n", ((int)velConNodes.size()));

	// Compute number of position and velocity continuity constraints
	int posVelCons = 3*(numNodes - 1) + 3*velConNodes.size();

	// Compute number of extra consraint functions to add
	int extraCons = 0;
	int numSlack = 0;
	bool foundDVCon = false;

	// Each entry holds the index of a constraint for later reference, the index of 
	// the entry itself tells which slack variable is associated with the constraint
	vector<int> slackAssignCon;
	
	for(int c = 0; c < set->getNumCons(); c++){
		adtk_constraint con = set->getConstraint(c);

		switch(con.getType()){
			case adtk_constraint::STATE:
				extraCons += con.countConstrainedStates();
				break;
			case adtk_constraint::MATCH_ALL:
				extraCons += 6;
				break;
			case adtk_constraint::MATCH_CUST:
				extraCons += con.countConstrainedStates();
				break;
			case adtk_constraint::SP:
				// extraCons += 3;
				break;
			case adtk_constraint::MAX_DIST:
			case adtk_constraint::MIN_DIST:
				X0.push_back(1e-4);
				slackAssignCon.push_back(c);
				numSlack++;
				// do NOT break here, continue on to do stuff for DIST as well
			case adtk_constraint::DIST:
				extraCons += 1;
				break;
			case adtk_constraint::MAX_DELTA_V:
			case adtk_constraint::DELTA_V:
				if(!foundDVCon){
					if(((int)velConNodes.size()) == numNodes-1){
						printErr("No velocity discontinuities are allowed, but a delta-V requirement is present.\n");
						printErr("The gramm matrix will likely be singular...\n");
					}

					extraCons += 1;
					foundDVCon = true;

					if(con.getType() == adtk_constraint::MAX_DELTA_V){
						/* Add a slack variable to the design vector and keep track
						 * of which constraint it is assigned to; value of slack
						 * variable will be recomputed later
						 */
						X0.push_back(1e-4);
						numSlack++;
						slackAssignCon.push_back(c);
					}
				}else{
					printErr("You can only apply ONE delta-V constraint!\n");
					throw;
				}
				break;
			default: break;
		}

		if(con.getNode() < 0 || con.getNode() > numNodes){
			printErr("Constraint #%d applies to a non-existsant node!\n", c);
			return;
		}
	}// end of loop through constraints

	// Define values for use in corrections loop
	double err = 1000000000;
	int count = 0;
	vector<double> X(X0);	// Copy X0 into vector X, which will be updated each iteration

	// Determine the number of free variables and time constraints based on the system type
	int totalFree = 0;
	int timeCons = 0;
	if(sysType == adtk_sys_data::CR3BP_SYS){
		totalFree = varTime ? (7*numNodes - 1 + numSlack) : (6*numNodes + numSlack);
		timeCons = 0;	// Autonomous system, no need to constrain time
	}else{
		totalFree = varTime ? (8*numNodes - 1 + numSlack) : (6*numNodes + numSlack);
		timeCons = varTime*(numNodes - 1);
	}

	int totalCons = posVelCons + timeCons + extraCons;

	printVerb(verbose, "  Pos + Vel Constraints: %d\n  Time Constraints: %d\n",
		posVelCons, timeCons);
	printVerb(verbose, "  Extra Constraints: %d\n  # Free: %d\n  # Constraints: %d\n",
		extraCons, totalFree, totalCons);
	printVerb(verbose, "  -> # Slack Variables: %d\n", numSlack);

	// create a simulation engine
	adtk_sys_data *sysData = set->getSysData();
	adtk_simulation_engine simEngine(sysData);
	simEngine.setVerbose(verbose);
	// TODO: Should I set the tolerance or use the default?

	if(findEvent){
		// To avoid using several steps, use fixed step size and only two steps
		simEngine.setVarStepSize(false);
		simEngine.setNumSteps(2);
		simEngine.clearEvents();	// don't use crash events when searching for an event
	}

	// Create containers for matrices/vectors that I need to access inside/outside the loop
	vector<double> FX;
	vector<double> DF;
	vector<double> deltaVs;
	vector<double> lastState;
	adtk_trajectory newSeg;
	adtk_bcr4bpr_traj bcNewSeg;
	adtk_matrix lastSTM(6,6);
	vector<double> last_dqdT;

	int pvConCount;		//rolling count of pos/vel constraints applied
	int conCount;		// rolling count of extra constraint (their rows)

	while( err > tol && count < maxIts){
		FX.clear();					// Clear vectors each iteration
		DF.clear();
		deltaVs.clear();
		FX.assign(totalCons, 0);	// Size the vectors and fill with zeros
		DF.assign(totalCons*totalFree, 0);
		deltaVs.assign(3*numNodes, 0);

		// printf("DF contains %d elements\n", ((int)DF.size()));
		// printf("FX contains %d elements\n", ((int)FX.size()));
		// printf("X contains %d elements\n", ((int)X.size()));
		// printf("DeltaVs contains %d elements\n", ((int)deltaVs.size()));

		conCount = posVelCons + timeCons;	// row where extra constraints will begin
		pvConCount = 0;		// reset to zero every iteration

		for(int n = 0; n < numNodes; n++){
			double tof = 0;
			double t0 = 0;

			// For all but the final node, integrate for the specified time
			if(n < numNodes-1){
				// Clear old variables to avoid any potential issues
				lastState.clear();
				last_dqdT.clear();
				lastSTM = adtk_matrix::Identity(6);
				
				// Determine TOF and t0
				tof = varTime ? X[6*numNodes+n] : set->getTOF(n);

				if(sysType == adtk_sys_data::CR3BP_SYS)
					t0 = 0;
				else
					t0 = varTime ? X[7*numNodes-1+n] : bcSet->getEpoch(n);
			
				// Copy IC for node n from free variable vector
				double *ic = &(X[6*n]);
				simEngine.setRevTime(tof < 0);
				simEngine.runSim(ic, t0, tof);
				newSeg = simEngine.getTraj();
				
				// Extract specific trajectory type for access to system-specific data
				if(sysType == adtk_sys_data::BCR4BPR_SYS){
					bcNewSeg = simEngine.getBCR4BPRTraj();
				}

				// Determine the final state and final STM on the integrated segment
				int segEnd = newSeg.getLength() - 1;
				lastState = newSeg.getState(segEnd);
				lastSTM = newSeg.getSTM(segEnd);

				// Add the default position constraint for this node
				FX[pvConCount+0] = lastState[0] - X[6*(n+1)+0];
				FX[pvConCount+1] = lastState[1] - X[6*(n+1)+1];
				FX[pvConCount+2] = lastState[2] - X[6*(n+1)+2];

				// Record delta-V info
				deltaVs[n*3+0] = lastState[3] - X[6*(n+1)+3];
				deltaVs[n*3+1] = lastState[4] - X[6*(n+1)+4];
				deltaVs[n*3+2] = lastState[5] - X[6*(n+1)+5];

				/* Add velocity constraint if applicable; we just integrated
				from node n to where (hopefully) node n+1 is. If node n+1 has
				a velocity constraint, we match the velocity at the end of newSeg
				to that node's velocity states */
				int maxRowCol = 3;	// if no vel constraint, we only need to loop through position values below
				if(std::find(velConNodes.begin(), velConNodes.end(), n+1) != velConNodes.end()){	// then node n+1 has a velocity constraint
					FX[pvConCount+3] = deltaVs[n*3+0];
					FX[pvConCount+4] = deltaVs[n*3+1];
					FX[pvConCount+5] = deltaVs[n*3+2];
					maxRowCol = 6;	// velocity constraint, use all 6 states in constraints below
				}

				// Create the partials for the default constraints. Rows correspond
				// to constraints, columns to free variables
				for(int r = 0; r < maxRowCol; r++){
					for(int c = 0; c < 6; c++){
						// put STM elements into DF matrix
						DF[totalFree*(pvConCount+r) + 6*n+c] = lastSTM.at(r,c);
						// Negative identity matrix
						if(r == c)
							DF[totalFree*(pvConCount+r) + 6*(n+1)+c] = -1;
					}

					// Columns of DF based on time constraints
					if(varTime){
						// Column of state derivatives: [vel; accel]
						DF[totalFree*(pvConCount+r) + 6*numNodes+n] = lastState[r+3];

						if(sysType == adtk_sys_data::BCR4BPR_SYS){
							if(r == 0){	// only extract dqdT once for each node
								last_dqdT = bcNewSeg.get_dqdT(segEnd);
							}
							// Epoch dependencies
							DF[totalFree*(pvConCount+r) + 7*numNodes-1+n] = last_dqdT[r];
						}
					}
				}

				pvConCount += maxRowCol;	// update count of applied constraints
				
				/* Add time-continuity constraints if applicable; we need to match
				the epoch time of node n+1 to the sum of node n's epoch and TOF */
				if(varTime && sysType == adtk_sys_data::BCR4BPR_SYS){
					FX[posVelCons+n] = X[7*numNodes+n] - (X[7*numNodes-1+n] + X[6*numNodes+n]);
					DF[totalFree*(posVelCons+n) + 6*numNodes+n] = -1;
					DF[totalFree*(posVelCons+n) + 7*numNodes-1+n] = -1;
					DF[totalFree*(posVelCons+n) + 7*numNodes-1+n+1] = 1;
				}
			}// End of if(n < numNodes-1)
			
			// Add extra constraints and form the partials for those constraints
			for(int c = 0; c < set->getNumCons(); c++){
				adtk_constraint con = set->getConstraint(c);

				if(con.getNode() == n){
					vector<double> conData = con.getData();

					switch(con.getType()){
						case adtk_constraint::STATE:
							// Allow user to constrain all 7 states
							for(int s = 0; s < con.getNodeSize(); s++){
								if(!isnan(conData[s])){
									if(s < 6){
										FX[conCount] = X[6*n+s] - conData[s];
										DF[totalFree*conCount + 6*n + s] = 1;
										conCount++;
									}else{
										if(s == 6){
											FX[conCount] = X[7*numNodes-1+n] - conData[s];
											DF[totalFree*conCount + 7*numNodes-1+n] = 1;
											conCount++;
										}else{
											printErr("State constraints must have <= 7 elements\n");
											return;
										}
									}
								}
							}
							break;
						case adtk_constraint::MATCH_ALL:
						{
							// Only allow matching 6 states, not epoch time (state 7)
							int cn = conData[0];
							for(int row = 0; row < 6; row++){
								// Constrain the states of THIS node to be equal to the node 
								// with index stored in conData[0]
								FX[conCount+row] = X[6*n+row] - X[6*cn+row];

								// Partial of this constraint wrt THIS node = I
								DF[totalFree*(conCount + row) + 6*n + row] = 1;

								// Partial of this constraint wrt other node = -I
								DF[totalFree*(conCount + row) + 6*cn+row] = -1;
							}
							conCount += 6;
							break;
						}
						case adtk_constraint::MATCH_CUST:
						{
							// Only allow matching 6 states, not epoch time (state 7)
							for(int s = 0; s < 6; s++){
								if(!isnan(conData[s])){
									int cn = conData[0];
									FX[conCount] = X[6*n+s] - X[6*cn+s];

									// partial of this constraint wrt THIS node = 1
									DF[totalFree*conCount + 6*n+s] = 1;

									// partial of this constraint wrt other node = -1
									DF[totalFree*conCount + 6*cn+s] = -1;

									conCount++;
								}
							}
							break;
						}
						case adtk_constraint::MAX_DIST:
						case adtk_constraint::MIN_DIST:
						case adtk_constraint::DIST:
						{
							int Pix = (int)(conData[0]);	// index of primary

							// Compute primary positions
							double primPos[9] = {0};
							switch(sysType){
								case adtk_sys_data::CR3BP_SYS:
								{
									adtk_cr3bp_sys_data *crSysData = static_cast<adtk_cr3bp_sys_data *>(sysData);
									primPos[0] = -1*crSysData->getMu();
									primPos[3] = 1 - crSysData->getMu();
									break;
								}
								case adtk_sys_data::BCR4BPR_SYS:
								{
									adtk_bcr4bpr_sys_data *bcSysData = static_cast<adtk_bcr4bpr_sys_data *>(sysData);
									bcr4bpr_getPrimaryPos(t0, *(bcSysData), primPos);
									break;
								}
								default:
									printErr("Cannot compute primary position for system type %s",
										sysData->getTypeStr().c_str());
							}

							// Get distance between node and primary in x, y, and z-coordinates
							double dx = X[6*n+0] - primPos[Pix*3+0];
							double dy = X[6*n+1] - primPos[Pix*3+1];
							double dz = X[6*n+2] - primPos[Pix*3+2];

							double h = sqrt(dx*dx + dy*dy + dz*dz); 	// true distance

							// Compute difference between desired distance and true distance
							FX[conCount] = h - conData[1];

							// Partials with respect to node position states
							DF[totalFree*conCount + 6*n + 0] = dx/h;
							DF[totalFree*conCount + 6*n + 1] = dy/h;
							DF[totalFree*conCount + 6*n + 2] = dz/h;

							// Extra stuff for inequality constraints
							if(con.getType() == adtk_constraint::MIN_DIST || 
								con.getType() == adtk_constraint::MAX_DIST ){
								// figure out which of the slack variables correspond to this constraint
								vector<int>::iterator slackIx = std::find(slackAssignCon.begin(), 
									slackAssignCon.end(), c);

								// which column of the DF matrix the slack variable is in
								int slackCol = totalFree - numSlack + (slackIx - slackAssignCon.begin());
								int sign = con.getType() == adtk_constraint::MAX_DIST ? 1 : -1;

								// Subtract squared slack variable from constraint
								FX[conCount] += sign*X[slackCol]*X[slackCol];

								// Partial with respect to slack variable
								DF[totalFree*conCount + slackCol] = sign*2*X[slackCol];

								// Epoch dependencies from primary positions
								if(sysType == adtk_sys_data::BCR4BPR_SYS){
									double dhdr_data[] = {-dx/h, -dy/h, -dz/h};
									double primVel[9] = {0};

									adtk_bcr4bpr_sys_data *bcSysData = static_cast<adtk_bcr4bpr_sys_data *>(sysData);
									bcr4bpr_getPrimaryVel(t0, *(bcSysData), primVel);
									adtk_matrix dhdr(1, 3, dhdr_data);
									adtk_matrix vel(3, 1, &(primVel[Pix*3]) );
									DF[totalFree*conCount + 7*numNodes - 1 + n];
								}
							}

							conCount++;
							break;
						}// End Min, Max, Equal DIST constraints
						case adtk_constraint::SP:
							printErr("SP constraint not implemented!\n");
							return;
						case adtk_constraint::DELTA_V:	// handled outside this loop
						case adtk_constraint::MAX_DELTA_V:
							break;
						default:
							printErr("Unrecognized consraint type\n");
							return;
					}// End switch/case
				}// End if(con.getNode() == n)

				/* If this node allows velocity discontinuity AND this constraint is the
				delta-V constraint; force the delta-V constraing to be the final one 
				(last row of FX and DF)*/
				if(n < numNodes-1 && 
					std::find(velConNodes.begin(), velConNodes.end(), n+1) == velConNodes.end() &&
					(con.getType() == adtk_constraint::DELTA_V || 
						con.getType() == adtk_constraint::MAX_DELTA_V)){
					
					// Compute deltaV magnitude, add to constraint function
					double dvMag = sqrt(deltaVs[n*3]*deltaVs[n*3] + deltaVs[n*3+1]*deltaVs[n*3+1] + 
						deltaVs[n*3+2]*deltaVs[n*3+2]);
					FX[totalCons-1] += dvMag;

					// Compute parial w.r.t. node n+1 (where velocity is discontinuous)
					double dFdq_ndf_data[] = {0, 0, 0, -1*deltaVs[n*3]/dvMag, 
						-1*deltaVs[n*3+1]/dvMag, -1*deltaVs[n*3+2]/dvMag};
					adtk_matrix dFdq_n2(1, 6, dFdq_ndf_data);

					// Partial w.r.t. integrated path (newSeg) from node n
					adtk_matrix dFdq_nf = -1*dFdq_n2*lastSTM;

					// Compute partial w.r.t. integration time n
					adtk_matrix state_dot(6, 1, &(lastState[3]));
					adtk_matrix dFdt_n = -1*dFdq_n2 * state_dot;

					for(int i = 0; i < 6; i++){
						DF[totalFree*(totalCons-1) + 6*(n+1) + i] = dFdq_n2.at(0, i);
						DF[totalFree*(totalCons-1) + 6*n + i] = dFdq_nf.at(0, i);
					}

					DF[totalFree*(totalCons-1) + 6*numNodes+n] = dFdt_n.at(0,0);

					if(sysType == adtk_sys_data::BCR4BPR_SYS){
						// Compute partial w.r.t. epoch time n
						adtk_matrix dqdT(6, 1, last_dqdT);
						adtk_matrix dFdT_n = -1*dFdq_n2 * dqdT;
						DF[totalFree*(totalCons-1) + 7*numNodes-1+n] = dFdT_n.at(0,0);
					}
				}
			}// End of loop through constraints
		}// end of loop through nodes

		// Once all nodes have been integrated, finish computing Delta-V constraints
		for(int c = 0; c < set->getNumCons(); c++){

			adtk_constraint con = set->getConstraint(c);

			switch(con.getType()){
				case adtk_constraint::DELTA_V:
					FX[totalCons-1] -= con.getData()[0];
					break;
				case adtk_constraint::MAX_DELTA_V:
				{
					// figure out which of the slack variables correspond to this constraint
					vector<int>::iterator slackIx = std::find(slackAssignCon.begin(), slackAssignCon.end(), c);

					// which column of the DF matrix the slack variable is in
					int slackCol = totalFree - numSlack + (slackIx - slackAssignCon.begin());

					/* FIRST ITERATION ONLY (before the targeter computes a value for beta)
					Set the slack variable such that the constraint will evaluate to zero if 
					the actual deltaV is less than the required deltaV */
					if(count == 0)
						X[slackCol] = sqrt(abs(con.getData()[0] - FX[totalCons-1]));

					FX[totalCons-1] -= con.getData()[0] + X[slackCol]*X[slackCol];
					DF[totalFree*(totalCons - 1) + slackCol] = 2*X[slackCol];
				}
				break;
				default: break;
			}
		}

		// Create matrices for X, Jacobian matrix DF, and constraint vector FX
		adtk_matrix oldX(totalFree, 1, X);
		adtk_matrix J(totalCons, totalFree, DF);
		adtk_matrix FX_mat(totalCons, 1, FX);
		
		// change sign for matrix multiplication
		FX_mat*=-1;

		// Create vector out of constraint vector FX
		gsl_vector_view b = gsl_vector_view_array(FX_mat.getDataPtr(), FX_mat.getRows());
		
		// Allocate memory for intermediate vector w
		gsl_vector *w;
		gsl_permutation *perm;
		adtk_matrix X_diff(totalFree, 1);
		int permSign;	// store sign (even/odd) of permutation matrix
		if(totalCons == totalFree){	// J is square, use regular inverse
			// Solve the system Jw = b
			w = gsl_vector_alloc(totalFree);
			perm = gsl_permutation_alloc(J.getRows());
			gsl_linalg_LU_decomp(J.getGSLMat(), perm, &permSign);
			gsl_linalg_LU_solve(J.getGSLMat(), perm, &(b.vector), w);

			// w, in this case, is X_diff
			X_diff = adtk_matrix(w, false);
		}else{
			if(totalCons < totalFree){	// Under-constrained
				if(count == 0){
					J.toCSV("J.csv");
					FX_mat.toCSV("FX.csv");
					oldX.toCSV("X.csv");
				}
				// Compute Gramm matrix
				adtk_matrix G = J*J.trans();

				/* Use LU decomposition to invert the Gramm matrix and find a vector
				w. Multiplying J^T by w yields the minimum-norm solution x, where x 
				lies in the column-space of J^T, or in the orthogonal complement of
				the nullspace of J.
				Source: <http://www.math.usm.edu/lambers/mat419/lecture15.pdf>
				 */
				// Solve the system Gw = b
				w = gsl_vector_alloc(totalCons);
				perm = gsl_permutation_alloc(G.getRows());
				gsl_linalg_LU_decomp(G.getGSLMat(), perm, &permSign);
				gsl_linalg_LU_solve(G.getGSLMat(), perm, &(b.vector), w);

				// Compute the optimal x from w
				adtk_matrix W(w, false);	// create column vector
				X_diff = J.trans()*W;	//X_diff = X_new - X_old
			}else{	// Over-constrained
				// dummy allocations to avoid errors when cleaning up
				perm = gsl_permutation_alloc(J.getRows());
				w = gsl_vector_alloc(J.getRows());
				printVerb(verbose, "System is over constrained... No solution implemented!\n");
			}
		}

		// Solve for X_new and copy into working vector X
		adtk_matrix newX = X_diff + oldX;
		X.clear();
		X.insert(X.begin(), newX.getDataPtr(), newX.getDataPtr()+totalFree);

		// Free up memory used to invert G or J
		gsl_permutation_free(perm);
		gsl_vector_free(w);

		// Compute error; norm of constraint vector
		err = FX_mat.norm();

		count++;
		printVerbColor(verbose, YELLOW, "Iteration %02d: ||F|| = %.4e\n", count, err);
	}// end of corrections loop

	if(err > tol){
		printVerb(verbose, "Corrections process did not converge\n");
		// throw;
	}else{
		printVerb(verbose, "Corrections processes SUCCEEDED\n");
	}

	createOutput(X, simEngine);
}//==========================================================


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
 *	@param X the final, corrected free variable vector
 *	@param engine the simulation engine used to perform corrections, which must still contain
 *	the final integrated trajectory so the full state data can be saved for event function use.
 *
 */
void adtk_correction_engine::createOutput(std::vector<double> X, adtk_simulation_engine engine){

	// get objects for the system type
	adtk_sys_data::system_t sysType = nodeset_in->getSysData()->getType();

	switch(sysType){
		case adtk_sys_data::CR3BP_SYS:
		{
			// Cast input nodeset to its specific type
			adtk_cr3bp_nodeset *nodeInCast = static_cast<adtk_cr3bp_nodeset *>(nodeset_in);
			// Copy that nodeset into the output guy
			nodeset_out = new adtk_cr3bp_nodeset(*nodeInCast);
			createdNodesetOut = true;
			break;
		}
		case adtk_sys_data::BCR4BPR_SYS:
		{
			// Cast input nodeset to its specific type
			adtk_bcr4bpr_nodeset *nodeInCast = static_cast<adtk_bcr4bpr_nodeset *>(nodeset_in);
			// Copy that nodeset into the output guy
			nodeset_out = new adtk_bcr4bpr_nodeset(*nodeInCast);
			createdNodesetOut = true;
			break;
		}
		default: 
			printErr("System type not supported\n");
			return;
	}

	nodeset_out->setNodeDistro(adtk_nodeset::NONE);
	int numNodes = nodeset_in->getNumNodes();

	// Get pointers to relevant data vectors
	vector<double> *nodes = nodeset_out->getNodes();
	vector<double> *tofs = nodeset_out->getTOFs();

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

		// Get an object containing data about the final segment
		adtk_trajectory lastSeg = engine.getTraj();
		int lastSegLen = lastSeg.getLength() - 1;	// number of data points in that segment
		adtk_matrix lastSTM = lastSeg.getSTM(lastSegLen);	// final STM

		// Append the 36 STM elements to the node vector
		nodes->insert(nodes->end(), lastSTM.getDataPtr(), lastSTM.getDataPtr()+36);

		if(sysType == adtk_sys_data::BCR4BPR_SYS){
			// Grab a specific BCR4BPR trajectory object so I can append the dqdT data
			adtk_bcr4bpr_traj bcLastSeg = engine.getBCR4BPRTraj();
			vector<double> last_dqdT = bcLastSeg.get_dqdT(lastSegLen);
			nodes->insert(nodes->end(), last_dqdT.begin(), last_dqdT.end());
		}
	}

	// Clear the epoch vector and copy in the final, corrected epochs (if applicable)
	if(sysType == adtk_sys_data::BCR4BPR_SYS){
		adtk_bcr4bpr_nodeset *nodeOutCast = static_cast<adtk_bcr4bpr_nodeset *>(nodeset_out);
		vector<double> *epochs = nodeOutCast->getEpochs();
		epochs->clear();
		epochs->insert(epochs->begin(), X.begin() + 7*numNodes-1, X.begin() + 8*numNodes-1);
	}
}//===========================


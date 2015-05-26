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

#include "adtk_constraint.hpp"
#include "adtk_cr3bp_nodeset.hpp"
#include "adtk_bcr4bpr_nodeset.hpp"
#include "adtk_matrix.hpp"
#include "adtk_nodeset.hpp"
#include "adtk_sys_data.hpp"

#include <iostream>
#include <vector>

using namespace std;

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------


//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

bool adtk_correction_engine::usesVarTime() const { return varTime; }
bool adtk_correction_engine::isVerbose() const { return verbose; }
int adtk_correction_engine::getMaxIts() const { return maxIts; }
double adtk_correction_engine::getTol() const { return tol; }

void adtk_correction_engine::setVarTime(bool b){ varTime = b; }
void adtk_correction_engine::setVerbose(bool b){ verbose = b; }
void adtk_correction_engine::setMaxIts(int i){ maxIts = i; }
void adtk_correction_engine::setTol(double d){ tol = d; }

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------


void adtk_correction_engine::correct_cr3bp(adtk_cr3bp_nodeset* set){
	nodeset = set;
	correct(nodeset);
}

void adtk_correction_engine::correct_bcr4bpr(adtk_bcr4bpr_nodeset* set){
	nodeset = set;
	correct(nodeset);
}

void adtk_correction_engine::correct(adtk_nodeset *set){
	cout << "Corrector:" << endl;

	int numNodes = set->getNumNodes();
	printf("  numNodes = %d\n", numNodes);

	adtk_sys_data::system_t sysType = set->getSysData()->getType();
	printf("  sysType = %s\n", set->getSysData()->getTypeStr().c_str());

	adtk_bcr4bpr_nodeset *bcSet;
	adtk_cr3bp_nodeset *crSet;

	switch(sysType){
		case adtk_sys_data::CR3BP_SYS:
			crSet = static_cast<adtk_cr3bp_nodeset *>(set);
			break;
		case adtk_sys_data::BCR4BPR_SYS:
			bcSet = static_cast<adtk_bcr4bpr_nodeset *>(set);
			break;
		default:
			cout << "Unrecognized system type " << set->getSysData()->getTypeStr() << endl;
			throw;
	}

	// Create the initial state vector
	vector<double> X0;

	// Copy in the state vectors
	X0.insert(X0.begin(), set->getNodes()->begin(), set->getNodes()->end());
	
	if(varTime){
		// Append the TOF values
		X0.insert(X0.end(), set->getTOFs()->begin(), set->getTOFs()->end());

		// for non-autonomous systems, append epoch at each node
		if(sysType == adtk_sys_data::BCR4BPR_SYS){
			X0.insert(X0.end(), bcSet->getEpochs()->begin(), bcSet->getEpochs()->end());
		}
	}

	//TODO: determine velocity continuity; for now, assume all nodes continuous
	vector<int> velConNodes(numNodes-1, 0);
	for(int n = 1; n < numNodes; n++){
		velConNodes.at(n-1) = n;	// First node is not continuous
	}

	// Compute number of position and velocity continuity constraints
	int posVelCons = 3*(numNodes - 1) + 3*velConNodes.size();

	// Compute number of extra consraint functions to add
	int extraCons = 0;
	int numSlack = 0;
	vector<int> slackAssignCon;
	bool foundDVCon = false;
	for(int c = 0; c < set->getNumCons(); c++){
		adtk_constraint con(6);

		switch(sysType){
			case adtk_sys_data::CR3BP_SYS:
				con = crSet->getConstraint(c);
				break;
			case adtk_sys_data::BCR4BPR_SYS:
				con = bcSet->getConstraint(c);
				break;
			default: break;
		}

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
				extraCons += 3;
			case adtk_constraint::MAX_DELTA_V:
			case adtk_constraint::DELTA_V:
				if(!foundDVCon){
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
					fprintf(stderr, "You can only apply ONE delta-V constraint!\n");
					throw;
				}
				break;
			default: break;
		}

		if(con.getNode() < 0 || con.getNode() > numNodes){
			fprintf(stderr, "Constraint #%d applies to a non-existsant node!\n", c);
			throw;
		}
	}// end of loop through constraints

	// Loop and correct until we reach maxIts or the error falls below tolerance
	double err = 1000000000;
	int count = 0;
	cout << "X0 size: " << X0.size() << endl;

	// Copy X0 into vector X, which will be updated each iteration
	vector<double> X(X0);
	
	int totalFree = 0;
	int timeCons = 0;

	// Determine the number of free variables and time constraints based on the system type
	if(sysType == adtk_sys_data::CR3BP_SYS){
		totalFree = varTime ? (7*numNodes - 1 + numSlack) : (6*numNodes + numSlack);
		timeCons = 0;	// Autonomous system, no need to constrain time
	}else{
		totalFree = varTime ? (8*numNodes - 1 + numSlack) : (6*numNodes + numSlack);
		timeCons = varTime*(numNodes - 1);
	}

	int totalCons = posVelCons + timeCons + extraCons;

	if(totalCons > totalFree){
		printf("There are more constraints than free variables\n");
	}

	printf("Corrections:\n  Pos + Vel Constraints: %d\n  Time Constraints: %d\n",
		posVelCons, timeCons);
	printf("  Extra Constraints: %d\n  # Free: %d\n  # Constraints: %d\n",
		extraCons, totalFree, totalCons);

	// create a simulation engine
	adtk_simulation_engine simEngine(set->getSysData());

	while( err > tol && count < maxIts){
		// Create FX and D(F(X)) vector and matrix
		vector<double> FX(totalCons, 0);
		vector<double> DF(totalCons*totalFree, 0);

		int pvConCount = 0;		//rolling count of pos/vel constraints applied
		int conCount = posVelCons + timeCons;	// rolling count of extra constraint (their rows)

		vector<double> deltaVs(3*(numNodes-1), 0);

		//Integrate from each node, grab the final STM and state
		for(int n = 0; n < numNodes; n++){
			double tf = 0;
			double t0 = 0;

			if(n < numNodes){
				tf = varTime ? X[6*numNodes+n] : set->getTOF(n);

				if(sysType == adtk_sys_data::CR3BP_SYS)
					t0 = 0;
				else
					t0 = varTime ? X[7*numNodes-1+n] : bcSet->getEpoch(n);
			

				double *ic = X.begin()+6*(n-1);
				simEngine.setRevTime(tf < 0);
				simEngine.runSim(ic, t0, tf);

				// Save trajectory segment
				adtk_trajectory *newSeg;
				adtk_bcr4bpr_traj *bcNewSeg;
				switch(sysType){
					case adtk_sys_data::CR3BP_SYS:
						newSeg = new adtk_cr3bp_traj(simEngine.getCR3BPTraj());
						break;
					case adtk_sys_data::BCR4BPR_SYS:
						newSeg = new adtk_bcr4bpr_traj(simEngine.getBCR4BPRTraj());
						bcNewSeg = static_cast<adtk_bcr4bpr_traj *>(newSeg);
						break;
					default:
						fprintf(stderr, "Unsupport system type\n");
						throw;
				}

				int segEnd = newSeg->getLength() - 1;
				vector<double> endState = newSeg->getState(segEnd);

				// Add the default position constraint for this node
				FX[pvConCount+0] = endState[0] - X[6*n+0];
				FX[pvConCount+1] = endState[1] - X[6*n+1];
				FX[pvConCount+2] = endState[2] - X[6*n+2];

				// Record delta-V info
				deltaVs[n*3+0] = endState[3] - X[6*n+3];
				deltaVs[n*3+1] = endState[4] - X[6*n+4];
				deltaVs[n*3+2] = endState[5] - X[6*n+5];

				// Add velocity constraint if applicable; we just integrated
				// from node n to where (hopefully) node n+1 is. If node n+1 has
				// a velocity constraint, we match the velocity at the end of newSeg
				// to that node's velocity states
				int *ptr = std::find(velConNodes.begin(), velConNodes.end(), n+1);
				int maxRowCol = 3;	// if no vel constraint, we only need to loop through position values below
				if(ptr != velConNodes.end()){	// then node n+1 has a velocity constraint
					FX[pvConCount+3] = deltaVs[n*3+0];
					FX[pvConCount+4] = deltaVs[n*3+1];
					FX[pvConCount+5] = deltaVs[n*3+2];
					maxRowCol = 6;	// velocity constraint, use all 6 states in constraints below
				}

				// Create the partials for the default constraints. Rows correspond
				// to constraints, columns to free variables
				adtk_matrix lastSTM = newSeg->getSTM(segEnd);
				vector<double> last_dqdT;
				for(int r = 0; r < maxRowCol; r++){
					for(int c = 0; c < maxRowCol; c++){
						// put STM elements into DF matrix
						DF[totalFree*(pvConCount+r) + 6*n+c] = lastSTM.at(r,c);
						// Negative identity matrix
						if(r == c)
							DF[totalFree*(pvConCount+r) + 6*n+c] = -1;
					}

					// Columns of DF based on time constraints
					if(varTime){
						// Column of state derivatives: [vel; accel]
						DF[totalFree*(pvConCount+r) + 6*numNodes+n] = endState[r+3];

						if(sysType == adtk_sys_data::BCR4BPR_SYS){
							if(r == 0)
								last_dqdT = bcNewSeg->get_dqdT(segEnd);

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

				// TODO: Save the STM for this segment
			}// End of if(n < numNodes)
			
			// Add extra constraints and form the partials for those constraints
			for(int c = 0; c < set->getNumCons(); c++){
				adtk_constraint con(6);

				switch(sysType){
					case adtk_sys_data::CR3BP_SYS:
						con = crSet->getConstraint(c);
						break;
					case adtk_sys_data::BCR4BPR_SYS:
						con = bcSet->getConstraint(c);
						break;
					default: break;
				}

				if(con.getNode() == n){
					vector<double> conData = con.getData();

					switch(con.getType()){
						case adtk_constraint::STATE:
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
											fprintf(stderr, "State constraints must have <= 7 elements\n");
											throw;
										}
									}
								}
							}
							break;
						case adtk_constraint::MATCH_ALL:
						{
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
							for(int s = 0; s < 6; s++){
								if(!isnan(conData[s])){
									int cn = conData[0];
									FX[conCount] = X[6*n+s] - X[6*cn+s];

									// partial of this constraint wrt THIS node = 1
									DF[totalFree*(conCount+s) + 6*n+s] = 1;

									// partial of this constraint wrt other node = -1
									DF[totalFree*(conCount+s) + 6*cn+s] = -1;

									conCount++;
								}
							}
							break;
						}
						case adtk_constraint::SP:
							// TODO: implement
							fprintf(stderr, "SP constraint not implemented!\n");
							throw;
						default:
							fprintf(stderr, "Unrecognized consraint type\n");
							throw;
					}// End switch/case
				}// End if(con.getNode() == n)
			}

			//clean-up
			delete newSeg;
			engine.reset();
		}// end of loop through nodes

	}// end of corrections loop

}//==========================================================
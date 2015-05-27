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
#include "adtk_cr3bp_traj.hpp"
#include "adtk_bcr4bpr_nodeset.hpp"
#include "adtk_bcr4bpr_traj.hpp"
#include "adtk_matrix.hpp"
#include "adtk_nodeset.hpp"
#include "adtk_simulation_engine.hpp"
#include "adtk_sys_data.hpp"
#include "adtk_trajectory.hpp"

#include <cmath>
#include <gsl/gsl_linalg.h>
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
			return;
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

	// Get the indices of the nodes that are continuous in velocity
	vector<int> velConNodes = set->getVelConNodes();

	printf("  Velocity-Continuous Nodes: %d\n", ((int)velConNodes.size()));

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
				// extraCons += 3;
				break;
			case adtk_constraint::MAX_DELTA_V:
			case adtk_constraint::DELTA_V:
				if(!foundDVCon){
					if(((int)velConNodes.size()) == numNodes-1){
						fprintf(stderr, "No velocity discontinuities are allowed, but a delta-V requirement is present.\n");
						fprintf(stderr, "The gramm matrix will likely be singular...\n");
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
					fprintf(stderr, "You can only apply ONE delta-V constraint!\n");
					throw;
				}
				break;
			default: break;
		}

		if(con.getNode() < 0 || con.getNode() > numNodes){
			fprintf(stderr, "Constraint #%d applies to a non-existsant node!\n", c);
			return;
		}
	}// end of loop through constraints

	// Loop and correct until we reach maxIts or the error falls below tolerance
	double err = 1000000000;
	int count = 0;
	cout << "  X0 size: " << X0.size() << endl;

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

	printf("  Pos + Vel Constraints: %d\n  Time Constraints: %d\n",
		posVelCons, timeCons);
	printf("  Extra Constraints: %d\n  # Free: %d\n  # Constraints: %d\n",
		extraCons, totalFree, totalCons);

	// create a simulation engine
	adtk_sys_data *sysData = set->getSysData();
	adtk_simulation_engine simEngine(sysData);

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
		FX.clear();
		DF.clear();
		deltaVs.clear();

		// Create FX and D(F(X)) vector and matrix
		FX.assign(totalCons, 0);
		DF.assign(totalCons*totalFree, 0);
		deltaVs.assign(3*numNodes, 0);

		// printf("DF contains %d elements\n", ((int)DF.size()));
		// printf("FX contains %d elements\n", ((int)FX.size()));
		// printf("X contains %d elements\n", ((int)X.size()));
		// printf("DeltaVs contains %d elements\n", ((int)deltaVs.size()));

		conCount = posVelCons + timeCons;
		pvConCount = 0;

		//Integrate from each node, grab the final STM and state
		for(int n = 0; n < numNodes; n++){
			double tof = 0;
			double t0 = 0;

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
			

				double *ic = &(X[6*n]);
				simEngine.setRevTime(tof < 0);
				simEngine.runSim(ic, t0, tof);

				newSeg = simEngine.getTraj();
				
				if(sysType == adtk_sys_data::BCR4BPR_SYS){
					bcNewSeg = simEngine.getBCR4BPRTraj();
				}

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

				// Add velocity constraint if applicable; we just integrated
				// from node n to where (hopefully) node n+1 is. If node n+1 has
				// a velocity constraint, we match the velocity at the end of newSeg
				// to that node's velocity states
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
							if(r == 0){
								last_dqdT = bcNewSeg.get_dqdT(segEnd);
							}

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
			}// End of if(n < numNodes-1)
			
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
											fprintf(stderr, "State constraints must have <= 7 elements\n");
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
									DF[totalFree*(conCount) + 6*n+s] = 1;

									// partial of this constraint wrt other node = -1
									DF[totalFree*(conCount) + 6*cn+s] = -1;

									conCount++;
								}
							}
							break;
						}
						case adtk_constraint::SP:
							// TODO: implement
							fprintf(stderr, "SP constraint not implemented!\n");
							return;
						case adtk_constraint::DELTA_V:
						case adtk_constraint::MAX_DELTA_V:
							// handled outside this loop
							break;
						default:
							fprintf(stderr, "Unrecognized consraint type\n");
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

			//clean-up
			simEngine.reset();
		}// end of loop through nodes

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
				case adtk_constraint::DELTA_V:
					FX[totalCons-1] -= con.getData()[0];
					break;
				case adtk_constraint::MAX_DELTA_V:
				{
					// figure out which of the slack variables correspond to this constraint
					vector<int>::iterator slackIx = std::find(slackAssignCon.begin(), slackAssignCon.end(), c);

					// which column of the DF matrix the slack variable is in
					int slackCol = totalFree - numSlack + (*slackIx);

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
		gsl_vector *w = gsl_vector_alloc(totalCons);

		// Compute Gramm matrix
		adtk_matrix G = J*J.trans();

		// Save matrices to CSV files for inspection; debugging
		// if(count == 0){
		// 	J.toCSV("J.csv");
		// 	Q.toCSV("Q.csv");
		// 	FX_mat.toCSV("FX.csv");
		// 	oldX.toCSV("oldX.csv");
		// }

		/* Use LU decomposition to invert the Gramm matrix and find a vector
		w. Multiplying J^T by w yields the minimum-norm solution x, where x 
		lies in the column-space of J^T, or in the orthogonal complement of
		the nullspace of J.
		Source: <http://www.math.usm.edu/lambers/mat419/lecture15.pdf>
		 */
		// Solve the system Gw = b for w
		int s;
		gsl_permutation *perm = gsl_permutation_alloc(G.getRows());
		gsl_linalg_LU_decomp(G.getGSLMat(), perm, &s);
		gsl_linalg_LU_solve(G.getGSLMat(), perm, &(b.vector), w);

		// Compute the optimal x from w
		adtk_matrix W(w, false);	// create column vector
		adtk_matrix X_diff = J.trans()*W;	//X_diff = X_new - X_old

		// Solve for X_new and copy into working vector X
		adtk_matrix newX = X_diff + oldX;
		X.clear();
		X.insert(X.begin(), newX.getDataPtr(), newX.getDataPtr()+totalFree);

		// Free up memory used to invert G
		gsl_permutation_free(perm);
		gsl_vector_free(w);

		// Compute error; norm of constraint vector
		err = FX_mat.norm();

		// save newly computed X to CSV for debugging
		// if(count == 0){
		// 	newX.toCSV("newX.csv");
		// }

		count++;
		printf("Iteration %02d: ||F|| = %.4e\n", count, err);
	}// end of corrections loop

	if(err > tol){
		fprintf(stderr, "Corrections process did not converge\n");
		// throw;
	}else{
		printf("Corrections processes SUCCEEDED\n");
	}

	// TODO: Output data
}//==========================================================
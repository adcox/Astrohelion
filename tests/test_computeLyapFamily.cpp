#include "tpat_calculations.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_cr3bp_nodeset.hpp"
#include "tpat_cr3bp_sys_data.hpp"
#include "tpat_cr3bp_traj.hpp"
#include "tpat_constraint.hpp"
#include "tpat_linear_motion_engine.hpp"
#include "tpat_utilities.hpp"

#include <cmath>
#include <iostream>
#include <vector>

int main(void){
	int numOrbits = 800;
	int numSimple = 3;
	double step_simple = 0.0005;
	double step_fitted = 0.0005;
	int curveFitMem = 5;
	int LPt = 1;
	double r0[] = {0.001,0,0};

	tpat_cr3bp_sys_data sys("earth", "moon");

	double LPt_data[] = {0,0,0};
	cr3bp_getEquilibPt(sys, LPt, 1e-12, LPt_data);

	// Begin solving - get linear approximation at ICs
	tpat_linear_motion_engine linEngine;
	tpat_cr3bp_traj linTraj = linEngine.getCR3BPLinear(LPt,r0,
		tpat_linear_motion_engine::ELLIP,"earth","moon");

	std::vector<double> linIC = linTraj.getState(0);
	double tof = linTraj.getTime(-1);

	int orbitCount = 0;
	int iterationCount = 0;
	std::vector<double> ydot_0_guess;
	std::vector<double> trueState;

	std::vector<tpat_cr3bp_traj> members;
	const int STATE_SIZE = 8;

	while(orbitCount < numOrbits){
		ydot_0_guess.push_back(linIC[4]);
		tpat_cr3bp_traj perOrbit;
		try{
			perOrbit = cr3bp_getPeriodic(sys, linIC, tof, MIRROR_XZ);
		}catch(tpat_diverge &e){
			break;
		}catch(tpat_linalg_err &e){
			printErr("There was a linear algebra error during family continuation...");
			break;
		}
		trueState.push_back(perOrbit.getState(0)[0]);
		trueState.push_back(perOrbit.getState(0)[4]);

		printf("Orbit %03d converged!\n", ((int)members.size()));

		bool leftFamily = false;

		// Check for large changes in period to detect leaving family
		if(orbitCount > 2){
			if(perOrbit.getTime(-1) > 1.75*members[members.size()-1].getTime(-1)){
				leftFamily = true;
				printWarn("Period jumped! Left the family! Exiting...");
				break;
			}
		}

		if(!leftFamily){
			members.push_back(perOrbit);
			orbitCount++;
			
			if(orbitCount < numSimple){
				r0[0] = perOrbit.getState(0)[0] + step_simple - LPt_data[0];

				// User linear approximation if distance from Lagrange point is
				// very small; else use stupid-simple continuation
				if(r0[0]*r0[0] + r0[1]*r0[1] + r0[2]*r0[2] < 0.001){
					linTraj = linEngine.getCR3BPLinear(LPt, r0,
						tpat_linear_motion_engine::ELLIP, "earth", "moon");
					linIC = linTraj.getState(0);
				}else{
					linIC[0] = perOrbit.getState(0)[0] + step_simple;
					linIC[3] = perOrbit.getState(0)[3];
					linIC[4] = perOrbit.getState(0)[4];
					linIC[5] = perOrbit.getState(0)[5];
				}
			}else{
				int first = ((int)members.size()) - curveFitMem < 0 ? 0 : ((int)members.size()) - curveFitMem;

				// Generate a matrix of previous states to pass to the least squares fcn
				std::vector<double> prevStates;
				for(size_t n = first; n < members.size(); n++){
					// printf("Member %d has %d states\n", ((int)n), ((int)members[n].getState()->size()));
					std::vector<double> ic = members[n].getState(0);
					prevStates.insert(prevStates.end(), ic.begin(), ic.begin()+6);
					prevStates.push_back(members[n].getTime(-1));
					prevStates.push_back(members[n].getJC(0));
				}

				// Update x0
				linIC[0] = perOrbit.getState(0)[0] + step_fitted;
				printf("New x0 = %f\n", linIC[0]);

				// Use least squares to predict a new ydot0 (state 4)
				std::vector<int> depVars;
				depVars.push_back(4);
				std::vector<double> newIC = familyCont_LS(0, linIC[0], depVars, prevStates);

				linIC[4] = newIC[4];
				// printf("New vy0 = %f\n", linIC[4]);
			}
		}//End if(!leftFamily)
		tof = perOrbit.getTime(-1);
		perOrbit.saveToMat("LyapCorrected.mat");
	}//END while loop

	tpat_matrix temp1(ydot_0_guess.size(),1,ydot_0_guess);
	tpat_matrix temp2(trueState.size()/2,2,trueState);
	temp1.toCSV("ydot_0_guess.csv");
	temp2.toCSV("correctedState.csv");
}
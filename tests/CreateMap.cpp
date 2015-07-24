/**
 *	Compute a Poincare Map
 */

#include "tpat/tpat_ascii_output.hpp"
#include "tpat/tpat_constants.hpp"
#include "tpat/tpat_cr3bp_sys_data.hpp"
#include "tpat/tpat_cr3bp_traj.hpp"
#include "tpat/tpat_simulation_engine.hpp"
#include "tpat/tpat_utilities.hpp"

#include "matio.h"
#include <cmath>
#include <vector>

using namespace std;

void computeMap(vector<double> ICs, int numReturns, double C, tpat_cr3bp_sys_data sysData){
	int numICs = ICs.size()/6;
	tpat_cr3bp_sys_data tempSysData(sysData);

	printf("Beginning map creation for C = %.2f\n", C);
	printf("  Will integrate %d orbits, %d total integrations\n", numICs, 
		numICs*numReturns);

	tpat_simulation_engine engine(&tempSysData);
	engine.setVerbose(false);
	engine.addEvent(tpat_event::XZ_PLANE, 1, true);

	for(int n = 0; n < numICs; n++){
		printf("  Computing orbit %03d (%d its)\n", n, numReturns);
		vector<double> IC(ICs.begin()+n*6, ICs.begin()+(n+1)*6);
		vector<double> returns;
		returns.reserve(numReturns*6);

		// if(n == 2){ engine.setVerbose(true); }else{ engine.setVerbose(false); }

		for(int i = 0; i < numReturns; i++){
			// printColor(BLUE, "IC = [%9.4f %9.4f %9.4f %9.4f %9.4f %9.4f]\n",
			// 	IC[0], IC[1], IC[2], IC[3], IC[4], IC[5]);

			engine.runSim(&(IC[0]), 6*PI);
			tpat_cr3bp_traj traj = engine.getCR3BPTraj();
			int trajEnd = traj.getLength()-1;
			vector<double> lastState = traj.getState(trajEnd);

			if(abs(lastState[1]) > 1e-10){
				printErr("Iteration %d did not return to map! err = %e\n", i, abs(lastState[1]));
				break;
			}else if(abs(lastState[0]) > 10){
				printWarn("State returned far from barycenter... exiting!\n");
				break;
			}else{
				returns.insert(returns.end(), lastState.begin(), lastState.end()-3);
				IC = lastState;
			}
		}

		// Save data to file
		char filename[32];
		sprintf(filename, "data/JC_%03d_Orb%04d.mat", (int)floor(C*100), n);

		mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
		if(NULL == matfp){
			printErr("Error creating MAT file\n");
		}else{
			// We store data in row-major order, but the Matlab file-writing algorithm takes data
			// in column-major order, so we transpose our vector and split it into two smaller ones
			int numPoints = (int)(returns.size()/6);
			vector<double> temp(numPoints*6);

			for(int r = 0; r < numPoints; r++){
				for(int c = 0; c < 6; c++){
					if(c < 6)
						temp[c*numPoints + r] = returns[r*6 + c];
				}
			}

			size_t dims[2] = {returns.size()/6, 6};
			matvar_t *matvar = Mat_VarCreate("data", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(temp[0]), MAT_F_DONT_COPY_DATA);
			saveVar(matfp, matvar, "data", MAT_COMPRESSION_NONE);
		}

		Mat_Close(matfp);

		// waitForUser();

	}// end of loop through ICs
}//==========================================

int main(void){
	tpat_cr3bp_sys_data sys("earth", "moon");
	double mu = sys.getMu();
	double y = 0;
	double C = 3.01;

	vector<double> ICs;
	for(double x = -1.55; x < 1.55; x += 0.1){
		double d = sqrt((x + mu)*(x + mu) + y*y);
		double r = sqrt((x -1 + mu)*(x -1 + mu) + y*y);
		double y_dot_squared = -1*C + 2*(1-mu)/d + 2*mu/r + x*x + y*y;
		
		if(y_dot_squared > 0 && y_dot_squared <= 25){
			double newIC[] = {x, 0, 0, 0, sqrt(y_dot_squared), 0};
			ICs.insert(ICs.end(), newIC, newIC+6);
		}
	}

	computeMap(ICs, 10, C, sys);
}




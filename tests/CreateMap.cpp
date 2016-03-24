/**
 *	Compute a Poincare Map
 */

#include "tpat/tpat_ascii_output.hpp"
#include "tpat/tpat_constants.hpp"
#include "tpat/tpat_simulation_engine.hpp"
#include "tpat/tpat_sys_data_cr3bp.hpp"
#include "tpat/tpat_traj_cr3bp.hpp"
#include "tpat/tpat_utilities.hpp"

#include "matio.h"
#include <cmath>
#include <vector>

using namespace std;

void computeMap(vector<double> ICs, int numReturns, double C, tpat_sys_data* sysData){
	int numICs = ICs.size()/6;

	printf("Beginning map creation for C = %.2f\n", C);
	printf("  Will integrate %d orbits, %d total integrations\n", numICs, 
		numICs*numReturns);

	// Create the simulation engine
	tpat_simulation_engine engine(sysData);
	engine.setVerbose(NO_MSG);

	// Create an event that fires at the map crossing
	tpat_event mapCross(sysData, tpat_event::XZ_PLANE, 1, true);
	mapCross.setStopCount(numReturns);
	engine.addEvent(mapCross);

	// Create another event that fires if the arc strays too far from P1
	double maxDistData[] {0, 4.0};
	tpat_event maxDist(sysData, tpat_event::DIST, 1, true, maxDistData);
	engine.addEvent(maxDist);

	// Save the events (includes automatically-generated crash events)
	std::vector<tpat_event> events = engine.getEvents();


	for(int n = 0; n < numICs; n++){
		printf("  Computing orbit %03d/%03d (%d its)\n", n, numICs, numReturns);
		vector<double> IC(ICs.begin()+n*6, ICs.begin()+(n+1)*6);
		
		// Create vector to store all map returns
		vector<double> returns;
		returns.reserve((numReturns+1)*6);		// Allocate space assuming all returns are reached
		returns.insert(returns.end(), IC.begin(), IC.end());	// Save IC as the first "return" to the map

		engine.runSim(IC, numReturns*6*PI);					// Run for up to 6*pi for each map return
		tpat_traj_cr3bp traj = engine.getCR3BP_Traj();		// Extract the full trajectory
		
		// Look at all the event occurences and save the states associated with the map crossing events
		std::vector<eventRecord> allEvents = engine.getEventRecords();
		for(size_t i = 0; i < allEvents.size(); i++){
			if(events[allEvents[i].eventIx] == mapCross){
				std::vector<double> state = traj.getState(allEvents[i].stepIx);
				returns.insert(returns.end(), state.begin(), state.end());
			}
		}

		// Save data to file
		char filename[32];
		sprintf(filename, "data/JC_%03d_Orb%04d.mat", (int)floor(C*100), n);
		saveMatrixToFile(filename, "data", returns, returns.size()/6, 6);
	}// end of loop through ICs
}//============================================================

int main(void){
	tpat_sys_data_cr3bp sys("earth", "moon");
	double mu = sys.getMu();
	double y = 0;
	double C = 3.144;

	vector<double> ICs;
	for(double x = -1.55; x < 1.55; x += 0.025){
		for(double x_dot = -2.5; x_dot < 2.5; x_dot += 0.025){
			double d = sqrt((x + mu)*(x + mu) + y*y);
			double r = sqrt((x -1 + mu)*(x -1 + mu) + y*y);
			double y_dot_squared = -1*C + 2*(1-mu)/d + 2*mu/r + x*x + y*y - x_dot*x_dot;
			
			if(y_dot_squared > 0 && y_dot_squared <= 10){
				double newIC[] = {x, 0, 0, x_dot, sqrt(y_dot_squared), 0};
				ICs.insert(ICs.end(), newIC, newIC+6);
			}
		}
	}

	computeMap(ICs, 1500, C, &sys);
	char filename[128];
	sprintf(filename, "data/JC_%03d_ICs.mat", (int)std::floor(C*100));
	saveMatrixToFile(filename, "ICs", ICs, ICs.size()/6, 6);
}//============================================================




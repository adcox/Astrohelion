/**
 *  Integrate the grid of initial conditions from IceCube_ManTubeForwardProp.cpp
 *
 *	
 */

#include "tpat_constants.hpp"
#include "tpat_event.hpp"
#include "tpat_matrix.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_utilities.hpp"

#include <cstdlib>
#include <vector>

/**
 *	@brief Integrate the grid of initial conditions for the
 * 	previously computed manifold tube with the specified energy
 *
 *	@param C desired Jacobi value
 */
int main(int argc, char *argv[]){
	double C = 0;
	if(argc > 1)
		C = atof(argv[1]);
	else{
		printf("Please enter a Jacobi constant for the manifold tube!\n");
		return EXIT_SUCCESS;
	}

	char filename[64];
	sprintf(filename, "data/gridIC_C%.3f.mat", C);
	tpat_matrix gridICs = readMatrixFromMat(filename, "gridICs");
	// gridICs = trans(gridICs);
	printf("Read matrix of size %d x %d\n", gridICs.getRows(), gridICs.getCols());

	tpat_sys_data_cr3bp emSys("earth", "moon");

	int maxPeriCount = 6;
	int periCount = 0;
	double maxLunarDist = 0.4;

	tpat_simulation_engine sim(&emSys);

	double periEvtData[] = {1};					// apse at P2
	double distEvtData[] = {1, maxLunarDist};	// dist from P2

	tpat_event maxDist(&emSys, tpat_event::DIST, 1, true, distEvtData);
	tpat_event periEvent(&emSys, tpat_event::APSE, 1, true, periEvtData);

	sim.addEvent(maxDist);
	sim.addEvent(periEvent);
	// sim.setVerbose(true);

	std::vector<double> allPeri;
	allPeri.reserve(8*gridICs.getRows());

	printf("Integrating segements to find periapses...\n");
	for(int r = 0; r < 5; r++){
	// for(int r = 0; r < gridICs.getRows(); r++){
		tpat_matrix row = gridICs.getRow(r);
		printf("row has size %d x %d\n", row.getRows(), row.getCols());
		printf("ic = [%.4e %.4e %.4e %.4e %.4e %.4e]\n", row.at(0), row.at(1), row.at(2), 
			row.at(3), row.at(4), row.at(5));
		std::vector<double> IC(row.getDataPtr(), row.getDataPtr() + gridICs.getCols());
		double ellapsed = 0;

		periCount = 0;
		while(periCount < maxPeriCount){
			sim.runSim(IC, 4*PI);

			bool keepGoing = true;
			std::vector<tpat_event> endEvents = sim.getEndEvents();
			for(size_t e = 0; e < endEvents.size(); e++){
				if(endEvents[e].getType() == tpat_event::APSE){
					// found a periapse
					periCount++;
					break;
				}else if(endEvents[e].getType() == tpat_event::DIST){
					// reached maximum distance from Moon, quit
					keepGoing = false;
					break;
				}
			}

			if(keepGoing){
				tpat_traj_cr3bp traj = sim.getCR3BP_Traj();
				IC = traj.getState(-1);
				ellapsed += traj.getTime(-1);

				allPeri.push_back(r);			// Save IC index
				allPeri.push_back(ellapsed);	// Save ellapsed time since original IC
				allPeri.insert(allPeri.end(), IC.begin(), IC.end());	// Save state at periapse
			}else{
				printf("  IC %07d had %d periapses\n", r, periCount);
				break;
			}
		}

		if(periCount == maxPeriCount)
			printf("  IC %07d had %d periapses\n", r, periCount);
	}

	// Save vector of periapses to an output file
	sprintf(filename, "data/apses_C%.3f.mat", C);
	saveMatrixToFile(filename, "apses", allPeri, allPeri.size()/8, 8);

	return EXIT_SUCCESS;
}
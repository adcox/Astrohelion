#include "tpat_calculations.hpp"
#include "tpat_constants.hpp"
#include "tpat_event.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_family_cr3bp.hpp"
#include "tpat_family_member_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_utilities.hpp"

#include "matio.h"

#include <algorithm>

int main(){
	int numApses = 2;
	// int numManifolds = 200;
	int numManifolds = 1;
	
	// Specify a range of desired jacobi
	// std::vector<double> desiredJC {3.169, 3.167, 3.165, 3.163, 3.161, 3.159, 3.157, 3.155, 3.153, 3.151, 3.149, 3.147, 3.145, 3.143, 3.141};
	// std::vector<double> desiredJC {3.151, 3.149, 3.147, 3.145, 4.143, 3.141, 3.139, 3.137, 3.135, 3.133, 3.131};
	// std::vector<double> desiredJC {3.152, 3.145, 3.14, 3.135, 3.13, 3.125, 3.12, 3.115, 3.11, 3.105, 3.10, 3.005, 3.01};
	// std::vector<double> desiredJC{3.152, 3.147, 3.143, 3.14, 3.13, 3.12, 3.11, 3.10};
	std::vector<double> desiredJC{3.05};


	// Load L2 Halo family data
	tpat_family_cr3bp L2Halos("../share/families/EM_L2_NHalo.mat");
	// tpat_family_cr3bp L2Halos("../share/families/EM_L2_Lyap.mat");
	tpat_sys_data_cr3bp sysData = L2Halos.getSysData();

	std::vector<double> apses;
	for(size_t j = 0; j < desiredJC.size(); j++){
		// Get a halo at the desired Jacobi
		std::vector<tpat_family_member_cr3bp> halos = L2Halos.getMemberByJacobi(desiredJC[j]);
		printf("JC %.4f: Found %zu matches!\n", desiredJC[j], halos.size());

		if(halos.size() == 0)
			continue;

		std::vector<double> haloIC = halos[0].getIC();
		// Flip the Halo to be a southern guy
		// haloIC[2] = -haloIC[2];	// z
		// haloIC[5] = -haloIC[5];	// z-dot

		// Integrate ICs to get a full orbit
		tpat_simulation_engine sim(&sysData);
		sim.runSim(haloIC, halos[0].getTOF());
		tpat_traj_cr3bp aHalo = sim.getCR3BP_Traj();
		aHalo.saveToMat("data/ManifoldBaseHalo.mat");

		// Use a very short integration time to get the ICs for a bunch of halo manifolds
		std::vector<tpat_traj_cr3bp> manifolds = getManifolds(MAN_U_M, &aHalo, numManifolds, 1e-5);

		// double evtData[] = {1,NAN,NAN,NAN,NAN,NAN}; // Specify periapse relative to P2
		double evtData[] = {1-sysData.getMu(),NAN,NAN,NAN,NAN,NAN};	// State at x = x_moon

		sim.setRevTime(false);	// Set true for Stable manifolds
		for(size_t m = 0; m < manifolds.size(); m++){
			// char manFile[128];
			// sprintf(manFile, "data/manifold%04zu.mat", m);
			// manifolds[m].saveToMat(manFile);

			// tpat_event endEvt(&sysData, tpat_event::APSE, 1, true, evtData);	// Stop at Periapsis
			tpat_event endEvt(&sysData, tpat_event::YZ_PLANE, 0, true, evtData);	// Stop at x = x_moon
			sim.addEvent(endEvt);

			std::vector<double> IC = manifolds[m].getState(0);

			for(int a = 0; a < numApses; a++){	
				try{		
					sim.runSim(IC, 2*PI);
				}catch(tpat_diverge &e){
					break;	// Move on to the next manifold
				}

				std::vector<tpat_event> endEvents = sim.getEndEvents();
				tpat_traj_cr3bp traj = sim.getCR3BP_Traj();

				// char filename[128];
				// sprintf(filename, "data/manSeg%04zu.mat", m);
				// traj.saveToMat(filename);

				if(std::find(endEvents.begin(), endEvents.end(), endEvt) != endEvents.end()){
					printf("Succesfully found periapsis for manifold #%02zu\n", m);
					std::vector<double> lastState = traj.getState(-1);
					apses.insert(apses.end(), lastState.begin(), lastState.end());
				}else{
					printf("Did not find periapsis for manifold #%02zu\n", m);
				}

				IC = traj.getState(-1);
			}
		}
	}

	printf("Apses size: %zu\n", apses.size());

	mat_t *matfp = Mat_CreateVer("data/Apses.mat", NULL, MAT_FT_DEFAULT);
	if(NULL != matfp){
		size_t dims[2] = {6, apses.size()/6};
		matvar_t *matvar = Mat_VarCreate("Apses", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(apses[0]), MAT_F_DONT_COPY_DATA);
		saveVar(matfp, matvar, "Apses", MAT_COMPRESSION_NONE);
	}else{
		printErr("Error creating mat file\n");
	}
	Mat_Close(matfp);

	return EXIT_SUCCESS;
}
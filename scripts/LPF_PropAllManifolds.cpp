/**
 *	Attempt to replicate LPF trajectory that patches a SE manifold
 *	to a 1:3 resonant orbit
 *
 *	To compile: g++ --std=c++11 -ltpat -Wall -pedantic LPF_1-3-Res.cpp -o a.out
 */

#include "tpat_all_includes.hpp"

#include <cmath>
#include <cstdio>

int main(void){
	int numHaloNodes = 8;
	int numManNodes = 10;
	int numManifolds = 100;

	tpat_correction_engine corrector = tpat_correction_engine();
	tpat_sys_data_cr3bp seSys("sun", "earth");
	tpat_sys_data_bcr4bpr bcSys("Sun", "Earth", "Moon");
	double leaveHaloEpoch = dateToEpochTime("2016/04/18");

	// Determine geometry at t0 so that we can start at T = 0;
	bcr4bpr_orientAtEpoch(leaveHaloEpoch, &bcSys);
	tpat_simulation_engine bcEngine(&bcSys);

	std::vector<double> haloIC {0.99168151973104, 0, -0.000804558927708935, 0, -0.00983269886098604, 0};

	// Correct the halo IC to be a periodic orbit
	tpat_traj_cr3bp halo = cr3bp_getPeriodic(&seSys, haloIC, 2*PI, MIRROR_XZ, 1e-14);

	// Compute manifolds from the halo for a short amount of time to get ICs
	std::vector<tpat_traj_cr3bp> manICs = getManifolds(MAN_U_P, &halo, numManifolds, 1e-4);

	// Create a simulation engine and two events used to propagate manifolds to XZ plane
	tpat_simulation_engine engine(&seSys);
	double earthX = 1 - seSys.getMu();
	tpat_event cross_EarthPlane(&seSys, tpat_event::YZ_PLANE, 1, true, &earthX);
	tpat_event cross_Y0(&seSys, tpat_event::XZ_PLANE, 1, true);

	double bcMaxDist_data[] = {1, 2};
	tpat_event bcMaxDist(&bcSys, tpat_event::DIST, 0, true, bcMaxDist_data);
	bcEngine.addEvent(bcMaxDist);

	std::vector<tpat_traj_cr3bp> manifolds;
	for(size_t n = 0; n < manICs.size(); n++){
		printf("Manifold %03zu:\n", n);
		// First, propagate to X = X_earth to leave the halo
		engine.clearEvents();
		engine.createCrashEvents();
		engine.addEvent(cross_EarthPlane);
		engine.runSim(manICs[n].getState(0), 2*PI);
		std::vector<tpat_event> endEvents = engine.getEndEvents();

		if(std::find(endEvents.begin(), endEvents.end(), cross_EarthPlane) == endEvents.end()){
			printErr("  Manifold %zu did not encounter the Earth-X plane...\n", n);
		}else{
			// If the manifold reached that plane, propagate to the XZ plane
			tpat_traj_cr3bp o = engine.getCR3BP_Traj();
			engine.clearEvents();
			engine.createCrashEvents();
			engine.addEvent(cross_Y0);
			engine.runSim(o.getState(-1), PI);
			endEvents = engine.getEndEvents();
			if(std::find(endEvents.begin(), endEvents.end(), cross_Y0) == endEvents.end()){
				printErr("  Manifold %zu did not reach the XZ plane...\n", n);
			}else{
				// If the XZ plane has been reached, concatenate the two arcs and save
				tpat_traj_cr3bp o2 = engine.getCR3BP_Traj();
				o += o2;
				
				// Turn the manifold into a nodeset
				tpat_nodeset_cr3bp manData(o, numManNodes);

				// Create some nodes to represent the halo
				tpat_nodeset_cr3bp haloData(o.getState(0), &seSys, -PI, numHaloNodes);
				double haloTOF = std::abs(haloData.getTotalTOF());

				// Reverse order, tack the manifold on to the end
				haloData.reverseOrder();
				haloData.deleteNode(-1);
				haloData += manData;

				// Transform to BCR4BP system
				double t0 = (0 - haloTOF*seSys.getCharT())/bcSys.getCharT();	// t = 0 occurs at the beginning of the manifold, so just subtract the halo TOF
				tpat_nodeset_bcr4bp bcNodes = bcr4bpr_SE2SEM(haloData, &bcSys, t0);

				// Constrain the final node to be on the XZ_Plane
				double fixY0_Data[] = {NAN, 0, NAN, NAN, NAN, NAN};
				tpat_constraint fixY0(tpat_constraint::STATE, bcNodes.getNumNodes()-1, fixY0_Data, 6);
				bcNodes.addConstraint(fixY0);

				// Correct to be continuous
				try{
					corrector.setTol(1e-10);
					corrector.multShoot(&bcNodes);

					// Now propagate forward for a while from the final state
					tpat_nodeset_bcr4bp bcCorrected = corrector.getBCR4BPR_Output();
					bcEngine.runSim(bcCorrected.getState(-1), bcCorrected.getEpoch(-1), 360*24*3600/bcSys.getCharT());
					tpat_traj_bcr4bp natArc = bcEngine.getBCR4BPR_Traj();
					tpat_traj_bcr4bp fullTraj = tpat_traj_bcr4bp::fromNodeset(bcCorrected);
					fullTraj += natArc;

					// Save the trajectory to file
					char name[32];
					sprintf(name, "data/LPF_4B_NaturalManifolds/Traj%03zu.mat", n);
					fullTraj.saveToMat(name);
				}catch(tpat_diverge &e){
					printErr("  Unable to correct manifold %03zu to be continuous in BC4BP... darn!\n", n);
				}
			}
		}
	}


	// // Generate nodesets for the halo manifold and the halo itself, both in reverse time
	// tpat_nodeset_cr3bp manData(IC, &seSys, tf, numManNodes);
	// tpat_nodeset_cr3bp haloData(manData.getNode(-1).getPosVelState(), &seSys, -PI, numHaloNodes);

	// // Concatenate nodesets;
	// manData.deleteNode(-1);
	// tpat_nodeset_cr3bp manToHaloData = manData + haloData;

	// // Change to BCR4BPR System
	// tpat_nodeset_bcr4bp manToHalo_bcNodes = bcr4bpr_SE2SEM(manToHaloData, &bcSys, t0);
	// manToHalo_bcNodes.reverseOrder();

	// // Make all nodes continuous in velocity
	// std::vector<int> notContinuous = std::vector<int>();
	// manToHalo_bcNodes.setVelConNodes_allBut(notContinuous);

	// // manToHalo_bcNodes.print();

	// // Correct to be continuous
	// tpat_correction_engine corrector = tpat_correction_engine();
	// corrector.setVerbose(SOME_MSG);
	// corrector.setTol(1e-10);
	// corrector.setScaleVars(false);
	// corrector.multShoot(&manToHalo_bcNodes);
	// tpat_nodeset_bcr4bp contNodes = corrector.getBCR4BPR_Output();

	// // Propagate past final node until a negative XZ plane crossing (hopefully near SP)
	// double maxTOF = 5.0/bcSys.getK();
	// engine.addEvent(tpat_event::XZ_PLANE, -1, true);
	// engine.runSim(contNodes.getNode(-1).getPosVelState(), t0, maxTOF);
	// tpat_traj_bcr4bp newArc = engine.getBCR4BPR_Traj();
	
	// // Create a nodeset from the arc segment
	// double newArcTOF = newArc.getTime(-1);
	// double numNewArcNodes = 5;
	// tpat_nodeset_bcr4bp evolveSegNodes(contNodes.getNode(-1).getPosVelState(), &bcSys, t0, newArcTOF, numNewArcNodes);

	// // Concatenate into one large nodeset
	// contNodes.deleteNode(-1);
	// tpat_nodeset_bcr4bp fullNodeset = contNodes + evolveSegNodes;
}
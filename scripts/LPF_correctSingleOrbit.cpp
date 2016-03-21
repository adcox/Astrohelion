/**
 *	Generate a bunch of manifold arcs and correct them to be natural in the BC4BP
 *
 *	To compile: g++ --std=c++11 -ltpat -Wall -pedantic LPF_correctSingleOrbit.cpp -o a.out
 */

#include "tpat_all_includes.hpp"

#include <cmath>
#include <cstdio>

int main(void){

	int orbNum = 83;

	double extraTOF = 360*24*3600;
	switch(orbNum){
		case 8: extraTOF = 180*24*3600; break;	// Enough time to fly by the SP multiple times
		case 66: extraTOF = 367*24*3600; break;
		case 81: extraTOF = 150*24*3600; break;
	}

	// ************************************************************
	// DO NOT EDIT THE CODE BELOW HERE UNTIL THE NEXT COMMENT BLOCK
	// ************************************************************
	int numHaloNodes = 8;
	int numManNodes = 10;
	int numManifolds = 100;

	tpat_correction_engine corrector = tpat_correction_engine();
	tpat_sys_data_cr3bp seSys("sun", "earth");
	tpat_sys_data_cr3bp emSys("earth", "moon");
	tpat_sys_data_bcr4bpr bcSys("Sun", "Earth", "Moon");
	double leaveHaloEpoch = dateToEpochTime("2016/04/18");
	tpat_nodeset_bcr4bp finalNodes(&bcSys);

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

	size_t n = (size_t)orbNum;
	
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
				finalNodes = bcCorrected;
				bcEngine.runSim(bcCorrected.getState(-1), bcCorrected.getEpoch(-1), extraTOF/bcSys.getCharT());
				tpat_traj_bcr4bp natArc = bcEngine.getBCR4BPR_Traj();
				tpat_traj_bcr4bp fullTraj = tpat_traj_bcr4bp::fromNodeset(bcCorrected);
				
				fullTraj.saveToMat("fullTraj.mat");
				natArc.saveToMat("natArc.mat");
				
				fullTraj += natArc;
				
				// Save the trajectory to file
				char name[32];
				sprintf(name, "data/Traj%03zu.mat", n);
				fullTraj.saveToMat(name);

				tpat_traj_cr3bp fullTraj_SE = bcr4bpr_SEM2SE(fullTraj, &seSys);
				sprintf(name, "data/Traj%03zu_SE.mat", n);
				fullTraj_SE.saveToMat(name);

				tpat_traj_cr3bp fullTraj_EM = cr3bp_SE2EM(fullTraj_SE, &emSys, bcSys.getTheta0(), bcSys.getPhi0(), bcSys.getGamma());
				sprintf(name, "data/Traj%03zu_EM.mat", n);
				fullTraj_EM.saveToMat(name);

				waitForUser();
			}catch(tpat_diverge &e){
				printErr("  Unable to correct manifold %03zu to be continuous in BC4BP... darn!\n", n);
			}
		}
	}
	
	// ************************************************************
	// YOU CAN EDIT BELOW THIS BLOCK!
	// ************************************************************
	switch(orbNum){
		case 8:
		{
			int nodeBeforeFlyby = numHaloNodes+numManNodes-8;	// -3 is the node before the lunar flyby
			// finalNodes.print();
			
			// Apse event doesn't appear to work for BC4BP, so use a hard-coded plane constraint for now
			double bcFlybyPlaneX = 378770/bcSys.getCharL();
			// tpat_event flybyEvent(&bcSys, tpat_event::YZ_PLANE, -1, true, &bcFlybyPlaneX);

			double moonIx = 2;
			tpat_event flybyEvent(&bcSys, tpat_event::APSE, 1, true, &moonIx);

			bcEngine.addEvent(flybyEvent);
			// bcEngine.setVerbose(ALL_MSG);
			bcEngine.runSim(finalNodes.getState(nodeBeforeFlyby), finalNodes.getEpoch(nodeBeforeFlyby),
				240*24*3600/bcSys.getCharT());

			tpat_traj_bcr4bp seg = bcEngine.getBCR4BPR_Traj();

			tpat_traj_cr3bp seg_SE = bcr4bpr_SEM2SE(seg, &seSys);
			tpat_traj_cr3bp seg_EM = cr3bp_SE2EM(seg_SE, &emSys, bcSys.getTheta0(), bcSys.getPhi0(), bcSys.getGamma());
			seg_EM.saveToMat("seg_EM.mat");
			seg.saveToMat("seg.mat");

			waitForUser();

			// Create nodeset for this first leg up to the lunar flyby
			tpat_nodeset_bcr4bp nodesToFlyby(seg.getState(0), &bcSys, seg.getTime(0),
				seg.getTime(-1) - seg.getTime(0), numManNodes);

			// Create a nodeset for post-flyby
			bcEngine.clearEvents();
			bcEngine.createCrashEvents();
			bcEngine.runSim(nodesToFlyby.getState(-1), nodesToFlyby.getEpoch(-1), 180*24*3600/bcSys.getCharT());
			seg = bcEngine.getBCR4BPR_Traj();

			tpat_nodeset_bcr4bp nodesPostFlyby(nodesToFlyby.getState(-1), &bcSys, nodesToFlyby.getEpoch(-1),
				180*24*3600/bcSys.getCharT(), 50);

			// Concatenate nodesets
			tpat_nodeset_bcr4bp fullSet = nodesToFlyby;
			fullSet.deleteNode(-1);
			fullSet += nodesPostFlyby;

			std::vector<int> dvNodes {8, 11, 20, 30};
			fullSet.allowDV_at(dvNodes);

			// Fix initial manifold state
			std::vector<double> initState = fullSet.getState(0);
			initState.push_back(fullSet.getEpoch(0));
			tpat_constraint fixState(tpat_constraint::STATE, 0, initState);

			// Planar didn't work
			// double planarState[] = {NAN, NAN, 0, NAN, NAN, NAN};
			// tpat_constraint makePlanar(tpat_constraint::STATE, 10, planarState, 6);

			// Encounter Saddle Point
			tpat_constraint spPass(tpat_constraint::SP, 54, &bcFlybyPlaneX, 1);
			
			// Fix maximum allowable delta-V
			double maxDVAmt = 250/1000*bcSys.getCharT()/bcSys.getCharL();
			tpat_constraint maxDV(tpat_constraint::MAX_DELTA_V, 0, &maxDVAmt, 1);


			fullSet.addConstraint(fixState);
			fullSet.addConstraint(spPass);
			fullSet.addConstraint(maxDV);

			// corrector.setMaxIts(100);
			// corrector.multShoot(&fullSet);
			// tpat_nodeset_bcr4bp cNodes = corrector.getBCR4BPR_Output();

			// bcEngine.clearEvents();
			// bcEngine.createCrashEvents();
			// bcEngine.addEvent(bcPerilune);
			// bcEngine.runSim(seg.getState(-1), seg.getTime(-1), 5*PI);

			// seg = bcEngine.getBCR4BPR_Traj();
			fullSet.saveToMat("fullSet.mat");
			// cNodes.saveToMat("cNodes.mat");

			break;
		}
	}
}
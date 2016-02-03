/**
 *	Generate a bunch of manifold arcs and correct them to be natural in the BC4BP
 *
 *	To compile: g++ --std=c++11 -ltpat -Wall -pedantic LPF_PropAllManifolds.cpp -o a.out
 */

#include "tpat_all_includes.hpp"

#include <cmath>
#include <cstdio>

int main(void){

	// Load trajectory from file
	int trajIx = 86;
	char filename[128];
	sprintf(filename, "data/LPF_4B_NaturalManifolds/Traj%03d_SEM.mat", trajIx);

	printf("Loading data from %s\n", filename);

	tpat_sys_data_bcr4bpr bcSys(filename);		// First, read system data from file
	tpat_traj_bcr4bp bcTraj(&bcSys);			// Create trajectory with the system data
	bcTraj.readFromMat(filename);				// Read trajectory data in from the file

	double tof = bcTraj.getTime(-1) - bcTraj.getTime(0);

	// Create simulation engine, events to stop trajectory near SP for easy node placement
	tpat_simulation_engine sim(&bcSys);
	double sp_xCoord = -0.173005958645270;		// Sun-Earth (no Moon) SP x-coord 
	tpat_event sp_yzPlane(&bcSys, tpat_event::YZ_PLANE, 1, true, &sp_xCoord);
	tpat_event sp_xzPlane(&bcSys, tpat_event::XZ_PLANE, 0, true);
	
	double BodyIx = 2;
	tpat_event apseCon(&bcSys, tpat_event::APSE, 0, true, &BodyIx);

	// First propagate to the SP YZ plane
	sim.addEvent(sp_yzPlane);
	sim.runSim(bcTraj.getState(0), bcTraj.getTime(0), tof);
	tpat_traj_bcr4bp arc = sim.getBCR4BPR_Traj();

	// Next, propagate to XZ Plane
	sim.clearEvents();
	sim.addEvent(sp_xzPlane);
	sim.runSim(arc.getState(-1), arc.getTime(-1), tof);
	tpat_traj_bcr4bp arc2 = sim.getBCR4BPR_Traj();

	arc += arc2;
	arc.saveToMat("arc.mat");
	tpat_nodeset_bcr4bp manNodes(arc.getState(0), &bcSys, arc.getTime(0), arc.getTOF(), 20);

	// Now propagate from this point to create a nodeset
	tpat_nodeset_bcr4bp capturedNodes(&bcSys);	
	
	std::vector<double> q0 = arc2.getState(-1);
	double T0 = arc2.getTime(-1);
	tof = 400*24*3600/bcSys.getCharT();

	for(int i = 0; i < 5; i++){
		// Integrate to XZ Plane crossing
		sim.runSim(q0, T0, tof);
		tpat_traj_bcr4bp tempTraj = sim.getBCR4BPR_Traj();
		tpat_nodeset_bcr4bp tempNodes(q0, &bcSys, T0, tempTraj.getTOF(), 4);
		
		if(i == 0){
			capturedNodes = tempNodes;
		}else{
			capturedNodes.deleteNode(-1);
			capturedNodes += tempNodes;
		}

		q0 = capturedNodes.getState(-1);
		T0 = capturedNodes.getEpoch(-1);
	}

	capturedNodes.saveToMat("CapturedNodes.mat");

	// Determine which nodes are nearest the saddle point
	for(int n = 0; n < capturedNodes.getNumNodes(); n++){
		std::vector<double> nodeState = capturedNodes.getState(n);
		Eigen::Vector3d scLoc = Eigen::Map<Eigen::Vector3d>(&(nodeState[0]));
		Eigen::Vector3d spLoc = bcr4bpr_getSPLoc(&bcSys, capturedNodes.getEpoch(n));
		Eigen::Vector3d diff = spLoc - scLoc;
		printf("Node %02d: SP Pass Distance = %.2f km\n", n, diff.norm()*bcSys.getCharL());
	}

	// waitForUser();
	// Manually ID the two closest passes
	int spNode1 = 9;
	int spNode2 = 15;

	tpat_constraint sp1Con(tpat_constraint::SP, spNode1, NULL, 0);
	tpat_constraint sp2Con(tpat_constraint::SP, spNode2, NULL, 0);

	std::vector<double> firstNode = capturedNodes.getState(0);
	firstNode.push_back(capturedNodes.getEpoch(0));
	tpat_constraint fixFirst(tpat_constraint::STATE, 0, firstNode);

	double maxDV = 1000/1000/bcSys.getCharL()*bcSys.getCharT();
	tpat_constraint maxDVCon(tpat_constraint::MAX_DELTA_V, 0, &maxDV, 1);

	capturedNodes.addConstraint(sp1Con);
	capturedNodes.addConstraint(sp2Con);
	capturedNodes.addConstraint(fixFirst);
	// capturedNodes.addConstraint(maxDVCon);

	// std::vector<int> dvNodes {5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65};
	std::vector<int> dvNodes {1, 7, 12};

	capturedNodes.allowDV_at(dvNodes);

	manNodes.deleteNode(-1);
	manNodes += capturedNodes;
	manNodes.print();

	waitForUser();

	tpat_correction_engine corrector;
	corrector.setTol(1e-11);
	corrector.setMaxIts(50);
	iterationData itData = corrector.multShoot(&capturedNodes);

	double totalDV = getTotalDV(&itData);
	double charV = bcSys.getCharL()/bcSys.getCharT();
	printf("Corrected trajectory with Delta-V = %.4f m/s\n", 1000*totalDV*charV);

	tpat_nodeset_bcr4bp corNodes = corrector.getBCR4BPR_Output();
	corNodes.saveToMat("CorrectedNodes.mat");

	// // Create SP constraints and add to nodeset
	// tpat_constraint sp1Con(tpat_constraint::SP, numNodes1-1, NULL, 0);
	// tpat_constraint sp2Con(tpat_constraint::SP, allNodes.getNumNodes()-1, NULL, 0);
	// allNodes.addConstraint(sp1Con);
	// allNodes.addConstraint(sp2Con);

	// // Allow velocity discontinuities
	// std::vector<int> dvNodes {numNodes1 - 2, allNodes.getNumNodes()-5};
	// allNodes.allowDV_at(dvNodes);

	// // Correct to encounter SP
	// tpat_correction_engine corrector;
	// corrector.setTol(1e-11);
	// iterationData itData = corrector.multShoot(&allNodes);
	// double totalDV = getTotalDV(&itData);
	// double charV = bcSys.getCharL()/bcSys.getCharT();
	// printf("Corrected trajectory with Delta-V = %.4f m/s\n", 1000*totalDV*charV);

	// tpat_nodeset_bcr4bp corNodes = corrector.getBCR4BPR_Output();
	// sprintf(filename, "data/LPF_4B_NaturalManifolds/Traj%03d_corrected.mat", trajIx);
	// corNodes.saveToMat(filename);

	// // Attempt to step DV down by iteratively correcting with a maxDV constraint
	// printColor(BOLDBLACK, "\nBeginning Delta-V Step-Down Loop\n");

	// // dvNodes.push_back(numNodes1);
	// // allNodes.allowDV_at(dvNodes);

	// double approxDV_dim = std::floor(totalDV*1000*charV)/1000;	// round to nearest 1 m/s
	// tpat_constraint maxDVCon(tpat_constraint::MAX_DELTA_V, 0, &approxDV_dim, 1);	// wrong dimensions on DV, updated in loop

	// for(double dvMax = approxDV_dim; dvMax > 0.01; dvMax -= 0.001){
	// 	double dvMax_nonDim = dvMax/charV;
		

	// 	maxDVCon.setData(&dvMax_nonDim, 1);
	// 	corNodes.clearConstraints();
	// 	corNodes.addConstraint(sp1Con);
	// 	corNodes.addConstraint(sp2Con);
	// 	corNodes.addConstraint(maxDVCon);
	// 	// corNodes.addConstraint(fixFirst);

	// 	try{
	// 		corrector.setMaxIts(40);
	// 		// corrector.setVerbose(NO_MSG);
	// 		// convergedNodes.print();
	// 		printf("Attempting to correct with DV = %.2f m/s\n", dvMax*1000);
	// 		itData = corrector.multShoot(&corNodes);
	// 		printf("Converged solution with DV = %.2f m/s\n", getTotalDV(&itData)*charV*1000);
	// 		corNodes = corrector.getBCR4BPR_Output();
			
	// 		// waitForUser();
	// 	}catch(tpat_diverge &e){
	// 		printf("Could not converge... exiting\n");
	// 		break;
	// 	}// End of DV loop try-catch
	// }// End of DV Loop

	// sprintf(filename, "data/LPF_4B_NaturalManifolds/Traj%03d_corrected_minDV.mat", trajIx);
	// corNodes.saveToMat(filename);
}
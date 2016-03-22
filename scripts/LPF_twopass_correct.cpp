/**
 *	Generate a bunch of manifold arcs and correct them to be natural in the BC4BP
 */

#include "tpat_all_includes.hpp"

#include <cmath>
#include <cstdio>

int main(void){

	tpat_sys_data_cr3bp seSys("sun", "earth");
	tpat_sys_data_cr3bp emSys("earth", "moon");

	// Load trajectory from file
	int trajIx = 39;
	char filename[128];
	sprintf(filename, "data/LPF_QH_4B_NaturalManifolds/Traj%03d_SEM.mat", trajIx);

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
	
	double earthIx = 1;
	tpat_event earthApse(&bcSys, tpat_event::APSE, 0, true, &earthIx);

	// propagate an arc to the SP YZ plane
	sim.addEvent(sp_yzPlane);
	sim.runSim(bcTraj.getState(0), bcTraj.getTime(0), tof);
	tpat_traj_bcr4bp arcToSP1 = sim.getBCR4BPR_Traj();

	// sim.clearEvents();
	// sim.addEvent(earthApse);
	// sim.runSim(arcToSP1.getState(-1), arcToSP1.getTime(-1), tof);
	// tpat_traj_bcr4bp arcToPerigee = sim.getBCR4BPR_Traj();

	// Propagate past Earth flyby to the SP YZ plane (again) or XZ plane
	sim.clearEvents();
	sim.addEvent(sp_yzPlane);
	sim.runSim(arcToSP1.getState(-1), arcToSP1.getTime(-1), tof);
	// sim.runSim(arcToPerigee.getState(-1), arcToPerigee.getTime(-1), tof);
	tpat_traj_bcr4bp arcToSP2 = sim.getBCR4BPR_Traj();

	// If y and y_dot have opposite signs, propagate to xz-plane
	if((arcToSP2.getState(-1)[1])*(arcToSP2.getState(-1)[4]) < 0){
		// Propagate a little further to XZ plane
		sim.clearEvents();
		sim.addEvent(sp_xzPlane);
		sim.runSim(arcToSP2.getState(-1), arcToSP2.getTime(-1), tof);
		tpat_traj_bcr4bp temp = sim.getBCR4BPR_Traj();
		arcToSP2 += temp;	// Add this extra arc to the old one
	}

	// arcToSP1.saveToMat("arcToSP1.mat");
	// arcToSP2.saveToMat("arcToSP2.mat");

	// Create nodesets from arcs
	int numNodes1 = 10, numNodes2 = 10;
	tpat_nodeset_bcr4bp nodesToSP1(arcToSP1.getState(0), &bcSys, arcToSP1.getTime(0), arcToSP1.getTOF(), numNodes1);
	// tpat_nodeset_bcr4bp nodesToPerigee(arcToPerigee.getState(0), &bcSys, arcToPerigee.getTime(0), arcToPerigee.getTOF(), 2);
	tpat_nodeset_bcr4bp nodesToSP2(arcToSP2.getState(0), &bcSys, arcToSP2.getTime(0), arcToSP2.getTOF(), numNodes2);

	// Concatenate nodes into one big set
	tpat_nodeset_bcr4bp allNodes = nodesToSP1;
	allNodes.deleteNode(-1);
	// allNodes += nodesToPerigee;
	// allNodes.deleteNode(-1);
	allNodes += nodesToSP2;

	sprintf(filename, "data/LPF_QH_4B_NaturalManifolds/Traj%03d_raw.mat", trajIx);
	allNodes.saveToMat(filename);

	// Create SP constraints and add to nodeset
	tpat_constraint sp1Con(tpat_constraint::SP, numNodes1-1, NULL, 0);
	tpat_constraint sp2Con(tpat_constraint::SP, allNodes.getNumNodes()-1, NULL, 0);
	
	std::vector<double> firstZ {NAN, NAN, allNodes.getState(0)[2], NAN, NAN, NAN};
	tpat_constraint fixFirst(tpat_constraint::STATE, 0, firstZ);

	allNodes.addConstraint(sp1Con);
	allNodes.addConstraint(sp2Con);

	// Allow velocity discontinuities
	std::vector<int> dvNodes {numNodes1 - 2, allNodes.getNumNodes()-5};
	allNodes.allowDV_at(dvNodes);

	// Correct to encounter SP
	tpat_correction_engine corrector;
	corrector.setTol(1e-11);
	iterationData itData = corrector.multShoot(&allNodes);
	double totalDV = getTotalDV(&itData);
	double charV = bcSys.getCharL()/bcSys.getCharT();
	printf("Corrected trajectory with Delta-V = %.4f m/s\n", 1000*totalDV*charV);

	tpat_nodeset_bcr4bp corNodes = corrector.getBCR4BPR_Output();
	sprintf(filename, "data/LPF_QH_4B_NaturalManifolds/Traj%03d_corrected.mat", trajIx);
	corNodes.saveToMat(filename);

	// Attempt to step DV down by iteratively correcting with a maxDV constraint
	printColor(BOLDBLACK, "\nBeginning Delta-V Step-Down Loop\n");

	// dvNodes.push_back(numNodes1);
	// allNodes.allowDV_at(dvNodes);

	double approxDV_dim = std::floor(totalDV*1000*charV)/1000;	// round to nearest 1 m/s
	tpat_constraint maxDVCon(tpat_constraint::MAX_DELTA_V, 0, &approxDV_dim, 1);	// wrong dimensions on DV, updated in loop

	// Step size is one order of magnitude less than the total DV
	double dvStep = pow(10, std::floor(log10(approxDV_dim)) - 1);
	double dvStepMin = 0.001;

	for(double dvMax = approxDV_dim; dvMax > 0.001; dvMax -= dvStep){
		double dvMax_nonDim = dvMax/charV;

		maxDVCon.setData(&dvMax_nonDim, 1);

		firstZ[2] = corNodes.getState(0)[2];
		fixFirst.setData(firstZ);
		corNodes.clearConstraints();
		corNodes.addConstraint(sp1Con);
		corNodes.addConstraint(sp2Con);
		corNodes.addConstraint(maxDVCon);
		corNodes.addConstraint(fixFirst);

		try{
			corrector.setMaxIts(40);
			// corrector.setVerbose(NO_MSG);
			// convergedNodes.print();
			printf("Attempting to correct with DV = %.2f m/s\n", dvMax*1000);
			itData = corrector.multShoot(&corNodes);
			printf("Converged solution with DV = %.2f m/s\n", getTotalDV(&itData)*charV*1000);
			corNodes = corrector.getBCR4BPR_Output();
			
			// waitForUser();
		}catch(tpat_diverge &e){

			if(dvStep/2 > dvStepMin){
				dvMax += dvStep;	// Reset
				dvStep /= 2;		// Reduce step size
			}else{
				printf("Could not converge... exiting\n");
				break;
			}
		}// End of DV loop try-catch

		while(dvMax - dvStep < 0.001){
			printf("Reducing step size:\n");
			printf("  dvMax is %.3f m/s\n", dvMax*1000);
			printf("  dvStep is %.4f m/s\n", dvStep*1000);
			
			if(dvStep/2 > dvStepMin){
				dvStep /= 2;
			}
			else{
				dvStep = dvStepMin;
				break;
			}
		}
	}// End of DV Loop

	sprintf(filename, "data/LPF_QH_4B_NaturalManifolds/Traj%03d_corrected_minDV.mat", trajIx);
	corNodes.saveToMat(filename);

	tpat_nodeset_bcr4bp corNodes_SE = bcr4bpr_SEM2SE(corNodes, &seSys);
	tpat_nodeset_cr3bp corNodes_EM = cr3bp_SE2EM(corNodes_SE, &emSys, bcSys.getEpoch0()*bcSys.getCharT()/seSys.getCharT(), 
		bcSys.getTheta0(), bcSys.getPhi0(), bcSys.getGamma());
	
	tpat_nodeset_cr3bp corNodes_ECI = cr3bp_rot2inert(corNodes_EM, 0);
	sprintf(filename, "data/LPF_QH_4B_NaturalManifolds/Traj%03d_corrected_minDV_ECI.mat", trajIx);
	corNodes_ECI.saveToMat(filename);
}
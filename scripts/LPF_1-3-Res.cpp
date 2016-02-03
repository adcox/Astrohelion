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
	tpat_sys_data_cr3bp seSys("sun", "earth");
	tpat_sys_data_cr3bp emSys("earth", "moon");
	tpat_sys_data_bcr4bpr bcSys("Sun", "Earth", "Moon");
	double leaveHaloEpoch = dateToEpochTime("2016/04/18");

	// double t0 = leaveHaloEpoch/bcSys.getCharT();	 // (Big number, let's avoid this)
	double t0 = 0;

	int numHaloNodes = 8;
	int numManNodes = 10;
	int numResNodes = 40;

	// Simulate from SE line backwards to beginning of manifold
	double IC[] = {1.001990975157441, 0, 0.000266392751420, 0.019779727575305, 0.041544455161967, 0.007099432234216};
	double tf = -2.618819626347190;

	double resIC[] = {0.661595270304142, 0.447542832021962, 0, 0.280878992941660, 0.759320138643845, 0};
	double resEpoch = 5.274288671852568e+08;
	double resTOF = 18.521059309999998;

	// Orient BC4BP system for the correct epoch
	bcr4bpr_orientAtEpoch(leaveHaloEpoch, &bcSys);

	// Create the resonant orbit via simulation
	tpat_simulation_engine engine(&emSys);
	// engine.runSim(resIC, 2.1*resTOF);
	// tpat_traj_cr3bp resTraj = engine.getCR3BP_Traj();
	// resTraj.shiftAllTimes((resEpoch - leaveHaloEpoch)/emSys.getCharT());
	// tpat_traj_cr3bp resTraj_SE = cr3bp_EM2SE(resTraj, &seSys, bcSys.getTheta0(), bcSys.getPhi0(), bcSys.getGamma());
	// resTraj.saveToMat("resTraj_EM.mat");
	// resTraj_SE.saveToMat("resTraj_SE.mat");

	tpat_nodeset_cr3bp resNodes(resIC, &emSys, 2*resTOF, numResNodes);	// Propagate in EM system
	tpat_nodeset_cr3bp resNodes_SE = cr3bp_EM2SE(resNodes, &seSys, (resEpoch-leaveHaloEpoch)/emSys.getCharT(),
		bcSys.getTheta0(), bcSys.getPhi0(), bcSys.getGamma());
	tpat_nodeset_bcr4bp resNodes_4B = bcr4bpr_SE2SEM(resNodes_SE, &bcSys, (resEpoch - leaveHaloEpoch)/bcSys.getCharT());

	// resNodes.saveToMat("resNodes_EM.mat");
	// resNodes_SE.saveToMat("resNodes_SE.mat");

	// waitForUser();

	// Generate nodesets for the halo manifold and the halo itself, both in reverse time
	tpat_nodeset_cr3bp manData(IC, &seSys, tf, numManNodes);
	tpat_nodeset_cr3bp haloData(manData.getNode(-1).getPosVelState(), &seSys, -PI, numHaloNodes);

	// Concatenate nodesets;
	manData.deleteNode(-1);
	tpat_nodeset_cr3bp manToHaloData = manData + haloData;

	// Change to BCR4BPR System
	tpat_nodeset_bcr4bp manToHalo_bcNodes = bcr4bpr_SE2SEM(manToHaloData, &bcSys, t0);
	manToHalo_bcNodes.reverseOrder();

	// Make all nodes continuous in velocity
	std::vector<int> notContinuous = std::vector<int>();
	manToHalo_bcNodes.allowDV_at(notContinuous);

	// Correct to be continuous
	tpat_correction_engine corrector = tpat_correction_engine();
	corrector.setVerbose(SOME_MSG);
	corrector.setTol(1e-10);
	corrector.setScaleVars(false);
	corrector.multShoot(&manToHalo_bcNodes);
	tpat_nodeset_bcr4bp contNodes = corrector.getBCR4BPR_Output();
	
	corrector.multShoot(&resNodes_4B);
	tpat_nodeset_bcr4bp resNodes_4BC = corrector.getBCR4BPR_Output();

	// double maxTOF = 5.0/bcSys.getK();
	// engine.addEvent(tpat_event::XZ_PLANE, -1, true);
	// engine.runSim(contNodes.getNode(-1).getPosVelState(), t0, maxTOF);
	// tpat_traj_bcr4bp newArc = engine.getBCR4BPR_Traj();
	
	// Concatenate resonant orbit nodes 
	contNodes.deleteNode(-1);
	tpat_nodeset_bcr4bp fullNodeset = contNodes + resNodes_4BC;

	// // Create a nodeset from the arc segment
	// double newArcTOF = newArc.getTime(-1);
	// double numNewArcNodes = 5;
	// tpat_nodeset_bcr4bp evolveSegNodes(contNodes.getNode(-1).getPosVelState(), &bcSys, t0, newArcTOF, numNewArcNodes);

	// // Concatenate into one large nodeset
	// contNodes.deleteNode(-1);
	// tpat_nodeset_bcr4bp fullNodeset = contNodes + evolveSegNodes;

	fullNodeset.saveToMat("temp.mat");
	waitForUser();

	// Get a subset of the full nodeset for corrections; don't want the halo nodes to change
	tpat_nodeset_bcr4bp workingNodes(fullNodeset, numHaloNodes-1, fullNodeset.getNumNodes());

	// Constraint to fix first node and intersect SP twice
	double maxR = 40;	// km
	double maxA = (5.358e-8*maxR - 2.3313e-8)/1000*bcSys.getCharT()*bcSys.getCharT()/bcSys.getCharL();
	std::vector<double> firstState = workingNodes.getNode(0).getPosVelState();
	firstState.push_back(workingNodes.getEpoch(0));
	tpat_constraint fixFirst(tpat_constraint::STATE, 0, firstState);
	int sp1_node = numManNodes-1;
	int sp2_node = workingNodes.getNumNodes()-1;
	// tpat_constraint sp1(tpat_constraint::SP_RANGE, sp1_node, &maxA, 1);
	// tpat_constraint sp2(tpat_constraint::SP_RANGE, sp2_node, &maxA, 1);
	// tpat_constraint sp1(tpat_constraint::SP, sp1_node, workingNodes.getNode(sp1_node).getPosVelState());
	// tpat_constraint sp2(tpat_constraint::SP, sp2_node, workingNodes.getNode(sp2_node).getPosVelState());
	double maxDist = 100/bcSys.getCharL();		// Set to 100 km to start with
	MatrixXRd coeff1(3,3);						// Empty matrix, for now
	MatrixXRd coeff2(3,3);						// Empty matrix, for now
	coeff1 = bcr4bpr_spLoc_polyFit(&bcSys, workingNodes.getEpoch(sp1_node));
	coeff2 = bcr4bpr_spLoc_polyFit(&bcSys, workingNodes.getEpoch(sp2_node));
	double d1[] = {maxDist, coeff1(0,0), coeff1(1,0), coeff1(2,0), coeff1(0,1), coeff1(1,1),
		coeff1(2,1), coeff1(0,2), coeff1(1,2), coeff1(2,2)};
	double d2[] = {maxDist, coeff2(0,0), coeff2(1,0), coeff2(2,0), coeff2(0,1), coeff2(1,1),
		coeff2(2,1), coeff2(0,2), coeff2(1,2), coeff2(2,2)};
	tpat_constraint sp1(tpat_constraint::SP_MAX_DIST, sp1_node, d1, 10);
	tpat_constraint sp2(tpat_constraint::SP_MAX_DIST, sp2_node, d2, 10);


	// Correct to satisfy new constraints
	workingNodes.addConstraint(fixFirst);
	workingNodes.addConstraint(sp1);
	workingNodes.addConstraint(sp2);
	notContinuous.push_back(sp1.getNode()-4);
	notContinuous.push_back(workingNodes.getNumNodes()-3);
	workingNodes.allowDV_at(notContinuous);
	workingNodes.print();
	workingNodes.saveToMat("WorkingNodes.mat");

	printColor(BOLDBLACK, "\nAdding SP constraints and correcting...\n");
	// corrector.setVerbose(ALL_MSG);
	corrector.setMaxIts(100);
	iterationData itBase = corrector.multShoot(&workingNodes);
	tpat_nodeset_bcr4bp convergedNodes = corrector.getBCR4BPR_Output();

	double totalDV = getTotalDV(&itBase);
	double charV = bcSys.getCharL()/bcSys.getCharT();
	printf("Converged solution has a total delta-V of %.4f m/s\n", 1000*totalDV*bcSys.getCharL()/bcSys.getCharT());

	convergedNodes.saveToMat("MinDV_LPF_Sample.mat");
	waitForUser();

	// Step Total Delta-V Down
	
	double approxDV_dim = std::floor(totalDV*100*charV)/100;	// round to nearest 10 m/s
	printf("Rounded dv = %.2f m/s\n", approxDV_dim*1000);

	// Create dummy distance inequality constraints for Saddle Point targeting
	double sp1ConData[] = {0,0,0,0,0,0,0,0,0,0};
	tpat_constraint sp1Ineq(tpat_constraint::SP_MAX_DIST, sp1_node, sp1ConData, 10);
	tpat_constraint sp2Ineq(tpat_constraint::SP_MAX_DIST, sp2_node, sp1ConData, 10);

	// Create dummy constraint for maxDV; need to put in maximum DV in NONDIM units
	tpat_constraint maxDVCon(tpat_constraint::MAX_DELTA_V, 0, &approxDV_dim, 1);
	
	// Initially, set active SP constraint to be the exact targeter
	tpat_constraint activeSP1_Con = sp1;
	tpat_constraint activeSP2_Con = sp2;

	maxDist = 0;							// Maximum distance from Saddle Point

	printColor(BOLDBLACK, "\nBeginning Delta-V Step-Down Loop\n");

	for(double dvMax = approxDV_dim; dvMax > 0.01; dvMax -= 0.001){
		double dvMax_nonDim = dvMax/charV;
		
		// Update coefficients for new epochs every iteration
		if(activeSP1_Con.getType() == tpat_constraint::SP_MAX_DIST){
			printf("Top Loop: Updating SP_MAX_DIST coefficients\n");
			coeff1 = bcr4bpr_spLoc_polyFit(&bcSys, convergedNodes.getEpoch(sp1_node));
			coeff2 = bcr4bpr_spLoc_polyFit(&bcSys, convergedNodes.getEpoch(sp2_node));
			double d1[] = {maxDist, coeff1(0,0), coeff1(1,0), coeff1(2,0), coeff1(0,1), coeff1(1,1),
				coeff1(2,1), coeff1(0,2), coeff1(1,2), coeff1(2,2)};
			double d2[] = {maxDist, coeff2(0,0), coeff2(1,0), coeff2(2,0), coeff2(0,1), coeff2(1,1),
				coeff2(2,1), coeff2(0,2), coeff2(1,2), coeff2(2,2)};
			activeSP1_Con.setData(d1, 10);
			activeSP2_Con.setData(d2, 10);
		}

		maxDVCon.setData(&dvMax_nonDim, 1);
		convergedNodes.clearConstraints();
		convergedNodes.addConstraint(activeSP1_Con);
		convergedNodes.addConstraint(activeSP2_Con);
		convergedNodes.addConstraint(maxDVCon);
		convergedNodes.addConstraint(fixFirst);
		// finiteDiff_checkMultShoot(&convergedNodes);

		try{
			// corrector.setMaxIts(3);
			// corrector.setVerbose(NO_MSG);
			// convergedNodes.print();
			printf("Attempting to correct with DV = %.2f m/s and SP Dist = %.2f km\n",
						dvMax*1000, maxDist*bcSys.getCharL());
			iterationData itData = corrector.multShoot(&convergedNodes);
			printf("Converged solution with DV = %.2f m/s and SP Dist = %.2f km\n",
						getTotalDV(&itData)*charV*1000, maxDist*bcSys.getCharL());
			convergedNodes = corrector.getBCR4BPR_Output();
			
			// waitForUser();
		}catch(tpat_diverge &e){
			
			printf("Could not converge... Attempting new SP Con OR increased distance\n");
			bool converged = false;
			while(!converged){

				// If corrections still utilize exact targeting, switch to inequality and update coefficients
				if(activeSP1_Con.getType() == tpat_constraint::SP){
					// printf("Switching SP constraint type to SP_MAX_DIST\n");
					maxDist = 10/bcSys.getCharL();	// Set to 10 km to start with
					coeff1 = bcr4bpr_spLoc_polyFit(&bcSys, convergedNodes.getEpoch(sp1_node));
					coeff2 = bcr4bpr_spLoc_polyFit(&bcSys, convergedNodes.getEpoch(sp2_node));
					double d1[] = {maxDist, coeff1(0,0), coeff1(1,0), coeff1(2,0), coeff1(0,1), coeff1(1,1),
						coeff1(2,1), coeff1(0,2), coeff1(1,2), coeff1(2,2)};
					double d2[] = {maxDist, coeff2(0,0), coeff2(1,0), coeff2(2,0), coeff2(0,1), coeff2(1,1),
						coeff2(2,1), coeff2(0,2), coeff2(1,2), coeff2(2,2)};
					activeSP1_Con = sp1Ineq;
					activeSP2_Con = sp2Ineq;
					activeSP1_Con.setData(d1, 10);
					activeSP2_Con.setData(d2, 10);
				}else{
					// If already using inequality, update maxDist
					printf("Inner Loop: Updating SP_MAX_DIST coefficients\n");
					coeff1 = bcr4bpr_spLoc_polyFit(&bcSys, convergedNodes.getEpoch(sp1_node));
					coeff2 = bcr4bpr_spLoc_polyFit(&bcSys, convergedNodes.getEpoch(sp2_node));
					double d1[] = {maxDist, coeff1(0,0), coeff1(1,0), coeff1(2,0), coeff1(0,1), coeff1(1,1),
						coeff1(2,1), coeff1(0,2), coeff1(1,2), coeff1(2,2)};
					double d2[] = {maxDist, coeff2(0,0), coeff2(1,0), coeff2(2,0), coeff2(0,1), coeff2(1,1),
						coeff2(2,1), coeff2(0,2), coeff2(1,2), coeff2(2,2)};
					activeSP1_Con.setData(d1, 10);
					activeSP2_Con.setData(d2, 10);
				}
				
				convergedNodes.clearConstraints();
				convergedNodes.addConstraint(activeSP1_Con);
				convergedNodes.addConstraint(activeSP2_Con);
				convergedNodes.addConstraint(maxDVCon);
				convergedNodes.addConstraint(fixFirst);
				
				try{
					// convergedNodes.print();
					// finiteDiff_checkMultShoot(&convergedNodes);
					// waitForUser();
					printf("Attempting to correct with DV = %.2f m/s and SP Dist = %.2f km\n",
						dvMax*1000, maxDist*bcSys.getCharL());
					iterationData itData = corrector.multShoot(&convergedNodes);
					printf("Converged solution with DV = %.2f m/s and SP Dist = %.2f km\n",
						getTotalDV(&itData)*charV*1000, maxDist*bcSys.getCharL());
					convergedNodes = corrector.getBCR4BPR_Output();
					converged = true;

					waitForUser();
				}catch(tpat_diverge &e){
					if(maxDist < 100/bcSys.getCharL()){
						maxDist += 5/bcSys.getCharL();
						printf("  Could not converge, increasing SP Dist to %.2f km\n", maxDist*bcSys.getCharL());
					}
					else{
						printErr("Could not increase SP pass distance any more... exiting!");
						break;
					}
				}
			}// End of convergence loop

			// If convergence loop was exited instead of satisfied, kill the entire corrections process
			if(!converged)
				break;
		}// End of DV loop try-catch
	}// End of DV Loop

	convergedNodes.saveToMat("MinDV_LPF_Sample.mat");
}
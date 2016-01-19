/**
 *	Attempt to replicate LPF trajectory
 *
 *	To compile: g++ --std=c++11 -ltpat -Wall -pedantic SampleLPFTraj.cpp -o a.out
 */

#include "tpat_all_includes.hpp"

#include <cmath>
#include <cstdio>

int main(void){
	tpat_sys_data_cr3bp crSys("sun", "earth");
	tpat_sys_data_bcr4bpr bcSys("Sun", "Earth", "Moon");
	double leaveHaloEpoch = dateToEpochTime("2016/04/18");
	

	// Determine geometry at t0 so that we can start at T = 0;
	bcr4bpr_orientAtEpoch(leaveHaloEpoch, &bcSys);
	// double t0 = leaveHaloEpoch/bcSys.getCharT();	 // (Big number, let's avoid this)
	double t0 = 0;

	int numHaloNodes = 8;
	int numManNodes = 10;

	// Simulate from SE line backwards to beginning of manifold
	double IC[] = {1.000287361135506, 0, 0.000009622676738, 0.043512678073230,
		0.132592205936539, 0.019177430186551};
	double tf = -2.566028939144172;

	// Generate nodesets for the halo manifold and the halo itself, both in reverse time
	tpat_nodeset_cr3bp manData(IC, &crSys, tf, numManNodes);
	tpat_nodeset_cr3bp haloData(manData.getNode(-1).getPosVelState(), &crSys, -PI, numHaloNodes);

	// manData.print();
	// haloData.print();

	// Concatenate nodesets;
	manData.deleteNode(-1);
	tpat_nodeset_cr3bp manToHaloData = manData + haloData;

	// manToHaloData.print();

	// Change to BCR4BPR System
	tpat_nodeset_bcr4bp manToHalo_bcNodes = bcr4bpr_SE2SEM(manToHaloData, &bcSys, t0);
	// manToHalo_bcNodes.print();
	manToHalo_bcNodes.reverseOrder();
	// manToHalo_bcNodes.print();

	// Get a node that (roughly) splits the last segment in half
	tpat_simulation_engine engine(&bcSys);
	double segTOF = 14.4*manToHalo_bcNodes.getTOF(-2)/16.0;
	engine.runSim(manToHalo_bcNodes.getNode(-2).getPosVelState(), manToHalo_bcNodes.getEpoch(-2), segTOF);
	tpat_traj_bcr4bp finalSeg = engine.getBCR4BPR_Traj();

	// Make the segment leading up to this new node shorter (TOF)
	tpat_node nodeBefore = manToHalo_bcNodes.getNode(-2);

	// Create the intermediate node
	tpat_node node(finalSeg.getState(-1), nodeBefore.getTOF() - segTOF);
	nodeBefore.setTOF(segTOF);

	// Delete second to last node, replace with edited nodeBefore
	manToHalo_bcNodes.deleteNode(-2);
	manToHalo_bcNodes.insertNode(-1, nodeBefore);

	// Insert intermediate node before the final node
	manToHalo_bcNodes.insertNode(-1, node);

	// Make all nodes continuous in velocity
	std::vector<int> notContinuous = std::vector<int>();
	manToHalo_bcNodes.setVelConNodes_allBut(notContinuous);

	// manToHalo_bcNodes.print();

	// Correct to be continuous
	tpat_correction_engine corrector = tpat_correction_engine();
	corrector.setVerbose(SOME_MSG);
	corrector.setTol(1e-10);
	corrector.setScaleVars(false);
	corrector.multShoot(&manToHalo_bcNodes);
	tpat_nodeset_bcr4bp contNodes = corrector.getBCR4BPR_Output();

	// Propagate past final node until a negative XZ plane crossing (hopefully near SP)
	double maxTOF = 5.0/bcSys.getK();
	engine.addEvent(tpat_event::XZ_PLANE, -1, true);
	engine.runSim(contNodes.getNode(-1).getPosVelState(), t0, maxTOF);
	tpat_traj_bcr4bp newArc = engine.getBCR4BPR_Traj();
	
	// Create a nodeset from the arc segment
	double newArcTOF = newArc.getTime(-1);
	double numNewArcNodes = 5;
	tpat_nodeset_bcr4bp evolveSegNodes(contNodes.getNode(-1).getPosVelState(), &bcSys, t0, newArcTOF, numNewArcNodes);

	// Concatenate into one large nodeset
	contNodes.deleteNode(-1);
	tpat_nodeset_bcr4bp fullNodeset = contNodes + evolveSegNodes;

	// Get a subset of the full nodeset for corrections; don't want the halo nodes to change
	tpat_nodeset_bcr4bp workingNodes(fullNodeset, numHaloNodes-1, fullNodeset.getNumNodes());

	// Constraint to fix first node and intersect SP twice
	// double maxR = 40;	// km
	// double maxA = (5.358e-8*maxR - 2.3313e-8)/1000*bcSys.getCharT()*bcSys.getCharT()/bcSys.getCharL();
	std::vector<double> firstState = workingNodes.getNode(0).getPosVelState();
	firstState.push_back(workingNodes.getEpoch(0));
	tpat_constraint fixFirst(tpat_constraint::STATE, 0, firstState);
	int sp1_node = numManNodes-1;
	int sp2_node = workingNodes.getNumNodes()-1;
	tpat_constraint sp1(tpat_constraint::SP, sp1_node, workingNodes.getNode(sp1_node).getPosVelState());
	tpat_constraint sp2(tpat_constraint::SP, sp2_node, workingNodes.getNode(sp2_node).getPosVelState());

	// Correct to satisfy new constraints
	workingNodes.addConstraint(fixFirst);
	workingNodes.addConstraint(sp1);
	workingNodes.addConstraint(sp2);
	notContinuous.push_back(sp1.getNode()-4);
	notContinuous.push_back(workingNodes.getNumNodes()-3);
	workingNodes.setVelConNodes_allBut(notContinuous);
	workingNodes.print();
	workingNodes.saveToMat("WorkingNodes.mat");

	printColor(BOLDBLACK, "\nAdding SP constraints and correcting...\n");
	// corrector.setVerbose(ALL_MSG);
	corrector.setMaxIts(100);
	iterationData itBase = corrector.multShoot(&workingNodes);
	
	double totalDV = getTotalDV(&itBase);
	double charV = bcSys.getCharL()/bcSys.getCharT();
	printf("Converged solution has a total delta-V of %.4f m/s\n", 1000*totalDV*bcSys.getCharL()/bcSys.getCharT());

	// Step Total Delta-V Down
	tpat_nodeset_bcr4bp convergedNodes = corrector.getBCR4BPR_Output();
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

	double maxDist = 0;							// Maximum distance from Saddle Point
	MatrixXRd coeff1(3,3);						// Empty matrix, for now
	MatrixXRd coeff2(3,3);						// Empty matrix, for now

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
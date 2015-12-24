#include "tpat_ascii_output.hpp"
#include "tpat_calculations.hpp"
#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_nodeset_bcr4bp.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_bcr4bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_utilities.hpp"

#include <vector>
#include <iostream>

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

/**
 *  @brief Determine if the difference bewteen two state vectors is less than 
 *  the desired tolerance
 * 
 *  @param data vector of state values from the corrections process
 *  @param correct array of state values from the constraint
 *  @param tol desired numerical tolerance
 *  @return whether or not <tt>data</tt> and <tt>correct</tt> are equal
 *  within the desired tolerance
 */
bool stateDiffBelowTol(std::vector<double> data, double *correct, double tol){
	double sum = 0;
	for(size_t i = 0; i < data.size(); i++){
		if(!isnan(correct[i]))
			sum += pow(data[i] - correct[i], 2);
	}
	return sqrt(sum) < tol;
}

bool stateDiffBelowTol(std::vector<double> data, std::vector<double> correct, double tol){
	return stateDiffBelowTol(data, &(correct[0]), tol);
}

void finiteDiff(tpat_nodeset *nodeset){
	printf("Finite Diff: Checking DF matrix... ");
	// Create multiple shooter that will only do 1 iteration
	tpat_correction_engine corrector;
	corrector.setMaxIts(1);
	corrector.setVerbose(NO_MSG);
	corrector.setIgnoreDiverge(true);

	// Run multiple shooter to get X, FX, and DF
	iterationData it;
	it = corrector.multShoot(nodeset);
	Eigen::VectorXd FX = Eigen::Map<Eigen::VectorXd>(&(it.FX[0]), it.totalCons, 1);
	MatrixXRd DF = Eigen::Map<MatrixXRd>(&(it.DF[0]), it.totalCons, it.totalFree);

	double pertSize = 1e-8;
	MatrixXRd DFest = MatrixXRd::Zero(it.totalCons, it.totalFree);
	for(int i = 0; i < it.totalFree; i++){
		std::vector<double> pertX = it.X0;
		pertX[i] += pertSize;
		it.X = pertX;
		iterationData pertIt = corrector.multShoot(it);
		Eigen::VectorXd newFX = Eigen::Map<Eigen::VectorXd>(&(pertIt.FX[0]), it.totalCons, 1);
		Eigen::VectorXd col = (newFX - FX)/std::abs(pertSize);
		DFest.block(0, i, it.totalCons, 1) = col;
	}

	MatrixXRd diff = DF - DFest;
	diff = diff.cwiseAbs();
	Eigen::VectorXd rowMax = diff.rowwise().maxCoeff();
	Eigen::RowVectorXd colMax = diff.colwise().maxCoeff();

	double rowMaxMax = rowMax.maxCoeff();
	double colMaxMax = colMax.maxCoeff();
	int errScalar = 3000;

	if(rowMaxMax < errScalar*pertSize && colMaxMax < errScalar*colMaxMax){
		printColor(BOLDGREEN, "No significant errors!\n");
	}else{
		printColor(BOLDRED, "Significant errors!\n");
		printf("Maximum Difference between computed DF and estimated DF\n");
		int conCount = 0;
		for(long r = 0; r < rowMax.size(); r++){
			if(r == 0 && it.totalCons > 0){
				printf("%s Constraint:\n", it.allCons[conCount].getTypeStr());
			}else if(conCount < it.totalCons && r >= it.conRows[conCount+1]){
				conCount++;
				printf("%s Constraint:\n", it.allCons[conCount].getTypeStr());
			}
			printColor(rowMax[r] > errScalar*pertSize ? RED : GREEN, "  row %03zu: %.6e\n", r, rowMax[r]);
		}
		for(long c = 0; c < colMax.size(); c++){
			printColor(colMax[c] > errScalar*pertSize ? RED : GREEN, "Free Var %03zu: %.6e\n", c, colMax[c]);
		}
	}
}//================================================

void testCR3BPCons(){
	tpat_sys_data_cr3bp sys("earth", "moon");
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period
	tpat_nodeset_cr3bp halfLyapNodeset(ic, &sys, T/1.25, 8);	// Create a nodeset with three nodes

	std::vector<double> initState, finalState;
	tpat_nodeset_cr3bp correctedSet(&sys);

	tpat_correction_engine corrector;
	printColor(BOLDBLACK, "Testing CR3BP Multiple Shooting Constraints\n");

	// STATE
	printColor(BOLDBLACK, "STATE Constraint\n");
	double stateConData[] = {0.9, 0.1, NAN, NAN, NAN, NAN};
	tpat_constraint stateCon(tpat_constraint::STATE, 7, stateConData, 6);
	halfLyapNodeset.addConstraint(stateCon);
	finiteDiff(&halfLyapNodeset);
	try{
		corrector.multShoot(&halfLyapNodeset);
		correctedSet = corrector.getCR3BP_Output();
		finalState = correctedSet.getState(stateCon.getNode());
		std::cout << "STATE Constraint: " << (stateDiffBelowTol(finalState, stateConData, 1e-12) ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "STATE Constraint: " << FAIL << std::endl;
	}

	// MATCH_ALL
	printColor(BOLDBLACK, "MATCH_ALL Constraint\n");
	double matchAllConData = 0;
	tpat_constraint matchAllCon(tpat_constraint::MATCH_ALL, 7, &matchAllConData, 1);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(matchAllCon);
	finiteDiff(&halfLyapNodeset);
	try{
		corrector.multShoot(&halfLyapNodeset);
		correctedSet = corrector.getCR3BP_Output();
		finalState = correctedSet.getState(matchAllCon.getNode());
		initState = correctedSet.getState(0);
		std::cout << "MATCH_ALL Constraint: " << (stateDiffBelowTol(finalState, initState, 1e-12) ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "MATCH_ALL Constraint: " << FAIL << std::endl;
	}

	// MATCH_CUST
	printColor(BOLDBLACK, "MATCH_CUST Constraint\n");
	double matchCustConData[] = {0,0,NAN,NAN,NAN,NAN};
	tpat_constraint matchCustCon(tpat_constraint::MATCH_CUST, 7, matchCustConData, 6);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(matchCustCon);
	finiteDiff(&halfLyapNodeset);
	try{
		corrector.multShoot(&halfLyapNodeset);
		correctedSet = corrector.getCR3BP_Output();
		finalState = correctedSet.getState(matchCustCon.getNode());
		initState = correctedSet.getState(0);
		finalState.erase(finalState.begin()+2, finalState.end());	// Erase entries 2 through 5; we're only comparing the first two
		initState.erase(initState.begin()+2, initState.end());
		std::cout << "MATCH_CUST Constraint: " << (stateDiffBelowTol(finalState, initState, 1e-12) ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "MATCH_CUST Constraint: " << FAIL << std::endl;
	}

	// DIST
	printColor(BOLDBLACK, "DIST Constraint\n");
	double matchDistConData[] = {1, 0.2};
	tpat_constraint matchDistCon(tpat_constraint::DIST, 5, matchDistConData, 2);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(matchDistCon);
	finiteDiff(&halfLyapNodeset);
	try{
		corrector.multShoot(&halfLyapNodeset);
		correctedSet = corrector.getCR3BP_Output();
		finalState = correctedSet.getState(matchDistCon.getNode());
		double dist = sqrt(pow(finalState[0] - 1 + sys.getMu(),2) + pow(finalState[1], 2));
		std::cout << "DIST Constraint: " << (std::abs(dist - matchDistConData[1]) < 1e-12 ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "DIST Constraint: " << FAIL << std::endl;
	}

	// MIN_DIST
	printColor(BOLDBLACK, "MIN_DIST Constraint\n");
	matchDistConData[1] = 0.1;
	matchDistCon.setData(matchDistConData, 2);
	matchDistCon.setType(tpat_constraint::MIN_DIST);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(matchDistCon);
	finiteDiff(&halfLyapNodeset);
	try{
		corrector.multShoot(&halfLyapNodeset);
		correctedSet = corrector.getCR3BP_Output();
		finalState = correctedSet.getState(matchDistCon.getNode());
		double dist = sqrt(pow(finalState[0] - 1 + sys.getMu(), 2) + pow(finalState[1], 2));
		std::cout << "MIN_DIST Constraint: " << (dist >= matchDistConData[1] ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "MIN_DIST Constraint: " << FAIL << std::endl;
	}

	// MAX_DIST
	printColor(BOLDBLACK, "MAX_DIST Constraint\n");
	matchDistConData[1] = 0.3;
	matchDistCon.setData(matchDistConData, 2);
	matchDistCon.setType(tpat_constraint::MAX_DIST);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(matchDistCon);
	finiteDiff(&halfLyapNodeset);
	try{
		corrector.multShoot(&halfLyapNodeset);
		correctedSet = corrector.getCR3BP_Output();
		finalState = correctedSet.getState(matchDistCon.getNode());
		double dist = sqrt(pow(finalState[0] - 1 + sys.getMu(), 2) + pow(finalState[1], 2));
		std::cout << "MAX_DIST Constraint: " << (dist <= matchDistConData[1] ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "MAX_DIST Constraint: " << FAIL << std::endl;
	}

	// MAX_DELTA_V
	printColor(BOLDBLACK, "MAX_DELTA_V Constraint\n");
	std::vector<double> state = halfLyapNodeset.getState(6);
	state[3] += 0.01;
	state[4] += 0.1;
	halfLyapNodeset.setState(6, state);	// Perturb the velocity of this state to create a discontinuity
	std::vector<int> dvNodes {6};
	halfLyapNodeset.setVelConNodes_allBut(dvNodes);	// Allow the perturbed node to have a delta-v
	double maxDVConData = 0.01*0.01;
	tpat_constraint dVCon(tpat_constraint::MAX_DELTA_V, 0, &maxDVConData, 1);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(dVCon);
	finiteDiff(&halfLyapNodeset);
	double totalDV = 0;
	try{
		iterationData itData = corrector.multShoot(&halfLyapNodeset);
		for(size_t i = 0; i < itData.deltaVs.size(); i+=3){
			totalDV += pow(itData.deltaVs[i], 2) + pow(itData.deltaVs[i+1], 2) +
				pow(itData.deltaVs[i+2], 2);
		}
		std::cout << "MAX_DELTA_V Constraint: " << (totalDV <= maxDVConData ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "MAX_DELTA_V Constraint: " << FAIL << std::endl;
	}

	// DELTA_V
	printColor(BOLDBLACK, "DELTA_V Constraint\n");
	maxDVConData = 0.02*0.02;
	dVCon.setData(&maxDVConData, 1);
	dVCon.setType(tpat_constraint::DELTA_V);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(dVCon);
	finiteDiff(&halfLyapNodeset);
	totalDV = 0;
	try{
		iterationData itData = corrector.multShoot(&halfLyapNodeset);
		for(size_t i = 0; i < itData.deltaVs.size(); i+=3){
			totalDV += pow(itData.deltaVs[i], 2) + pow(itData.deltaVs[i+1], 2) +
				pow(itData.deltaVs[i+2], 2);
		}
		std::cout << "DELTA_V Constraint: " << (std::abs(totalDV - maxDVConData) < 1e-12 ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "DELTA_V Constraint: " << FAIL << std::endl;
		printf("  Total DV = %.4e\n", totalDV);
	}
	
	// JC
	printColor(BOLDBLACK, "JC Constraint\n");
	double jacobiData = 3.1149;
	tpat_constraint jacobiCon(tpat_constraint::JC, 0, &jacobiData, 1);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(jacobiCon);
	finiteDiff(&halfLyapNodeset);
	try{
		corrector.multShoot(&halfLyapNodeset);
		correctedSet = corrector.getCR3BP_Output();
		double jacobi = correctedSet.getJacobi(0);
		std::cout << "JC Constraint: " << (std::abs(jacobi - jacobiData) < 1e-12 ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "JC Constraint: " << FAIL << std::endl;
	}

	// TOF
	printColor(BOLDBLACK, "TOF Constraint\n");
	double tofData = 2.5;
	tpat_constraint tofCon(tpat_constraint::TOF, 0, &tofData, 1);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(tofCon);
	finiteDiff(&halfLyapNodeset);
	try{
		corrector.multShoot(&halfLyapNodeset);
		correctedSet = corrector.getCR3BP_Output();
		double totalTOF = correctedSet.getTotalTOF();
		std::cout << "TOF Constraint: " << (std::abs(totalTOF - tofData) < 1e-12 ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "TOF Constraint: " << FAIL << std::endl;
	}

	// APSE
	printColor(BOLDBLACK, "APSE Constraint\n");
	double apseData = 1;
	tpat_constraint apseCon(tpat_constraint::APSE, 7, &apseData, 1);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(apseCon);
	finiteDiff(&halfLyapNodeset);
	try{
		corrector.multShoot(&halfLyapNodeset);
		correctedSet = corrector.getCR3BP_Output();
		finalState = correctedSet.getState(apseCon.getNode());
		tpat_model *model = sys.getModel();
		std::vector<double> primPos = model->getPrimPos(0, &sys);
		double dx = finalState[0] - primPos[apseData*3 + 0];
		double dy = finalState[1] - primPos[apseData*3 + 1];
		double dz = finalState[2] - primPos[apseData*3 + 2];
		double rdot = dx*finalState[3] + dy*finalState[4] + dz*finalState[5];
		std::cout << "APSE Constraint: " << (std::abs(rdot) < 1e-12 ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "APSE Constraint: " << FAIL << std::endl;
	}
}

void testBCR4BPCons(){
	tpat_sys_data_bcr4bpr sys("sun", "earth", "moon");
	double ic[] = {99.18886503808, 0.00267453509779252, 6.82940156780485e-05, -0.000123170516283809, -0.0112246703481597, -2.38334057821979e-06};	// SE L1
	double T = 205.1445;	// SE L1 Period
	tpat_nodeset_bcr4bp halfLyapNodeset(ic, &sys, 0, T, 8);	// Create a nodeset
	std::vector<double> initState, finalState;
	tpat_nodeset_bcr4bp correctedSet(&sys);

	tpat_correction_engine corrector;
	printColor(BOLDBLACK, "Testing BCR4BP Multiple Shooting Constraints\n");

	// // STATE
	// printColor(BOLDBLACK, "STATE Constraint\n");
	// double stateConData[] = {98.98, 0.45, NAN, NAN, NAN, NAN};
	// tpat_constraint stateCon(tpat_constraint::STATE, 7, stateConData, 6);
	// halfLyapNodeset.addConstraint(stateCon);
	// finiteDiff(&halfLyapNodeset);
	// try{
	// 	corrector.multShoot(&halfLyapNodeset);
	// 	correctedSet = corrector.getBCR4BPR_Output();
	// 	finalState = correctedSet.getState(stateCon.getNode());
	// 	std::cout << "STATE Constraint: " << (stateDiffBelowTol(finalState, stateConData, 1e-12) ? PASS : FAIL) << std::endl;
	// }catch(tpat_diverge &e){
	// 	std::cout << "STATE Constraint: " << FAIL << std::endl;
	// }

	// // MATCH_ALL
	// printColor(BOLDBLACK, "MATCH_ALL Constraint\n");
	// double matchAllConData = 0;
	// tpat_constraint matchAllCon(tpat_constraint::MATCH_ALL, 7, &matchAllConData, 1);
	// halfLyapNodeset.clearConstraints();
	// halfLyapNodeset.addConstraint(matchAllCon);
	// finiteDiff(&halfLyapNodeset);
	// try{
	// 	corrector.multShoot(&halfLyapNodeset);
	// 	correctedSet = corrector.getBCR4BPR_Output();
	// 	finalState = correctedSet.getState(matchAllCon.getNode());
	// 	initState = correctedSet.getState(0);
	// 	std::cout << "MATCH_ALL Constraint: " << (stateDiffBelowTol(finalState, initState, 1e-12) ? PASS : FAIL) << std::endl;
	// }catch(tpat_diverge &e){
	// 	std::cout << "MATCH_ALL Constraint: " << FAIL << std::endl;
	// }


	// // MATCH_CUST
	// printColor(BOLDBLACK, "MATCH_CUST Constraint\n");
	// double matchCustConData[] = {0,0,NAN,NAN,NAN,NAN};
	// tpat_constraint matchCustCon(tpat_constraint::MATCH_CUST, 7, matchCustConData, 6);
	// halfLyapNodeset.clearConstraints();
	// halfLyapNodeset.addConstraint(matchCustCon);
	// finiteDiff(&halfLyapNodeset);
	// try{
	// 	corrector.multShoot(&halfLyapNodeset);
	// 	correctedSet = corrector.getBCR4BPR_Output();
	// 	finalState = correctedSet.getState(matchCustCon.getNode());
	// 	initState = correctedSet.getState(0);
	// 	finalState.erase(finalState.begin()+2, finalState.end());	// Erase entries 2 through 5; we're only comparing the first two
	// 	initState.erase(initState.begin()+2, initState.end());
	// 	std::cout << "MATCH_CUST Constraint: " << (stateDiffBelowTol(finalState, initState, 1e-12) ? PASS : FAIL) << std::endl;
	// }catch(tpat_diverge &e){
	// 	std::cout << "MATCH_CUST Constraint: " << FAIL << std::endl;
	// }

	// // DIST
	// printColor(BOLDBLACK, "DIST Constraint\n");
	// double matchDistConData[] = {1, 1.0};
	// tpat_constraint matchDistCon(tpat_constraint::DIST, 5, matchDistConData, 2);
	// halfLyapNodeset.clearConstraints();
	// halfLyapNodeset.addConstraint(matchDistCon);
	// finiteDiff(&halfLyapNodeset);
	// try{
	// 	corrector.multShoot(&halfLyapNodeset);
	// 	correctedSet = corrector.getBCR4BPR_Output();
	// 	finalState = correctedSet.getState(matchDistCon.getNode());
	// 	std::vector<double> primPos = sys.getModel()->getPrimPos(correctedSet.getEpoch(matchDistCon.getNode()), &sys);
	// 	double dist = sqrt(pow(finalState[0] - primPos[3] ,2) + pow(finalState[1] - primPos[4], 2) + pow(finalState[2] - primPos[5], 2));
	// 	std::cout << "DIST Constraint: " << (std::abs(dist - matchDistConData[1]) < 1e-12 ? PASS : FAIL) << std::endl;
	// }catch(tpat_diverge &e){
	// 	std::cout << "DIST Constraint: " << FAIL << std::endl;
	// }

	// // MIN_DIST
	// printColor(BOLDBLACK, "MIN_DIST Constraint\n");
	// matchDistConData[1] = 1.1;
	// matchDistCon.setData(matchDistConData, 2);
	// matchDistCon.setType(tpat_constraint::MIN_DIST);
	// halfLyapNodeset.clearConstraints();
	// halfLyapNodeset.addConstraint(matchDistCon);
	// finiteDiff(&halfLyapNodeset);
	// try{
	// 	corrector.multShoot(&halfLyapNodeset);
	// 	correctedSet = corrector.getBCR4BPR_Output();
	// 	finalState = correctedSet.getState(matchDistCon.getNode());
	// 	std::vector<double> primPos = sys.getModel()->getPrimPos(correctedSet.getEpoch(matchDistCon.getNode()), &sys);
	// 	double dist = sqrt(pow(finalState[0] - primPos[3] ,2) + pow(finalState[1] - primPos[4], 2) + pow(finalState[2] - primPos[5], 2));
	// 	std::cout << "MIN_DIST Constraint: " << (dist >= matchDistConData[1] ? PASS : FAIL) << std::endl;
	// }catch(tpat_diverge &e){
	// 	std::cout << "MIN_DIST Constraint: " << FAIL << std::endl;
	// }

	// // waitForUser();

	// // MAX_DIST
	// printColor(BOLDBLACK, "MAX_DIST Constraint\n");
	// matchDistConData[1] = 0.9;
	// // matchDistConData[1] = 2;
	// matchDistCon.setData(matchDistConData, 2);
	// matchDistCon.setType(tpat_constraint::MAX_DIST);
	// halfLyapNodeset.clearConstraints();
	// halfLyapNodeset.addConstraint(matchDistCon);
	// finiteDiff(&halfLyapNodeset);
	// try{
	// 	corrector.multShoot(&halfLyapNodeset);
	// 	correctedSet = corrector.getBCR4BPR_Output();
	// 	finalState = correctedSet.getState(matchDistCon.getNode());
	// 	std::vector<double> primPos = sys.getModel()->getPrimPos(correctedSet.getEpoch(matchDistCon.getNode()), &sys);
	// 	double dist = sqrt(pow(finalState[0] - primPos[3] ,2) + pow(finalState[1] - primPos[4], 2) + pow(finalState[2] - primPos[5], 2));
	// 	std::cout << "MAX_DIST Constraint: " << (dist <= matchDistConData[1] ? PASS : FAIL) << std::endl;
	// }catch(tpat_diverge &e){
	// 	std::cout << "MAX_DIST Constraint: " << FAIL << std::endl;
	// }


	// // MAX_DELTA_V
	// printColor(BOLDBLACK, "MAX_DELTA_V Constraint\n");
	// std::vector<double> state = halfLyapNodeset.getState(6);
	// state[3] += 0.01;
	// state[4] += 0.1;
	// halfLyapNodeset.setState(6, state);	// Perturb the velocity of this state to create a discontinuity
	// std::vector<int> dvNodes {6};
	// halfLyapNodeset.setVelConNodes_allBut(dvNodes);	// Allow the perturbed node to have a delta-v
	// double maxDVConData = 0.03*0.03;
	// tpat_constraint dVCon(tpat_constraint::MAX_DELTA_V, 0, &maxDVConData, 1);
	// halfLyapNodeset.clearConstraints();
	// halfLyapNodeset.addConstraint(dVCon);
	// finiteDiff(&halfLyapNodeset);
	// double totalDV = 0;
	// try{
	// 	iterationData itData = corrector.multShoot(&halfLyapNodeset);
	// 	for(size_t i = 0; i < itData.deltaVs.size(); i+=3){
	// 		totalDV += pow(itData.deltaVs[i], 2) + pow(itData.deltaVs[i+1], 2) +
	// 			pow(itData.deltaVs[i+2], 2);
	// 	}
	// 	std::cout << "MAX_DELTA_V Constraint: " << (totalDV <= maxDVConData ? PASS : FAIL) << std::endl;
	// }catch(tpat_diverge &e){
	// 	std::cout << "MAX_DELTA_V Constraint: " << FAIL << std::endl;
	// }

	// // DELTA_V
	// printColor(BOLDBLACK, "DELTA_V Constraint\n");
	// maxDVConData = 0.02*0.02;
	// dVCon.setData(&maxDVConData, 1);
	// dVCon.setType(tpat_constraint::DELTA_V);
	// halfLyapNodeset.clearConstraints();
	// halfLyapNodeset.addConstraint(dVCon);
	// finiteDiff(&halfLyapNodeset);
	// totalDV = 0;
	// try{
	// 	iterationData itData = corrector.multShoot(&halfLyapNodeset);
	// 	correctedSet = corrector.getBCR4BPR_Output();
	// 	for(size_t i = 0; i < itData.deltaVs.size(); i+=3){
	// 		totalDV += pow(itData.deltaVs[i], 2) + pow(itData.deltaVs[i+1], 2) +
	// 			pow(itData.deltaVs[i+2], 2);
	// 	}
	// 	std::cout << "DELTA_V Constraint: " << (std::abs(totalDV - maxDVConData) < 1e-12 ? PASS : FAIL) << std::endl;
	// }catch(tpat_diverge &e){
	// 	std::cout << "DELTA_V Constraint: " << FAIL << std::endl;
	// 	printf("  Total DV = %.4e\n", totalDV);
	// }

	// // TOF
	// printColor(BOLDBLACK, "TOF Constraint\n");
	// double tofData = 2.5;
	// tpat_constraint tofCon(tpat_constraint::TOF, 0, &tofData, 1);
	// halfLyapNodeset.clearConstraints();
	// halfLyapNodeset.addConstraint(tofCon);
	// finiteDiff(&halfLyapNodeset);
	// try{
	// 	corrector.multShoot(&halfLyapNodeset);
	// 	correctedSet = corrector.getBCR4BPR_Output();
	// 	double totalTOF = correctedSet.getTotalTOF();
	// 	std::cout << "TOF Constraint: " << (std::abs(totalTOF - tofData) < 1e-12 ? PASS : FAIL) << std::endl;
	// }catch(tpat_diverge &e){
	// 	std::cout << "TOF Constraint: " << FAIL << std::endl;
	// }

	// // APSE
	// printColor(BOLDBLACK, "APSE Constraint\n");
	// double apseData = 1;
	// tpat_constraint apseCon(tpat_constraint::APSE, 7, &apseData, 1);
	// halfLyapNodeset.clearConstraints();
	// halfLyapNodeset.addConstraint(apseCon);
	// finiteDiff(&halfLyapNodeset);
	// try{
	// 	corrector.multShoot(&halfLyapNodeset);
	// 	correctedSet = corrector.getBCR4BPR_Output();
	// 	finalState = correctedSet.getState(apseCon.getNode());
	// 	tpat_model *model = sys.getModel();
	// 	std::vector<double> primPos = model->getPrimPos(0, &sys);
	// 	double dx = finalState[0] - primPos[apseData*3 + 0];
	// 	double dy = finalState[1] - primPos[apseData*3 + 1];
	// 	double dz = finalState[2] - primPos[apseData*3 + 2];
	// 	double rdot = dx*finalState[3] + dy*finalState[4] + dz*finalState[5];
	// 	std::cout << "APSE Constraint: " << (std::abs(rdot) < 1e-12 ? PASS : FAIL) << std::endl;
	// }catch(tpat_diverge &e){
	// 	std::cout << "APSE Constraint: " << FAIL << std::endl;
	// }

	// Saddle Point, Exact
	printColor(BOLDBLACK, "SP Constraint\n");
	double IC[] = {98.529563503083, 1.18065951206182, 0.00114414803935654, 0.270235338496738, -0.222413502966165, 1.71361901677108e-05};
	double t0 = 51.23235;
	double tof = 5;
	tpat_nodeset_bcr4bp nodes0(IC, &sys, t0, tof, 2);

	double spData = 0;
	tpat_constraint spCon(tpat_constraint::SP, 1, &spData, 1);
	nodes0.addConstraint(spCon);
	finiteDiff(&nodes0);
	try{
		corrector.multShoot(&nodes0);
		correctedSet = corrector.getBCR4BPR_Output();
		finalState = correctedSet.getState(spCon.getNode());
		finalState.erase(finalState.begin()+3, finalState.end());
		Eigen::Vector3d spPos = bcr4bpr_getSPLoc(&sys, correctedSet.getEpoch(spCon.getNode()));
		double diff = sqrt(pow(spPos(0) - finalState[0], 2) + pow(spPos(1) - finalState[1], 2) + pow(spPos(2) - finalState[2], 2));
		std::cout << "SP Constraint: " << (diff < 1e-10 ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "SP Constraint: " << FAIL << std::endl;
	}

	// Saddle Point, Range
	printColor(BOLDBLACK, "SP Range Constraint\n");
	double spConData = 100/sys.getCharL();
	tpat_constraint spConRange(tpat_constraint::SP_RANGE, 1, &spConData, 1);
	nodes0.clearConstraints();
	nodes0.addConstraint(spConRange);
	finiteDiff(&nodes0);
	try{
		corrector.multShoot(&nodes0);
		correctedSet = corrector.getBCR4BPR_Output();
		finalState = correctedSet.getState(spCon.getNode());
		finalState.erase(finalState.begin()+3, finalState.end());
		Eigen::Vector3d spPos = bcr4bpr_getSPLoc(&sys, correctedSet.getEpoch(spCon.getNode()));
		double diff = sqrt(pow(spPos(0) - finalState[0], 2) + pow(spPos(1) - finalState[1], 2) + pow(spPos(2) - finalState[2], 2));
		std::cout << "SP Range Constraint: " << (diff <= spConData ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "SP Range Constraint: " << FAIL << std::endl;
	}
}

/**
 *  @brief Test all constraint types available to ensure they converge correctly
 */
int main(void){
	// testCR3BPCons();
	testBCR4BPCons();

	return EXIT_SUCCESS;
}
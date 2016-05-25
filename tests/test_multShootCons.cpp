#include "tpat_ascii_output.hpp"
#include "tpat_calculations.hpp"
#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_multShoot_data.hpp"
#include "tpat_nodeset_bc4bp.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_bc4bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_bc4bp.hpp"
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
		if(!std::isnan(correct[i]))
			sum += pow(data[i] - correct[i], 2);
	}
	return sqrt(sum) < tol;
}

bool stateDiffBelowTol(std::vector<double> data, std::vector<double> correct, double tol){
	return stateDiffBelowTol(data, &(correct[0]), tol);
}

void testCR3BP_SE_Cons(bool equalArcTime){
	TPAT_Sys_Data_CR3BP sys("sun", "earth");
	double ic[] = {0.993986593871357, 0, 0, 0, -0.022325793891591, 0};	// SE L1
	double T = 3.293141367224790;
	TPAT_Nodeset_CR3BP halfLyapNodeset(ic, &sys, T/1.25, 8);	// Create a nodeset
	std::vector<double> initState, finalState;
	TPAT_Nodeset_CR3BP correctedSet(&sys);

	TPAT_Correction_Engine corrector;
	corrector.setEqualArcTime(equalArcTime);
	printColor(BOLDBLACK, "Testing Sun-Earth CR3BP Multiple Shooting Constraints\n");

	// TPAT_Constraint_Tp::STATE
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::STATE Constraint\n");
	double stateConData[] = {0.9934, 0.001, NAN, NAN, NAN, NAN};
	TPAT_Constraint stateCon(TPAT_Constraint_Tp::STATE, 7, stateConData, 6);
	halfLyapNodeset.addConstraint(stateCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_CR3BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		finalState = correctedSet.getState(stateCon.getID());
		std::cout << "TPAT_Constraint_Tp::STATE Constraint: " << (stateDiffBelowTol(finalState, stateConData, 1e-12) ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::STATE Constraint: " << FAIL << std::endl;
	}	
}//====================================================

void testCR3BP_EM_Cons(bool equalArcTime){

	TPAT_Sys_Data_CR3BP sys("earth", "moon");
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period
	TPAT_Nodeset_CR3BP halfLyapNodeset(ic, &sys, T, 6);	// Create a nodeset

	std::vector<double> initState, finalState;
	TPAT_Nodeset_CR3BP correctedSet(&sys);

	TPAT_Correction_Engine corrector;
	corrector.setEqualArcTime(equalArcTime);
	printColor(BOLDBLACK, "Testing Earth-Moon CR3BP Multiple Shooting Constraints\n");

	// TPAT_Constraint_Tp::STATE
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::STATE Constraint\n");
	double stateConData[] = {0.9, 0.1, NAN, NAN, NAN, NAN};
	TPAT_Constraint stateCon(TPAT_Constraint_Tp::STATE, 4, stateConData, 6);
	halfLyapNodeset.addConstraint(stateCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		// corrector.setVerbose(TPAT_Verbosity_Tp::ALL_MSG);
		TPAT_Nodeset_CR3BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		finalState = correctedSet.getState(stateCon.getID());
		std::cout << "TPAT_Constraint_Tp::STATE Constraint: " << (stateDiffBelowTol(finalState, stateConData, 1e-12) ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::STATE Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::MATCH_ALL
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::MATCH_ALL Constraint\n");
	double matchAllConData = 0;
	TPAT_Constraint matchAllCon(TPAT_Constraint_Tp::MATCH_ALL, 4, &matchAllConData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(matchAllCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_CR3BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		finalState = correctedSet.getState(matchAllCon.getID());
		initState = correctedSet.getStateByIx(0);
		std::cout << "TPAT_Constraint_Tp::MATCH_ALL Constraint: " << (stateDiffBelowTol(finalState, initState, 1e-12) ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::MATCH_ALL Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::MATCH_CUST
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::MATCH_CUST Constraint\n");
	double matchCustConData[] = {0,0,NAN,NAN,NAN,NAN};
	TPAT_Constraint matchCustCon(TPAT_Constraint_Tp::MATCH_CUST, 4, matchCustConData, 6);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(matchCustCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_CR3BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		finalState = correctedSet.getState(matchCustCon.getID());
		initState = correctedSet.getStateByIx(0);
		finalState.erase(finalState.begin()+2, finalState.end());	// Erase entries 2 through 5; we're only comparing the first two
		initState.erase(initState.begin()+2, initState.end());
		std::cout << "TPAT_Constraint_Tp::MATCH_CUST Constraint: " << (stateDiffBelowTol(finalState, initState, 1e-12) ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::MATCH_CUST Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::DIST
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::DIST Constraint\n");
	double matchDistConData[] = {1, 0.2};
	TPAT_Constraint matchDistCon(TPAT_Constraint_Tp::DIST, 3, matchDistConData, 2);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(matchDistCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_CR3BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		correctedSet.saveToMat("CR_EM_DIST_corrected.mat");
		finalState = correctedSet.getState(matchDistCon.getID());
		double dist = sqrt(pow(finalState[0] - 1 + sys.getMu(),2) + pow(finalState[1], 2));
		std::cout << "TPAT_Constraint_Tp::DIST Constraint: " << (std::abs(dist - matchDistConData[1]) < 1e-12 ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::DIST Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::MIN_DIST
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::MIN_DIST Constraint\n");
	matchDistConData[1] = 0.4;
	matchDistCon.setData(matchDistConData, 2);
	matchDistCon.setType(TPAT_Constraint_Tp::MIN_DIST);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(matchDistCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_CR3BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		correctedSet.saveToMat("CR_EM_MIN_DIST_corrected.mat");
		finalState = correctedSet.getState(matchDistCon.getID());
		double dist = sqrt(pow(finalState[0] - 1 + sys.getMu(), 2) + pow(finalState[1], 2));
		std::cout << "TPAT_Constraint_Tp::MIN_DIST Constraint: " << (dist >= matchDistConData[1] ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::MIN_DIST Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::MAX_DIST
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::MAX_DIST Constraint\n");
	matchDistConData[1] = 0.15;
	matchDistCon.setData(matchDistConData, 2);
	matchDistCon.setType(TPAT_Constraint_Tp::MAX_DIST);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(matchDistCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_CR3BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		correctedSet.saveToMat("CR_EM_MAX_DIST_corrected.mat");
		finalState = correctedSet.getState(matchDistCon.getID());
		double dist = sqrt(pow(finalState[0] - 1 + sys.getMu(), 2) + pow(finalState[1], 2));
		std::cout << "TPAT_Constraint_Tp::MAX_DIST Constraint: " << (dist <= matchDistConData[1] ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::MAX_DIST Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::MAX_DELTA_V
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::MAX_DELTA_V Constraint\n");
	std::vector<int> dvNodes {3};
	std::vector<double> state = halfLyapNodeset.getState(dvNodes[0]);
	state[3] += 0.01;
	state[4] += 0.1;
	halfLyapNodeset.setState(dvNodes[0], state);	// Perturb the velocity of this state to create a discontinuity
	halfLyapNodeset.allowDV_at(dvNodes);	// Allow the perturbed node to have a delta-v
	double maxDVConData = 0.01;
	TPAT_Constraint dVCon(TPAT_Constraint_Tp::MAX_DELTA_V, 0, &maxDVConData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(dVCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_CR3BP correctedSet(&sys);
		TPAT_MultShoot_Data itData = corrector.multShoot(&halfLyapNodeset, &correctedSet);
		correctedSet.saveToMat("CR_EM_MAX_DV_corrected.mat");
		double totalDV = getTotalDV(&itData);
		std::cout << "TPAT_Constraint_Tp::MAX_DELTA_V Constraint: " << (totalDV <= maxDVConData ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::MAX_DELTA_V Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::DELTA_V
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::DELTA_V Constraint\n");
	maxDVConData = 0.01;
	dVCon.setData(&maxDVConData, 1);
	dVCon.setType(TPAT_Constraint_Tp::DELTA_V);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(dVCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_CR3BP correctedSet(&sys);
		TPAT_MultShoot_Data itData = corrector.multShoot(&halfLyapNodeset, &correctedSet);
		correctedSet.saveToMat("CR_EM_DV_corrected.mat");
		double totalDV = getTotalDV(&itData);
		std::cout << "TPAT_Constraint_Tp::DELTA_V Constraint: " << (std::abs(totalDV - maxDVConData) < 1e-10 ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::DELTA_V Constraint: " << FAIL << std::endl;
	}
	
	// JC
	printColor(BOLDBLACK, "JC Constraint\n");
	double jacobiData = 3.1149;
	TPAT_Constraint jacobiCon(TPAT_Constraint_Tp::JC, 0, &jacobiData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(jacobiCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_CR3BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		double jacobi = correctedSet.getJacobiByIx(0);
		std::cout << "JC Constraint: " << (std::abs(jacobi - jacobiData) < 1e-12 ? PASS : FAIL) << std::endl;
		printf("Desired value = %.4f, actual value = %.4f, diff = %.4e\n", jacobiData, jacobi, std::abs(jacobi - jacobiData));
	}catch(TPAT_Diverge &e){
		std::cout << "JC Constraint: " << FAIL << std::endl;
	}

	// TOF
	printColor(BOLDBLACK, "TOF Constraint\n");
	double tofData = 2.5;
	TPAT_Constraint tofCon(TPAT_Constraint_Tp::TOF, 0, &tofData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(tofCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_CR3BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		double totalTOF = correctedSet.getTotalTOF();
		std::cout << "TOF Constraint: " << (std::abs(totalTOF - tofData) < 1e-12 ? PASS : FAIL) << std::endl;
		printf("Desired TOF = %.4f, actual value = %.4f, diff = %.4e\n", tofData, totalTOF, std::abs(totalTOF - tofData));
	}catch(TPAT_Diverge &e){
		std::cout << "TOF Constraint: " << FAIL << std::endl;
	}

	// APSE
	printColor(BOLDBLACK, "APSE Constraint\n");
	double apseData = 1;
	TPAT_Constraint apseCon(TPAT_Constraint_Tp::APSE, 4, &apseData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(apseCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_CR3BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		correctedSet.saveToMat("CR_EM_APSE_corrected.mat");
		finalState = correctedSet.getState(apseCon.getID());
		const TPAT_Model *model = sys.getModel();
		std::vector<double> primPos = model->getPrimPos(0, &sys);
		double dx = finalState[0] - primPos[apseData*3 + 0];
		double dy = finalState[1] - primPos[apseData*3 + 1];
		double dz = finalState[2] - primPos[apseData*3 + 2];
		double rdot = dx*finalState[3] + dy*finalState[4] + dz*finalState[5];
		std::cout << "APSE Constraint: " << (std::abs(rdot) < 1e-12 ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "APSE Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::SEG_CONT_PV
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::SEG_CONT_PV Constraint\n");

	TPAT_Nodeset_CR3BP forwardArc(ic, &sys, T/2, 4);
	TPAT_Nodeset_CR3BP reverseArc(ic, &sys, -T/2.1, 4);
	TPAT_Nodeset_CR3BP doubleSrcLyap = forwardArc;
	doubleSrcLyap.appendSetAtNode(&reverseArc, 0, 0, 0);
	// doubleSrcLyap.print();
	// doubleSrcLyap.printInChrono();

	double contData[] = {4, 4, NAN, 4, 4, NAN};
	TPAT_Constraint contCon(TPAT_Constraint_Tp::SEG_CONT_PV, 2, contData, 6);
	doubleSrcLyap.addConstraint(contCon);

	try{
		finiteDiff_checkMultShoot(&doubleSrcLyap, corrector);
		try{
			TPAT_Nodeset_CR3BP correctedSet(&sys);
			TPAT_MultShoot_Data it = corrector.multShoot(&doubleSrcLyap, &correctedSet);
			TPAT_Traj forwardTraj = it.propSegs[correctedSet.getSegIx(contCon.getID())];
			TPAT_Traj reverseTraj = it.propSegs[correctedSet.getSegIx(contData[0])];
			std::vector<double> for_lastState = forwardTraj.getStateByIx(-1);
			std::vector<double> rev_lastState = reverseTraj.getStateByIx(-1);
			double sum = 0;
			for(int i = 0; i < 6; i++){
				if(!std::isnan(contData[i]))
					sum += pow(for_lastState[i] - rev_lastState[i], 2);
			}
			double dist = std::sqrt(sum);
			std::cout << "TPAT_Constraint_Tp::SEG_CONT_PV Constraint: " << (std::abs(dist) < 1e-12 ? PASS : FAIL) << std::endl;
		}catch(TPAT_Diverge &e){
			std::cout << "TPAT_Constraint_Tp::SEG_CONT_PV Constraint: " << FAIL << std::endl;
		}

		if(equalArcTime)
			std::cout << "TPAT_Constraint_Tp::SEG_CONT_PV Constraint (equalArcTime): " << FAIL << std::endl;
	}catch(TPAT_Exception &e){
		std::cout << "TPAT_Constraint_Tp::SEG_CONT_PV Constraint (equalArcTime): " << PASS << std::endl;
	}
}//====================================================

void testBCR4BPCons(bool equalArcTime){
	TPAT_Sys_Data_BC4BP sys("sun", "earth", "moon");
	double ic[] = {-0.745230328320519, 7.22625684942683e-04, 7.45549413286038e-05, -7.30710697247992e-06, -0.0148897145134465, -1.23266135281459e-06};
	double T = 313;	// SE L1 Period
	TPAT_Nodeset_BC4BP halfLyapNodeset(ic, &sys, 0, T, 6);	// Create a nodeset
	std::vector<double> initState, finalState;
	TPAT_Nodeset_BC4BP correctedSet(&sys);

	TPAT_Correction_Engine corrector;
	corrector.setEqualArcTime(equalArcTime);
	printColor(BOLDBLACK, "Testing BCR4BP Multiple Shooting Constraints\n");

	// TPAT_Constraint_Tp::STATE
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::STATE Constraint\n");
	double stateConData[] = {-0.77, 0.5, NAN, NAN, NAN, NAN};
	TPAT_Constraint stateCon(TPAT_Constraint_Tp::STATE, 4, stateConData, 6);
	halfLyapNodeset.addConstraint(stateCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_BC4BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		finalState = correctedSet.getState(stateCon.getID());
		std::cout << "TPAT_Constraint_Tp::STATE Constraint: " << (stateDiffBelowTol(finalState, stateConData, 1e-12) ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::STATE Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::MATCH_ALL
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::MATCH_ALL Constraint\n");
	double matchAllConData = 0;
	TPAT_Constraint matchAllCon(TPAT_Constraint_Tp::MATCH_ALL, 4, &matchAllConData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(matchAllCon);
	// halfLyapNodeset.print();
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_BC4BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		finalState = correctedSet.getState(matchAllCon.getID());
		initState = correctedSet.getStateByIx(0);
		std::cout << "TPAT_Constraint_Tp::MATCH_ALL Constraint: " << (stateDiffBelowTol(finalState, initState, 1e-12) ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::MATCH_ALL Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::MATCH_CUST
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::MATCH_CUST Constraint\n");
	double matchCustConData[] = {0,0,NAN,NAN,NAN,NAN};
	TPAT_Constraint matchCustCon(TPAT_Constraint_Tp::MATCH_CUST, 4, matchCustConData, 6);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(matchCustCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_BC4BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		finalState = correctedSet.getState(matchCustCon.getID());
		initState = correctedSet.getStateByIx(0);
		finalState.erase(finalState.begin()+2, finalState.end());	// Erase entries 2 through 5; we're only comparing the first two
		initState.erase(initState.begin()+2, initState.end());
		std::cout << "TPAT_Constraint_Tp::MATCH_CUST Constraint: " << (stateDiffBelowTol(finalState, initState, 1e-12) ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::MATCH_CUST Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::DIST
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::DIST Constraint\n");
	double matchDistConData[] = {1, 1.0};
	TPAT_Constraint matchDistCon(TPAT_Constraint_Tp::DIST, 3, matchDistConData, 2);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(matchDistCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_BC4BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		finalState = correctedSet.getState(matchDistCon.getID());
		std::vector<double> primPos = sys.getModel()->getPrimPos(correctedSet.getEpoch(matchDistCon.getID()), &sys);
		double dist = sqrt(pow(finalState[0] - primPos[3] ,2) + pow(finalState[1] - primPos[4], 2) + pow(finalState[2] - primPos[5], 2));
		std::cout << "TPAT_Constraint_Tp::DIST Constraint: " << (std::abs(dist - matchDistConData[1]) < 1e-12 ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::DIST Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::MIN_DIST
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::MIN_DIST Constraint\n");
	matchDistConData[1] = 1.1;
	matchDistCon.setData(matchDistConData, 2);
	matchDistCon.setType(TPAT_Constraint_Tp::MIN_DIST);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(matchDistCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_BC4BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		finalState = correctedSet.getState(matchDistCon.getID());
		std::vector<double> primPos = sys.getModel()->getPrimPos(correctedSet.getEpoch(matchDistCon.getID()), &sys);
		double dist = sqrt(pow(finalState[0] - primPos[3] ,2) + pow(finalState[1] - primPos[4], 2) + pow(finalState[2] - primPos[5], 2));
		std::cout << "TPAT_Constraint_Tp::MIN_DIST Constraint: " << (dist >= matchDistConData[1] ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::MIN_DIST Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::MAX_DIST
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::MAX_DIST Constraint\n");
	matchDistConData[1] = 0.9;
	// matchDistConData[1] = 2;
	matchDistCon.setData(matchDistConData, 2);
	matchDistCon.setType(TPAT_Constraint_Tp::MAX_DIST);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(matchDistCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_BC4BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		finalState = correctedSet.getState(matchDistCon.getID());
		std::vector<double> primPos = sys.getModel()->getPrimPos(correctedSet.getEpoch(matchDistCon.getID()), &sys);
		double dist = sqrt(pow(finalState[0] - primPos[3] ,2) + pow(finalState[1] - primPos[4], 2) + pow(finalState[2] - primPos[5], 2));
		std::cout << "TPAT_Constraint_Tp::MAX_DIST Constraint: " << (dist <= matchDistConData[1] ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::MAX_DIST Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::MAX_DELTA_V
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::MAX_DELTA_V Constraint\n");
	std::vector<double> state = halfLyapNodeset.getStateByIx(3);
	state[3] += 0.01;
	state[4] += 0.1;
	state[5] += 0.001;
	halfLyapNodeset.setState(3, state);	// Perturb the velocity of this state to create a discontinuity
	std::vector<int> dvSegs {2};
	halfLyapNodeset.allowDV_at(dvSegs);	// Allow the perturbed node to have a delta-v
	double maxDVConData = 0.03;
	TPAT_Constraint dVCon(TPAT_Constraint_Tp::MAX_DELTA_V, 0, &maxDVConData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(dVCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_BC4BP correctedSet(&sys);
		TPAT_MultShoot_Data itData = corrector.multShoot(&halfLyapNodeset, &correctedSet);
		double totalDV = getTotalDV(&itData);
		std::cout << "TPAT_Constraint_Tp::MAX_DELTA_V Constraint: " << (totalDV <= maxDVConData ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::MAX_DELTA_V Constraint: " << FAIL << std::endl;
	}

	// TPAT_Constraint_Tp::DELTA_V
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::DELTA_V Constraint\n");
	maxDVConData = 0.02*0.02;
	dVCon.setData(&maxDVConData, 1);
	dVCon.setType(TPAT_Constraint_Tp::DELTA_V);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(dVCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_BC4BP correctedSet(&sys);
		TPAT_MultShoot_Data itData = corrector.multShoot(&halfLyapNodeset, &correctedSet);
		double totalDV = getTotalDV(&itData);
		std::cout << "TPAT_Constraint_Tp::DELTA_V Constraint: " << (std::abs(totalDV - maxDVConData) < 1e-12 ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TPAT_Constraint_Tp::DELTA_V Constraint: " << FAIL << std::endl;
	}

	// TOF
	printColor(BOLDBLACK, "TOF Constraint\n");
	double tofData = 340.0;
	TPAT_Constraint tofCon(TPAT_Constraint_Tp::TOF, 0, &tofData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(tofCon);
	// halfLyapNodeset.print();
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_BC4BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		double totalTOF = correctedSet.getTotalTOF();
		std::cout << "TOF Constraint: " << (std::abs(totalTOF - tofData) < 1e-12 ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "TOF Constraint: " << FAIL << std::endl;
	}

	// APSE
	printColor(BOLDBLACK, "APSE Constraint\n");
	double apseData = 2;
	TPAT_Constraint apseCon(TPAT_Constraint_Tp::APSE, 4, &apseData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(apseCon);
	finiteDiff_checkMultShoot(&halfLyapNodeset, corrector);
	try{
		TPAT_Nodeset_BC4BP correctedSet(&sys);
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
		finalState = correctedSet.getState(apseCon.getID());
		const TPAT_Model *model = sys.getModel();
		double rdot = model->getRDot(apseData, correctedSet.getEpoch(apseCon.getID()), &(finalState[0]), &sys);
		std::cout << "APSE Constraint: " << (std::abs(rdot) < 1e-12 ? PASS : FAIL) << std::endl;
		printf("rdot = %.4e\n", rdot);
	}catch(TPAT_Diverge &e){
		std::cout << "APSE Constraint: " << FAIL << std::endl;
	}


	// Saddle Point, Exact
	printColor(BOLDBLACK, "SP Constraint\n");
	double IC[] = {-0.71200455, 0.16675922, 0.02755461, 0.01186449, -0.00004723, -0.0010737};
	double t0 = 10207.19;
	// double t0 = 0;
	double tof = 51.32;
	TPAT_Nodeset_BC4BP nodes0(IC, &sys, t0, tof, 5);
	corrector.setTol(1e-11);
	// nodes0.saveToMat("SP_TestCase.mat");

	double spData = 0;
	TPAT_Constraint spCon(TPAT_Constraint_Tp::SP, 2, &spData, 1);
	nodes0.addConstraint(spCon);
	finiteDiff_checkMultShoot(&nodes0, corrector);
	try{
		TPAT_Nodeset_BC4BP correctedSet(&sys);
		corrector.multShoot(&nodes0, &correctedSet);
		finalState = correctedSet.getState(spCon.getID());
		finalState.erase(finalState.begin()+3, finalState.end());

		Eigen::Vector3d spPos = bcr4bpr_getSPLoc(&sys, correctedSet.getEpoch(spCon.getID()));
		double diff = sqrt(pow(spPos(0) - finalState[0], 2) + pow(spPos(1) - finalState[1], 2) + pow(spPos(2) - finalState[2], 2));
		std::cout << "SP Constraint: " << (diff < 1e-10 ? PASS : FAIL) << std::endl;
	}catch(TPAT_Diverge &e){
		std::cout << "SP Constraint: " << FAIL << std::endl;
	}
	
	// // Saddle Point, Range
	// printColor(BOLDBLACK, "SP Range Constraint\n");
	// double maxR = 40;	// km
	// double maxA = (5.358e-8*maxR - 2.3313e-8)/1000*sys.getCharT()*sys.getCharT()/sys.getCharL();
	// printf("Maximum Accel  = %e (non-dim)\n", maxA);
	// TPAT_Constraint spConRange(TPAT_Constraint_Tp::SP_RANGE, 2, &maxA, 1);
	// nodes0.clearAllConstraints();
	// nodes0.addConstraint(spConRange);
	// finiteDiff_checkMultShoot(&nodes0);
	// try{
	// 	TPAT_Nodeset_BC4BP correctedSet(&sys);
	// 	corrector.multShoot(&nodes0, &correctedSet);
	// 	finalState = correctedSet.getState(spCon.getID());

	// 	// Compute the state derivative at the node to get acceleration
	// 	double stateDot[6];
	// 	TPAT_Model_BC4BP::simpleEOMs(correctedSet.getEpoch(spCon.getID()), &(finalState[0]), stateDot, &sys);

	// 	// Adjust derivative from EOMs to remove terms from rotating frame
	// 	stateDot[3] += -2*sys.getK()*finalState[4] - sys.getK()*sys.getK()*finalState[0] - sys.getK()*sys.getK()*(1/sys.getK() - sys.getMu());
	// 	stateDot[4] += 2*sys.getK()*finalState[3] - sys.getK()*sys.getK()*finalState[1];
	// 	double accel = sqrt(stateDot[3]*stateDot[3] + stateDot[4]*stateDot[4] + stateDot[5]*stateDot[5]);

	// 	// Compute the location of the SP
	// 	finalState.erase(finalState.begin()+3, finalState.end());
	// 	Eigen::Vector3d spPos = bcr4bpr_getSPLoc(&sys, correctedSet.getEpoch(spCon.getID()));
	// 	double diff = sqrt(pow(spPos(0) - finalState[0], 2) + pow(spPos(1) - finalState[1], 2) + pow(spPos(2) - finalState[2], 2));

	// 	std::cout << "SP Range Constraint: " << (accel <= maxA ? PASS : FAIL) << std::endl;
	// 	printf("Accel = %e, Max Accel = %e\n", accel, maxA);
	// 	printf("Distance from SP = %.4f km\n", diff*sys.getCharL());
	// }catch(TPAT_Diverge &e){
	// 	std::cout << "SP Range Constraint: " << FAIL << std::endl;
	// }

	// //Saddle Point, Approx. Location, Distance Range
	// printColor(BOLDBLACK, "SP TPAT_Constraint_Tp::DIST Constraint\n");
	// MatrixXRd coeff = bcr4bpr_spLoc_polyFit(&sys, nodes0.getEpochByIx(2));
	// toCSV(coeff, "PolyFitCoeff.csv");
	// printf("Epoch at SP Poly Fit is %.4f\n", nodes0.getEpochByIx(2));
	// double maxDist = 1000.0/sys.getCharL();
	// double spDistConData[] = {maxDist, coeff(0,0), coeff(1,0), coeff(2,0), coeff(0,1), coeff(1,1),
	// 	coeff(2,1), coeff(0,2), coeff(1,2), coeff(2,2)};
	// TPAT_Constraint spDistCon(TPAT_Constraint_Tp::SP_DIST, 2, spDistConData, 10);
	// nodes0.clearAllConstraints();
	// nodes0.addConstraint(spDistCon);
	// finiteDiff_checkMultShoot(&nodes0);
	// try{
	// 	TPAT_Nodeset_BC4BP correctedSet(&sys);
	// 	corrector.multShoot(&nodes0, &correctedSet);
	// 	finalState = correctedSet.getState(spDistCon.getID());
	// 	double T = correctedSet.getEpoch(spDistCon.getID());

	// 	// Compute the approximate location of SP
	// 	Eigen::RowVector3d indVarMat;
	// 	indVarMat << T*T, T, 1;
	// 	Eigen::RowVector3d spPos_approx = indVarMat*coeff;

	// 	// Compute the location of the SP
	// 	Eigen::RowVector3d spPos = bcr4bpr_getSPLoc(&sys, T);
		
	// 	// Compute errors
	// 	double diff = sqrt(pow(spPos_approx(0) - finalState[0], 2) + pow(spPos_approx(1) - finalState[1], 2) + pow(spPos_approx(2) - finalState[2], 2));
	// 	double epochShift = nodes0.getEpoch(spDistCon.getID()) - T;
	// 	Eigen::Vector3d approxDiff = spPos_approx - spPos;

	// 	std::cout << "SP Dist Constraint: " << (std::abs(diff - maxDist) < 1e-12 ? PASS : FAIL) << std::endl;
	// 	printf("Node distance from approximated SP = %.4f km\n", diff*sys.getCharL());
	// 	printf("Approximated SP error = %.4f km\n", approxDiff.norm()*sys.getCharL());
	// 	printf("Epoch shifted by %.4e sec = %.4e days\n", epochShift*sys.getCharT(), epochShift*sys.getCharT()/24/3600);
	// }catch(TPAT_Diverge &e){
	// 	std::cout << "SP Dist Constraint: " << FAIL << std::endl;
	// }

	// // Saddle Point, Approx. Location, Distance Range
	// printColor(BOLDBLACK, "SP MAX TPAT_Constraint_Tp::DIST Constraint\n");
	// TPAT_Constraint spMaxDistCon(TPAT_Constraint_Tp::SP_MAX_DIST, 2, spDistConData, 10);	// Same dist as previous, but with inequality now
	// nodes0.clearAllConstraints();
	// nodes0.addConstraint(spMaxDistCon);
	// finiteDiff_checkMultShoot(&nodes0);
	// try{
	// 	TPAT_Nodeset_BC4BP correctedSet(&sys);
	// 	corrector.multShoot(&nodes0, &correctedSet);
	// 	finalState = correctedSet.getState(spMaxDistCon.getID());
	// 	double T = correctedSet.getEpoch(spMaxDistCon.getID());

	// 	// Compute the approximate location of SP
	// 	Eigen::RowVector3d indVarMat;
	// 	indVarMat << T*T, T, 1;
	// 	Eigen::RowVector3d spPos_approx = indVarMat*coeff;

	// 	// Compute the location of the SP
	// 	Eigen::RowVector3d spPos = bcr4bpr_getSPLoc(&sys, T);
		
	// 	// Compute errors
	// 	double diff = sqrt(pow(spPos_approx(0) - finalState[0], 2) + pow(spPos_approx(1) - finalState[1], 2) + pow(spPos_approx(2) - finalState[2], 2));
	// 	double epochShift = nodes0.getEpoch(spMaxDistCon.getID()) - T;
	// 	Eigen::Vector3d approxDiff = spPos_approx - spPos;

	// 	std::cout << "SP Max Dist Constraint: " << (std::abs(diff) <= maxDist ? PASS : FAIL) << std::endl;
	// 	printf("Node distance from approximated SP = %.4f km\n", diff*sys.getCharL());
	// 	printf("Approximated SP error = %.4f km\n", approxDiff.norm()*sys.getCharL());
	// 	printf("Epoch shifted by %.4e sec = %.4e days\n", epochShift*sys.getCharT(), epochShift*sys.getCharT()/24/3600);
	// }catch(TPAT_Diverge &e){
	// 	std::cout << "SP Max Dist Constraint: " << FAIL << std::endl;
	// }

	// TPAT_Constraint_Tp::SEG_CONT_PV
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::SEG_CONT_PV Constraint\n");

	TPAT_Nodeset_BC4BP forwardArc(ic, &sys, 0, T/2, 4);
	forwardArc.deleteNode(3);
	TPAT_Nodeset_BC4BP reverseArc(ic, &sys, T, -T/2, 4);
	reverseArc.deleteNode(3);
	TPAT_Nodeset_BC4BP doubleSrcLyap = forwardArc;
	doubleSrcLyap.concatArcset(&reverseArc);
	// doubleSrcLyap.appendSetAtNode(&reverseArc, 0, 0, 0);

	double contData[] = {4, 4, NAN, 4, 4, NAN};
	TPAT_Constraint contCon(TPAT_Constraint_Tp::SEG_CONT_PV, 2, contData, 6);
	doubleSrcLyap.addConstraint(contCon);

	// doubleSrcLyap.print();
	// doubleSrcLyap.printInChrono();

	try{
		finiteDiff_checkMultShoot(&doubleSrcLyap, corrector);
		try{
			TPAT_Nodeset_BC4BP correctedSet(&sys);
			TPAT_MultShoot_Data it = corrector.multShoot(&doubleSrcLyap, &correctedSet);
			TPAT_Traj forwardTraj = it.propSegs[correctedSet.getSegIx(contCon.getID())];
			TPAT_Traj reverseTraj = it.propSegs[correctedSet.getSegIx(contData[0])];
			std::vector<double> for_lastState = forwardTraj.getStateByIx(-1);
			std::vector<double> rev_lastState = reverseTraj.getStateByIx(-1);
			double sum = 0;
			for(int i = 0; i < 6; i++){
				if(!std::isnan(contData[i]))
					sum += pow(for_lastState[i] - rev_lastState[i], 2);
			}
			double dist = std::sqrt(sum);
			std::cout << "TPAT_Constraint_Tp::SEG_CONT_PV Constraint: " << (std::abs(dist) < 1e-12 ? PASS : FAIL) << std::endl;
		}catch(TPAT_Diverge &e){
			std::cout << "TPAT_Constraint_Tp::SEG_CONT_PV Constraint: " << FAIL << std::endl;
		}

		if(equalArcTime)
			std::cout << "TPAT_Constraint_Tp::SEG_CONT_PV Constraint (equalArcTime): " << FAIL << std::endl;
	}catch(TPAT_Exception &e){
		std::cout << "TPAT_Constraint_Tp::SEG_CONT_PV Constraint (equalArcTime): " << PASS << std::endl;
	}

	// TPAT_Constraint_Tp::SEG_CONT_EX
	printColor(BOLDBLACK, "TPAT_Constraint_Tp::SEG_CONT_EX Constraint\n");
	doubleSrcLyap.clearAllConstraints();

	double contExData[] = {4, 0};
	TPAT_Constraint extraContCon(TPAT_Constraint_Tp::SEG_CONT_EX, 2, contExData, 2);
	doubleSrcLyap.addConstraint(extraContCon);

	try{
		finiteDiff_checkMultShoot(&doubleSrcLyap, corrector);
		try{
			TPAT_Nodeset_BC4BP correctedSet(&sys);
			TPAT_MultShoot_Data it = corrector.multShoot(&doubleSrcLyap, &correctedSet);
			TPAT_Traj forwardTraj = it.propSegs[correctedSet.getSegIx(contCon.getID())];
			TPAT_Traj reverseTraj = it.propSegs[correctedSet.getSegIx(contData[0])];
			std::cout << "TPAT_Constraint_Tp::SEG_CONT_EX Constraint: " << 
				(std::abs(forwardTraj.getEpochByIx(-1) - reverseTraj.getEpochByIx(-1)) < 1e-12 ? PASS : FAIL) <<
				std::endl;
		}catch(TPAT_Diverge &e){
			std::cout << "TPAT_Constraint_Tp::SEG_CONT_EX Constraint: " << FAIL << std::endl;
		}

		if(equalArcTime)
			std::cout << "TPAT_Constraint_Tp::SEG_CONT_PV Constraint (equalArcTime): " << FAIL << std::endl;
	}catch(TPAT_Exception &e){
		std::cout << "TPAT_Constraint_Tp::SEG_CONT_PV Constraint (equalArcTime): " << PASS << std::endl;
	}
}//====================================================

/**
 *  @brief Test all constraint types available to ensure they converge correctly
 */
int main(void){
	// First, run with equalArcTime = false
	testCR3BP_SE_Cons(false);
	testCR3BP_EM_Cons(false);
	testBCR4BPCons(false);

	testCR3BP_SE_Cons(true);
	testCR3BP_EM_Cons(true);
	testBCR4BPCons(true);
	return EXIT_SUCCESS;
}
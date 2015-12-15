#include "tpat_ascii_output.hpp"
#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
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

/**
 *  @brief Test all constraint types available to ensure they converge correctly
 */
int main(void){

	tpat_sys_data_cr3bp sys("earth", "moon");
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period
	tpat_nodeset_cr3bp halfLyapNodeset(ic, &sys, T/1.25, 8);	// Create a nodeset with three nodes

	std::vector<double> initState, finalState;
	tpat_nodeset_cr3bp correctedSet(&sys);

	tpat_correction_engine corrector;
	printColor(BOLDBLACK, "Testing CR3BP Multiple Shooting Constraints\n");

	// STATE
	double stateConData[] = {0.9, 0.1, NAN, NAN, NAN, NAN};
	tpat_constraint stateCon(tpat_constraint::STATE, 7, stateConData, 6);
	halfLyapNodeset.addConstraint(stateCon);
	try{
		corrector.multShoot(&halfLyapNodeset);
		correctedSet = corrector.getCR3BP_Output();
		finalState = correctedSet.getState(stateCon.getNode());
		std::cout << "STATE Constraint: " << (stateDiffBelowTol(finalState, stateConData, 1e-12) ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "STATE Constraint: " << FAIL << std::endl;
	}

	// MATCH_ALL
	double matchAllConData = 0;
	tpat_constraint matchAllCon(tpat_constraint::MATCH_ALL, 7, &matchAllConData, 1);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(matchAllCon);
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
	double matchCustConData[] = {0,0,NAN,NAN,NAN,NAN};
	tpat_constraint matchCustCon(tpat_constraint::MATCH_CUST, 7, matchCustConData, 6);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(matchCustCon);
	try{
		corrector.multShoot(&halfLyapNodeset);
		correctedSet = corrector.getCR3BP_Output();
		finalState = correctedSet.getState(matchCustCon.getNode());
		initState = correctedSet.getState(0);
		finalState.erase(finalState.begin()+2, finalState.begin()+5);	// Erase entries 2 through 5; we're only comparing the first two
		initState.erase(initState.begin()+2, initState.begin()+5);
		std::cout << "MATCH_CUST Constraint: " << (stateDiffBelowTol(finalState, initState, 1e-12) ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "MATCH_CUST Constraint: " << FAIL << std::endl;
	}

	// DIST
	double matchDistConData[] = {1, 0.2};
	tpat_constraint matchDistCon(tpat_constraint::DIST, 5, matchDistConData, 2);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(matchDistCon);
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
	matchDistConData[1] = 0.1;
	matchDistCon.setData(matchDistConData, 2);
	matchDistCon.setType(tpat_constraint::MIN_DIST);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(matchDistCon);
	try{
		corrector.multShoot(&halfLyapNodeset);
		correctedSet = corrector.getCR3BP_Output();
		finalState = correctedSet.getState(matchDistCon.getNode());
		double dist = sqrt(pow(finalState[0] - + sys.getMu(), 2) + pow(finalState[1], 2));
		std::cout << "MIN_DIST Constraint: " << (dist >= matchDistConData[1] ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "MIN_DIST Constraint: " << FAIL << std::endl;
	}

	// MAX_DIST
	matchDistConData[1] = 0.3;
	matchDistCon.setData(matchDistConData, 2);
	matchDistCon.setType(tpat_constraint::MAX_DIST);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(matchDistCon);
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
	std::vector<double> state = halfLyapNodeset.getState(4);
	state[3] += 0.01;
	state[4] += 0.01;
	halfLyapNodeset.setState(4, state);
	std::vector<int> dvNodes {4};
	halfLyapNodeset.setVelConNodes_allBut(dvNodes);
	double maxDVConData = 0.01;
	tpat_constraint dVCon(tpat_constraint::MAX_DELTA_V, 0, &maxDVConData, 1);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(dVCon);
	// halfLyapNodeset.saveToMat("tempNodes.mat");
	// halfLyapNodeset.print();
	double totalDV = 0;
	try{
		iterationData itData = corrector.multShoot(&halfLyapNodeset);
		for(size_t i = 0; i < itData.deltaVs.size()/3; i+=3){
			totalDV += sqrt(pow(itData.deltaVs[i], 2) + pow(itData.deltaVs[i+1], 2) +
				pow(itData.deltaVs[i+2], 2));
		}
		std::cout << "MAX_DELTA_V Constraint: " << (totalDV <= maxDVConData ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "MAX_DELTA_V Constraint: " << FAIL << std::endl;
		printf("  Total DV = %.4e\n", totalDV);
	}

	// DELTA_V
	maxDVConData = 0.02;
	dVCon.setData(&maxDVConData, 1);
	dVCon.setType(tpat_constraint::DELTA_V);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(dVCon);
	totalDV = 0;
	try{
		iterationData itData = corrector.multShoot(&halfLyapNodeset);
		for(size_t i = 0; i < itData.deltaVs.size()/3; i+=3){
			totalDV += sqrt(pow(itData.deltaVs[i], 2) + pow(itData.deltaVs[i+1], 2) +
				pow(itData.deltaVs[i+2], 2));
		}
		std::cout << "DELTA_V Constraint: " << (std::abs(totalDV - maxDVConData) < 1e-12 ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "DELTA_V Constraint: " << FAIL << std::endl;
		printf("  Total DV = %.4e\n", totalDV);
	}

	// SP
	// SP_RANGE
	
	// JC
	double jacobiData = 3.1149;
	tpat_constraint jacobiCon(tpat_constraint::JC, 0, &jacobiData, 1);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(jacobiCon);
	try{
		corrector.multShoot(&halfLyapNodeset);
		correctedSet = corrector.getCR3BP_Output();
		double jacobi = correctedSet.getJacobi(0);
		std::cout << "JC Constraint: " << (std::abs(jacobi - jacobiData) < 1e-12 ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "JC Constraint: " << FAIL << std::endl;
	}

	// TOF
	double tofData = 2.5;
	tpat_constraint tofCon(tpat_constraint::TOF, 0, &tofData, 1);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(tofCon);
	try{
		corrector.multShoot(&halfLyapNodeset);
		correctedSet = corrector.getCR3BP_Output();
		double totalTOF = correctedSet.getTotalTOF();
		std::cout << "TOF Constraint: " << (std::abs(totalTOF - tofData) < 1e-12 ? PASS : FAIL) << std::endl;
	}catch(tpat_diverge &e){
		std::cout << "TOF Constraint: " << FAIL << std::endl;
	}

	// APSE
	double apseData = 1;
	tpat_constraint apseCon(tpat_constraint::APSE, 7, &apseData, 1);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(apseCon);
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
	return EXIT_SUCCESS;
}
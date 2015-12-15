#include "tpat_ascii_output.hpp"
#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
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

	tpat_correction_engine corrector;
	printColor(BOLDBLACK, "Testing CR3BP Multiple Shooting Constraints\n");

	// STATE
	double stateConData[] = {0.9, 0.1, NAN, NAN, NAN, NAN};
	tpat_constraint stateCon(tpat_constraint::STATE, 7, stateConData, 6);
	halfLyapNodeset.addConstraint(stateCon);
	corrector.multShoot(&halfLyapNodeset);
	tpat_nodeset_cr3bp correctedSet = corrector.getCR3BP_Output();
	std::vector<double> finalState = correctedSet.getState(-1);
	std::cout << "STATE Constraint: " << (stateDiffBelowTol(finalState, stateConData, 1e-12) ? PASS : FAIL) << std::endl;

	// MATCH_ALL
	double matchAllConData = 0;
	tpat_constraint matchAllCon(tpat_constraint::MATCH_ALL, 7, &matchAllConData, 1);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(matchAllCon);
	corrector.multShoot(&halfLyapNodeset);
	correctedSet = corrector.getCR3BP_Output();
	finalState = correctedSet.getState(-1);
	std::vector<double> initState = correctedSet.getState(0);
	std::cout << "MATCH_ALL Constraint: " << (stateDiffBelowTol(finalState, initState, 1e-10) ? PASS : FAIL) << std::endl;

	// MATCH_CUST
	double matchCustConData[] = {0,0,NAN,NAN,NAN,NAN};
	tpat_constraint matchCustCon(tpat_constraint::MATCH_CUST, 7, matchCustConData, 6);
	halfLyapNodeset.clearConstraints();
	halfLyapNodeset.addConstraint(matchCustCon);
	corrector.multShoot(&halfLyapNodeset);
	correctedSet = corrector.getCR3BP_Output();
	finalState = correctedSet.getState(-1);
	initState = correctedSet.getState(0);
	finalState.erase(finalState.begin()+2, finalState.begin()+5);	// Erase entries 2 through 5; we're only comparing the first two
	initState.erase(initState.begin()+2, initState.begin()+5);

	std::cout << "MATCH_CUST Constraint: " << (stateDiffBelowTol(finalState, initState, 1e-10) ? PASS : FAIL) << std::endl;

	// DIST
	// MIN_DIST
	// MAX_DIST
	// MAX_DELTA_V
	// DELTA_V
	// SP
	// SP_RANGE
	// JC
	// TOF
	// APSE
	return EXIT_SUCCESS;
}
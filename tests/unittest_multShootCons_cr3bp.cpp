#define BOOST_TEST_MODULE MultipleShootingConstraints_CR3BP

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <vector>

#include "Arcset_cr3bp.hpp"
#include "Calculations.hpp"
#include "Constraint.hpp"
#include "MultShootEngine.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "SimEngine.hpp"
#include "SysData_cr3bp.hpp"
#include "Utilities.hpp"

using namespace astrohelion;

double se_lyap_ic[] = {0.993986593871357, 0, 0, 0, -0.022325793891591, 0};	// SE L1
double em_lyap_ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0};	// EM L1
double se_lyap_T = 3.293141367224790;
double em_lyap_T = 3.02796323553149;

bool stateDiffBelowTol(std::vector<double>, double*, double);
bool stateDiffBelowTol(std::vector<double>, std::vector<double>, double);

/**
 *  \brief Determine if the difference bewteen two state vectors is less than 
 *  the desired tolerance
 * 
 *  \param data vector of state values from the corrections process
 *  \param correct array of state values from the constraint
 *  \param tol desired numerical tolerance
 *  \return whether or not <tt>data</tt> and <tt>correct</tt> are equal
 *  within the desired tolerance
 */
bool stateDiffBelowTol(std::vector<double> data, double *correct, double tol){
	double sum = 0;
	for(unsigned int i = 0; i < data.size(); i++){
		if(!std::isnan(correct[i]))
			sum += pow(data[i] - correct[i], 2);
	}
	return sqrt(sum) < tol;
}//====================================================

bool stateDiffBelowTol(std::vector<double> data, std::vector<double> correct, double tol){
	return stateDiffBelowTol(data, &(correct[0]), tol);
}//====================================================

//************************************************************
//* CR3BP Sun-Earth Constraints
//************************************************************

BOOST_AUTO_TEST_SUITE(CR3BP_SunEarth)

BOOST_AUTO_TEST_CASE(CR3BP_SE_STATE){
	SysData_cr3bp sys("sun", "earth");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(se_lyap_T < 0);
	sim.runSim_manyNodes(se_lyap_ic, se_lyap_T/1.25, 8, &halfLyapNodeset);

	std::vector<double> initState, finalState;

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);

	double stateConData[] = {0.9934, 0.001, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 7, stateConData, 6);
	halfLyapNodeset.addConstraint(stateCon);
	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));

	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_SE_STATE_EQUAL_ARC){
	SysData_cr3bp sys("sun", "earth");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(se_lyap_T < 0);
	sim.runSim_manyNodes(se_lyap_ic, se_lyap_T/1.25, 8, &halfLyapNodeset);
	std::vector<double> initState, finalState;

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);

	double stateConData[] = {0.9934, 0.001, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 7, stateConData, 6);
	halfLyapNodeset.addConstraint(stateCon);
	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));

	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_AUTO_TEST_SUITE_END()

//************************************************************
//* CR3BP Earth-Moon Constraints
//************************************************************

BOOST_AUTO_TEST_SUITE(CR3BP_EarthMoon)

BOOST_AUTO_TEST_CASE(CR3BP_EM_STATE){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);

	double stateConData[] = {0.9, 0.1, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 4, stateConData, 6);
	halfLyapNodeset.addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_STATE_EQUAL_ARC){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);

	double stateConData[] = {0.9, 0.1, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 4, stateConData, 6);
	halfLyapNodeset.addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_STATE_ENDSEG){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T/2.1, 3, &halfLyapNodeset);
	halfLyapNodeset.deleteNode(2);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);

	double stateConData[] = {NAN, 0, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::ENDSEG_STATE, 1, stateConData, 6);
	halfLyapNodeset.addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> fullState = correctedSet.getSegByIx(stateCon.getID()).getStateByRow(-1);
	std::vector<double> finalState(fullState.begin(), fullState.begin()+6);
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_STATE_ENDSEG_EQUAL_ARC){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T/2.1, 3, &halfLyapNodeset);
	halfLyapNodeset.deleteNode(2);

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);

	double stateConData[] = {NAN, 0, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::ENDSEG_STATE, 1, stateConData, 6);
	halfLyapNodeset.addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> fullState = correctedSet.getSegByIx(stateCon.getID()).getStateByRow(-1);
	std::vector<double> finalState(fullState.begin(), fullState.begin()+6);
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_STATE_RMSTATE){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T/2.1, 4, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);

	double stateConData[] = {NAN, 0, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 3, stateConData, 6);
	halfLyapNodeset.addConstraint(stateCon);

	// Remove the initial state from the free variable vector
	Constraint rmState(Constraint_tp::RM_STATE, 0, nullptr, 0);
	halfLyapNodeset.addConstraint(rmState);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> initState = correctedSet.getState(rmState.getID());
	std::vector<double> finalState = correctedSet.getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
	BOOST_CHECK(stateDiffBelowTol(initState, halfLyapNodeset.getState(rmState.getID()), 1e-12));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_MATCH_ALL){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);

	// Constraint_tp::MATCH_ALL
	double matchAllConData = 0;
	Constraint matchAllCon(Constraint_tp::MATCH_ALL, 4, &matchAllConData, 1);
	halfLyapNodeset.addConstraint(matchAllCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(matchAllCon.getID());
	std::vector<double> initState = correctedSet.getStateByIx(0);
	BOOST_CHECK(stateDiffBelowTol(finalState, initState, 1e-12));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_MATCH_ALL_EQUAL_ARC){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);

	// Constraint_tp::MATCH_ALL
	double matchAllConData = 0;
	Constraint matchAllCon(Constraint_tp::MATCH_ALL, 4, &matchAllConData, 1);
	halfLyapNodeset.addConstraint(matchAllCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(matchAllCon.getID());
	std::vector<double> initState = correctedSet.getStateByIx(0);
	BOOST_CHECK(stateDiffBelowTol(finalState, initState, 1e-12));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_MATCH_CUST){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);

	// Constraint_tp::MATCH_CUST
	double matchCustConData[] = {0,0,NAN,NAN,NAN,NAN};
	Constraint matchCustCon(Constraint_tp::MATCH_CUST, 4, matchCustConData, 6);
	halfLyapNodeset.addConstraint(matchCustCon);
	
	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(matchCustCon.getID());
	std::vector<double> initState = correctedSet.getStateByIx(0);
	finalState.erase(finalState.begin()+2, finalState.end());	// Erase entries 2 through 5; we're only comparing the first two
	initState.erase(initState.begin()+2, initState.end());

	BOOST_CHECK(stateDiffBelowTol(finalState, initState, 1e-12));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_MATCH_CUST_EQUAL_ARC){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);

	// Constraint_tp::MATCH_CUST
	double matchCustConData[] = {0,0,NAN,NAN,NAN,NAN};
	Constraint matchCustCon(Constraint_tp::MATCH_CUST, 4, matchCustConData, 6);
	halfLyapNodeset.addConstraint(matchCustCon);
	
	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(matchCustCon.getID());
	std::vector<double> initState = correctedSet.getStateByIx(0);
	finalState.erase(finalState.begin()+2, finalState.end());	// Erase entries 2 through 5; we're only comparing the first two
	initState.erase(initState.begin()+2, initState.end());

	BOOST_CHECK(stateDiffBelowTol(finalState, initState, 1e-12));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_DIST){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);

	// Constraint_tp::DIST
	double matchDistConData[] = {1, 0.2};
	Constraint matchDistCon(Constraint_tp::DIST, 3, matchDistConData, 2);
	halfLyapNodeset.addConstraint(matchDistCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(matchDistCon.getID());
	double dist = sqrt(pow(finalState[0] - 1 + sys.getMu(),2) + pow(finalState[1], 2));
	BOOST_CHECK_SMALL(dist - matchDistConData[1], 1e-12);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_DIST_EQUAL_ARC){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);

	// Constraint_tp::DIST
	double matchDistConData[] = {1, 0.2};
	Constraint matchDistCon(Constraint_tp::DIST, 3, matchDistConData, 2);
	halfLyapNodeset.addConstraint(matchDistCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(matchDistCon.getID());
	double dist = sqrt(pow(finalState[0] - 1 + sys.getMu(),2) + pow(finalState[1], 2));
	BOOST_CHECK_SMALL(dist - matchDistConData[1], 1e-12);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_MIN_DIST){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);

	// Constraint_tp::MIN_DIST
	double matchDistConData[] = {1, 0.4};
	Constraint matchDistCon(Constraint_tp::MIN_DIST, 3, matchDistConData, 2);
	halfLyapNodeset.addConstraint(matchDistCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(matchDistCon.getID());
	double dist = sqrt(pow(finalState[0] - 1 + sys.getMu(), 2) + pow(finalState[1], 2));

	BOOST_CHECK_GE(dist, matchDistConData[1]);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_MIN_DIST_EQUAL_ARC){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);

	// Constraint_tp::MIN_DIST
	double matchDistConData[] = {1, 0.4};
	Constraint matchDistCon(Constraint_tp::MIN_DIST, 3, matchDistConData, 2);
	halfLyapNodeset.addConstraint(matchDistCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(matchDistCon.getID());
	double dist = sqrt(pow(finalState[0] - 1 + sys.getMu(), 2) + pow(finalState[1], 2));

	BOOST_CHECK_GE(dist, matchDistConData[1]);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_MAX_DIST){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);

	// Constraint_tp::MIN_DIST
	double matchDistConData[] = {1, 0.15};
	Constraint matchDistCon(Constraint_tp::MAX_DIST, 3, matchDistConData, 2);
	halfLyapNodeset.addConstraint(matchDistCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(matchDistCon.getID());
	double dist = sqrt(pow(finalState[0] - 1 + sys.getMu(), 2) + pow(finalState[1], 2));

	BOOST_CHECK_LE(dist, matchDistConData[1]);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_MAX_DIST_EQUAL_ARC){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);

	// Constraint_tp::MIN_DIST
	double matchDistConData[] = {1, 0.15};
	Constraint matchDistCon(Constraint_tp::MAX_DIST, 3, matchDistConData, 2);
	halfLyapNodeset.addConstraint(matchDistCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(matchDistCon.getID());
	double dist = sqrt(pow(finalState[0] - 1 + sys.getMu(), 2) + pow(finalState[1], 2));

	BOOST_CHECK_LE(dist, matchDistConData[1]);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_MAX_DELTA_V){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);

	// Constraint_tp::MAX_DELTA_V
	std::vector<int> dvNodes {3};
	std::vector<double> state = halfLyapNodeset.getState(dvNodes[0]);
	state[3] += 0.01;
	state[4] += 0.1;
	halfLyapNodeset.setState(dvNodes[0], state);	// Perturb the velocity of this state to create a discontinuity
	halfLyapNodeset.allowDV_at(dvNodes);	// Allow the perturbed node to have a delta-v
	double maxDVConData = 0.01;
	Constraint dVCon(Constraint_tp::MAX_DELTA_V, 0, &maxDVConData, 1);
	halfLyapNodeset.addConstraint(dVCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	MultShootData itData(&correctedSet);
	BOOST_CHECK_NO_THROW(itData = corrector.multShoot(&halfLyapNodeset, &correctedSet));

	double totalDV = MultShootEngine::getTotalDV(&itData);
	BOOST_CHECK_LE(totalDV, maxDVConData);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_MAX_DELTA_V_EQUAL_ARC){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);

	// Constraint_tp::MAX_DELTA_V
	std::vector<int> dvNodes {3};
	std::vector<double> state = halfLyapNodeset.getState(dvNodes[0]);
	state[3] += 0.01;
	state[4] += 0.1;
	halfLyapNodeset.setState(dvNodes[0], state);	// Perturb the velocity of this state to create a discontinuity
	halfLyapNodeset.allowDV_at(dvNodes);	// Allow the perturbed node to have a delta-v
	double maxDVConData = 0.01;
	Constraint dVCon(Constraint_tp::MAX_DELTA_V, 0, &maxDVConData, 1);
	halfLyapNodeset.addConstraint(dVCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	MultShootData itData(&correctedSet);
	BOOST_CHECK_NO_THROW(itData = corrector.multShoot(&halfLyapNodeset, &correctedSet));

	double totalDV = MultShootEngine::getTotalDV(&itData);
	BOOST_CHECK_LE(totalDV, maxDVConData);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_DELTA_V){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);

	// Constraint_tp::DELTA_V
	std::vector<int> dvNodes {3};
	std::vector<double> state = halfLyapNodeset.getState(dvNodes[0]);
	state[3] += 0.01;
	state[4] += 0.1;
	halfLyapNodeset.setState(dvNodes[0], state);	// Perturb the velocity of this state to create a discontinuity
	halfLyapNodeset.allowDV_at(dvNodes);	// Allow the perturbed node to have a delta-v
	double maxDVConData = 0.01;
	Constraint dVCon(Constraint_tp::DELTA_V, 0, &maxDVConData, 1);
	halfLyapNodeset.addConstraint(dVCon);
	
	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	MultShootData itData(&correctedSet);
	BOOST_CHECK_NO_THROW(itData = corrector.multShoot(&halfLyapNodeset, &correctedSet));

	double totalDV = MultShootEngine::getTotalDV(&itData);
	BOOST_CHECK_SMALL(totalDV - maxDVConData, 1e-10);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_DELTA_V_EQUAL_ARC){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);

	// Constraint_tp::DELTA_V
	std::vector<int> dvNodes {3};
	std::vector<double> state = halfLyapNodeset.getState(dvNodes[0]);
	state[3] += 0.01;
	state[4] += 0.1;
	halfLyapNodeset.setState(dvNodes[0], state);	// Perturb the velocity of this state to create a discontinuity
	halfLyapNodeset.allowDV_at(dvNodes);	// Allow the perturbed node to have a delta-v
	double maxDVConData = 0.01;
	Constraint dVCon(Constraint_tp::DELTA_V, 0, &maxDVConData, 1);
	halfLyapNodeset.addConstraint(dVCon);
	
	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	MultShootData itData(&correctedSet);
	BOOST_CHECK_NO_THROW(itData = corrector.multShoot(&halfLyapNodeset, &correctedSet));

	double totalDV = MultShootEngine::getTotalDV(&itData);
	BOOST_CHECK_SMALL(totalDV - maxDVConData, 1e-10);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_JACOBI){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);

	double jacobiData = 3.1149;
	Constraint jacobiCon(Constraint_tp::JC, 0, &jacobiData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(jacobiCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));
	BOOST_CHECK_SMALL(correctedSet.getJacobiByIx(0) - jacobiData, 1e-12);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_JACOBI_EQUAL_ARC){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);

	double jacobiData = 3.1149;
	Constraint jacobiCon(Constraint_tp::JC, 0, &jacobiData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(jacobiCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));
	BOOST_CHECK_SMALL(correctedSet.getJacobiByIx(0) - jacobiData, 1e-12);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_TOF){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);

	double tofData = 2.5;
	Constraint tofCon(Constraint_tp::TOF, 0, &tofData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(tofCon);
	
	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));
	BOOST_CHECK_SMALL(correctedSet.getTotalTOF() - tofData, 1e-12);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_TOF_EQUAL_ARC){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);

	double tofData = 2.5;
	Constraint tofCon(Constraint_tp::TOF, 0, &tofData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(tofCon);
	
	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));
	BOOST_CHECK_SMALL(correctedSet.getTotalTOF() - tofData, 1e-12);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_APSE){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);

	double apseData = 1;
	Constraint apseCon(Constraint_tp::APSE, 4, &apseData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(apseCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(apseCon.getID());
	const DynamicsModel *model = sys.getDynamicsModel();
	std::vector<double> primPos = model->getPrimPos(0, &sys);
	double dx = finalState[0] - primPos[apseData*3 + 0];
	double dy = finalState[1] - primPos[apseData*3 + 1];
	double dz = finalState[2] - primPos[apseData*3 + 2];
	double rdot = dx*finalState[3] + dy*finalState[4] + dz*finalState[5];

	BOOST_CHECK_SMALL(rdot, 1e-12);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_APSE_EQUAL_ARC){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setRevTime(em_lyap_T < 0);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T, 6, &halfLyapNodeset);

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);

	// APSE
	double apseData = 1;
	Constraint apseCon(Constraint_tp::APSE, 4, &apseData, 1);
	halfLyapNodeset.clearAllConstraints();
	halfLyapNodeset.addConstraint(apseCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(apseCon.getID());
	const DynamicsModel *model = sys.getDynamicsModel();
	std::vector<double> primPos = model->getPrimPos(0, &sys);
	double dx = finalState[0] - primPos[apseData*3 + 0];
	double dy = finalState[1] - primPos[apseData*3 + 1];
	double dz = finalState[2] - primPos[apseData*3 + 2];
	double rdot = dx*finalState[3] + dy*finalState[4] + dz*finalState[5];

	BOOST_CHECK_SMALL(rdot, 1e-12);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_SEG_CONT_PV){
	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp forwardArc(&sys), reverseArc(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T/2, 4, &forwardArc);
	sim.setRevTime(true);
	sim.runSim_manyNodes(em_lyap_ic, em_lyap_T/2, 4, &reverseArc);
	Arcset_cr3bp doubleSrcLyap = forwardArc;
	doubleSrcLyap.appendSetAtNode(&reverseArc, 0, 0, 0);
	// doubleSrcLyap.print();
	// doubleSrcLyap.printInChrono();

	MultShootEngine corrector;
	// corrector.setVerbosity(Verbosity_tp::DEBUG);
	corrector.setEqualArcTime(false);
	

	double contData[] = {4, 4, NAN, 4, 4, NAN};
	Constraint contCon(Constraint_tp::SEG_CONT_PV, 2, contData, 6);
	doubleSrcLyap.addConstraint(contCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&doubleSrcLyap, corrector, Verbosity_tp::NO_MSG));
	MultShootData it(&correctedSet);
	BOOST_CHECK_NO_THROW(it = corrector.multShoot(&doubleSrcLyap, &correctedSet));

	Arcset forwardTraj = it.propSegs[correctedSet.getSegIx(contCon.getID())];
	Arcset reverseTraj = it.propSegs[correctedSet.getSegIx(contData[0])];
	std::vector<double> for_lastState = forwardTraj.getStateByIx(-1);
	std::vector<double> rev_lastState = reverseTraj.getStateByIx(-1);
	double sum = 0;
	for(int i = 0; i < 6; i++){
		if(!std::isnan(contData[i]))
			sum += pow(for_lastState[i] - rev_lastState[i], 2);
	}
	double dist = std::sqrt(sum);

	BOOST_CHECK_SMALL(dist, 1e-12);
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_DOUBLE_SOURCE){
	SysData_cr3bp sys("earth", "moon");

	// EM L2 Butterfly Orbit
	std::vector<double> ic {1.0639767173456007, 0, 0.1644973017995331, 0, -0.0311246472806882, 0};
	double tof = 3.1256890778;

	Arcset_cr3bp halfPlus(&sys), halfMinus(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(ic, tof, 4, &halfPlus);
	sim.setRevTime(true);
	sim.runSim_manyNodes(ic, tof, 4, &halfMinus);
	halfPlus.appendSetAtNode(&halfMinus, 0, 0, 0);

	MultShootEngine corrector;

	// Constraint_tp::STATE
	double stateConData[] = {1.1, 0, NAN, 0, NAN, 0};
	Constraint stateCon(Constraint_tp::STATE, 3, stateConData, 6);
	halfPlus.addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfPlus, corrector, Verbosity_tp::NO_MSG));

	Arcset_cr3bp correctedSet(&sys);
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfPlus, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_AUTO_TEST_SUITE_END()
#define BOOST_TEST_MODULE MultipleShootingConstraints_CR3BP_LT

#include <boost/test/unit_test.hpp>

#include <vector>

#include "Arcset_cr3bp_lt.hpp"
#include "Calculations.hpp"
#include "Constraint.hpp"
#include "MultShootEngine.hpp"
#include "ControlLaw_cr3bp_lt.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "SimEngine.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Utilities.hpp"

using namespace astrohelion;

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
//* CR3BP Earth-Moon Constraints
//************************************************************

BOOST_AUTO_TEST_SUITE(CR3BP_LT_EarthMoon)

BOOST_AUTO_TEST_CASE(CR3BP_LT_EM_CONT_Law1){
	SysData_cr3bp_lt sys("earth", "moon", 14);
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::Law_tp::CONST_C_2D_LEFT, 12e-3, 1500);
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period
	
	Arcset_cr3bp_lt halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(ic, T, 2, &halfLyapNodeset, &law);

	MultShootEngine corrector;
	// corrector.setVerbosity(Verbosity_tp::SOME_MSG);
	corrector.setEqualArcTime(false);
	// corrector.setVerbosity(Verbosity_tp::ALL_MSG);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG, true));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_EM_CONT_Law2){
	SysData_cr3bp_lt sys("earth", "moon", 14);
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::Law_tp::CONST_C_2D_RIGHT, 12e-3, 1500);
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period
	
	Arcset_cr3bp_lt halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(ic, T, 2, &halfLyapNodeset, &law);

	MultShootEngine corrector;
	// corrector.setVerbosity(Verbosity_tp::SOME_MSG);
	corrector.setEqualArcTime(false);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));
}//====================================================

// BOOST_AUTO_TEST_CASE(CR3BP_LT_EM_CONT_Law3){
// 	SysData_cr3bp_lt sys("earth", "moon", 12e-3, 1500, 14);
// 	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
// 	double T = 3.02796323553149;	// EM L1 Period
// 	Arcset_cr3bp_lt halfLyapNodeset(&sys, ic, T, 6, Arcset::TIME, ControlLaw_cr3bp_lt::Law_tp::PRO_VEL);	// Create a nodeset
// 	Arcset_cr3bp_lt correctedSet(&sys);

// 	MultShootEngine corrector;
// 	corrector.setEqualArcTime(false);

// 	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
// 	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));
// }//====================================================

// BOOST_AUTO_TEST_CASE(CR3BP_LT_EM_CONT_Law4){
// 	SysData_cr3bp_lt sys("earth", "moon", 12e-3, 1500, 14);
// 	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
// 	double T = 3.02796323553149;	// EM L1 Period
// 	Arcset_cr3bp_lt halfLyapNodeset(&sys, ic, T, 6, Arcset::TIME, ControlLaw_cr3bp_lt::Law_tp::ANTI_VEL);	// Create a nodeset
// 	Arcset_cr3bp_lt correctedSet(&sys);

// 	MultShootEngine corrector;
// 	corrector.setEqualArcTime(false);

// 	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
// 	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));
// }//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_EM_STATE){
	SysData_cr3bp_lt sys("earth", "moon", 14);
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::Law_tp::CONST_C_2D_LEFT, 12e-3, 1500);
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period

	Arcset_cr3bp_lt halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(ic, T, 3, &halfLyapNodeset, &law);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);
	// corrector.setVerbosity(Verbosity_tp::ALL_MSG);

	double stateConData[] = {-0.1, -0.3, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 2, stateConData, 7);
	halfLyapNodeset.addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));

	// correctedSet.saveToMat("data/lt_correctedSet.mat");
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_EM_STATE_EQUAL_ARC){
	SysData_cr3bp_lt sys("earth", "moon", 14);
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::Law_tp::CONST_C_2D_LEFT, 12e-3, 1500);
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period

	Arcset_cr3bp_lt halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(ic, T, 3, &halfLyapNodeset, &law);

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);
	// corrector.setVerbosity(Verbosity_tp::ALL_MSG);

	double stateConData[] = {-0.1, -0.4, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 2, stateConData, 7);
	halfLyapNodeset.addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));

	// correctedSet.saveToMat("data/lt_correctedSet.mat");
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_EM_STATE_ENDSEG){
	SysData_cr3bp_lt sys("earth", "moon", 14);
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::Law_tp::CONST_C_2D_LEFT, 12e-3, 1500);
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period

	Arcset_cr3bp_lt halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(ic, T, 3, &halfLyapNodeset, &law);
	halfLyapNodeset.deleteNode(2);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);
	// corrector.setVerbosity(Verbosity_tp::ALL_MSG);

	double stateConData[] = {-0.1, -0.4, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::ENDSEG_STATE, 1, stateConData, 7);
	halfLyapNodeset.addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> fullState = correctedSet.getSegByIx(stateCon.getID()).getStateByRow(-1, 56);
	std::vector<double> finalState(fullState.begin(), fullState.begin()+7);
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_EM_STATE_ENDSEG_EQUAL_ARC){
	SysData_cr3bp_lt sys("earth", "moon", 14);
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::Law_tp::CONST_C_2D_LEFT, 12e-3, 1500);
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period

	Arcset_cr3bp_lt halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(ic, T, 3, &halfLyapNodeset, &law);
	halfLyapNodeset.deleteNode(2);

	MultShootEngine corrector;
	corrector.setEqualArcTime(true);
	// corrector.setVerbosity(Verbosity_tp::ALL_MSG);

	double stateConData[] = {-0.1, -0.4, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::ENDSEG_STATE, 1, stateConData, 7);
	halfLyapNodeset.addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> fullState = correctedSet.getSegByIx(stateCon.getID()).getStateByRow(-1, 56);
	std::vector<double> finalState(fullState.begin(), fullState.begin()+7);
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_EM_STATE_RMSTATE){
	SysData_cr3bp_lt sys("earth", "moon", 14);
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::Law_tp::CONST_C_2D_LEFT, 12e-3, 1500);
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period

	Arcset_cr3bp_lt halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(ic, T/2.8, 3, &halfLyapNodeset, &law);

	MultShootEngine corrector;
	corrector.setEqualArcTime(false);
	// corrector.setVerbosity(Verbosity_tp::ALL_MSG);

	double stateConData[] = {NAN, 0, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 2, stateConData, 7);
	halfLyapNodeset.addConstraint(stateCon);

	Constraint rmState(Constraint_tp::RM_STATE, 0, nullptr, 0);
	halfLyapNodeset.addConstraint(rmState);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));

	// correctedSet.saveToMat("data/lt_correctedSet.mat");
}//====================================================

BOOST_AUTO_TEST_SUITE_END()
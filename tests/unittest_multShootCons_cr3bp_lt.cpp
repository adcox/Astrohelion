/**
 * Useful info about Boost Unit Tests:
 * <http://www.boost.org/doc/libs/1_62_0/libs/test/doc/html/index.html>
 */
#define BOOST_TEST_MODULE MultipleShootingConstraints_CR3BP_LT

#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <vector>

#include "Arcset_cr3bp_lt.hpp"
#include "Calculations.hpp"
#include "Constraint.hpp"
#include "MultShootEngine.hpp"
#include "ControlLaw_cr3bp_lt.hpp"
#include "Exceptions.hpp"
#include "SimEngine.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Utilities.hpp"

using namespace astrohelion;
using namespace boost::unit_test;
using ltlaw = astrohelion::ControlLaw_cr3bp_lt;

std::vector<double> emL1Lyap_ic {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
std::vector<double> emL3Lyap_ic {-0.628097117249632, 0, 0, 0, -0.87229410151913, 0, 1};	// EM L3
double emL1Lyap_T = 3.02797;	// EM L1 Period
double emL3Lyap_T = 6.2238;

// All combinations of the control laws that require 0 control states
std::vector<unsigned int> law_bases = { 
	ltlaw::VEL_PT | ltlaw::VEL_PRO, ltlaw::VEL_PT | ltlaw::VEL_ANTI, 
	ltlaw::CONST_C_2D | ltlaw::VEL_LEFT, ltlaw::CONST_C_2D | ltlaw::VEL_RIGHT};

// All variable thrust parameterizations
std::vector<unsigned int> law_varF = {ltlaw::VAR_F_BND, ltlaw::VAR_F_UBND};

// Define the basic control law for every single case
const unsigned int basicLaw = ltlaw::CSI_VAR_M;

// All the different ways to parameterize time in the multiple shooting algorithm
std::vector<MSTOF_tp> tofTypes {MSTOF_tp::VAR_FREE, MSTOF_tp::VAR_FIXSIGN,
	MSTOF_tp::VAR_EQUALARC};

bool stateDiffBelowTol(std::vector<double>, double*, double);
bool stateDiffBelowTol(std::vector<double>, std::vector<double>, double);

/**
 *  @brief Determine if the difference bewteen two state vectors is less than 
 *  the desired tolerance
 * 
 *  @param data vector of state values from the corrections process
 *  @param correct array of state values from the constraint
 *  @param tol desired numerical tolerance
 *  @return whether or not `data` and `correct` are equal
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

bool stateDiffBelowTol(std::vector<double> data, std::vector<double> correct, 
	double tol){

	return stateDiffBelowTol(data, &(correct[0]), tol);
}//====================================================

//************************************************************
//* CR3BP Earth-Moon Constraints, Constant Thrust Laws
//************************************************************
BOOST_AUTO_TEST_SUITE(CONST_F)

BOOST_DATA_TEST_CASE(test_continuity, 
	data::make(law_bases) * data::make(tofTypes), baseType, tofType){

	// Use the specified base and set thrust to be constant
	unsigned int lawType = basicLaw | baseType | ltlaw::CONST_F;

	SysData_cr3bp_lt sys("earth", "moon", 14);
	std::vector<double> ltParams {0.3, 1500};
	ControlLaw_cr3bp_lt law(lawType, ltParams);

	Arcset_cr3bp_lt halfLyapArcset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL3Lyap_ic, emL3Lyap_T, 2, &halfLyapArcset, &law);

	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::NO_MSG);
	corrector.setTOFType(tofType);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapArcset, 
		corrector, Verbosity_tp::NO_MSG, false));

	try{
		corrector.multShoot(&halfLyapArcset, &correctedSet);
	}catch(const Exception &e){
		printErr("Error (%s, %s): %s\n", 
			ControlLaw_cr3bp_lt::typeToString(lawType).c_str(),
			MSTOF_tp_cStr(tofType), e.what());
		char filename[128];
		snprintf(filename, 128, "%s_%s_input.mat", 
			ControlLaw_cr3bp_lt::typeToString(lawType).c_str(),
			MSTOF_tp_cStr(tofType));
		halfLyapArcset.saveToMat(filename);

		snprintf(filename, 128, "%s_%s_corrected.mat", 
			ControlLaw_cr3bp_lt::typeToString(lawType).c_str(),
			MSTOF_tp_cStr(tofType));
		correctedSet.saveToMat(filename);
	}
}//====================================================

BOOST_DATA_TEST_CASE(test_stateConstraint, data::make(law_bases) * 
	data::make(tofTypes), baseType, tofType){

	// Use the specified base and set thrust to be constant
	unsigned int lawType = basicLaw | baseType | ltlaw::CONST_F;

	SysData_cr3bp_lt sys("earth", "moon", 500);
	std::vector<double> ltParams {9e-3, 1500};
	ControlLaw_cr3bp_lt law(lawType, ltParams);

	Arcset_cr3bp_lt halfLyapArcset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, emL1Lyap_T/2, 3, &halfLyapArcset, &law);

	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::NO_MSG);
	corrector.setTOFType(tofType);
	// corrector.setVerbosity(Verbosity_tp::NO_MSG);

	double stateConData[] = {0.78, 0, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 2, stateConData, 7);
	halfLyapArcset.addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapArcset, 
		corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapArcset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_DATA_TEST_CASE(test_endSegCon, data::make(law_bases) * 
	data::make(tofTypes), baseType, tofType){

	// Use the specified base and set thrust to be constant
	unsigned int lawType = basicLaw | baseType | ltlaw::CONST_F;

	SysData_cr3bp_lt sys("earth", "moon", 500);
	std::vector<double> ltParams {9e-3, 1500};
	ControlLaw_cr3bp_lt law(lawType, ltParams);


	Arcset_cr3bp_lt halfLyapArcset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, emL1Lyap_T/2, 3, &halfLyapArcset, &law);
	halfLyapArcset.deleteNode(2);

	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::NO_MSG);
	corrector.setTOFType(tofType);
	// corrector.setVerbosity(Verbosity_tp::NO_MSG);

	double stateConData[] = {0.78, 0, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::ENDSEG_STATE, 1, stateConData, 7);
	halfLyapArcset.addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapArcset, 
		corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapArcset, &correctedSet));

	std::vector<double> fullState = correctedSet.getSegByIx(
		stateCon.getID()).getStateByRow(-1);
	std::vector<double> finalState(fullState.begin(), fullState.begin()+7);
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_DATA_TEST_CASE(test_rmState, data::make(law_bases) * 
	data::make(tofTypes), baseType, tofType){

	// Use the specified base and set thrust to be constant
	unsigned int lawType = basicLaw | baseType | ltlaw::CONST_F;

	SysData_cr3bp_lt sys("earth", "moon", 14);
	std::vector<double> ltParams {1e-2, 1500};
	ControlLaw_cr3bp_lt law(lawType, ltParams);

	Arcset_cr3bp_lt halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, emL1Lyap_T/2.8, 3, &halfLyapNodeset, &law);

	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::NO_MSG);
	corrector.setTOFType(tofType);
	corrector.setTol(1e-11);
	// corrector.setVerbosity(Verbosity_tp::NO_MSG);

	double stateConData[] = {NAN, 0, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 2, stateConData, 7);
	halfLyapNodeset.addConstraint(stateCon);

	Constraint rmState(Constraint_tp::RM_STATE, 0, nullptr, 0);
	halfLyapNodeset.addConstraint(rmState);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, 
		corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));

	// correctedSet.saveToMat("data/lt_correctedSet.mat");
}//====================================================

BOOST_DATA_TEST_CASE(test_HLTConstraint, data::make(law_bases), baseType){

	// Use the specified base and set thrust to be constant
	unsigned int lawType = basicLaw | baseType | ltlaw::CONST_F;

	SysData_cr3bp_lt sys("earth", "moon", 500);
	std::vector<double> ltParams {9e-3, 1500};
	ControlLaw_cr3bp_lt law(lawType, ltParams);

	Arcset_cr3bp_lt halfLyapArcset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, emL1Lyap_T/2, 3, &halfLyapArcset, &law);

	MultShootEngine corrector;
	// corrector.setVerbosity(Verbosity_tp::SOME_MSG);
	corrector.setVerbosity(Verbosity_tp::NO_MSG);

	double hltConData[] = {-1.6};
	Constraint hltCon(Constraint_tp::HLT, 2, hltConData, 1);
	halfLyapArcset.addConstraint(hltCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapArcset,
		corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapArcset, &correctedSet));

	std::vector<double> q = correctedSet.getState(hltCon.getID());
	double H = DynamicsModel_cr3bp_lt::getHamiltonian(0, &(q[0]), &sys, &law);
	BOOST_CHECK_SMALL(H - hltConData[0], 1e-12);
}//====================================================

BOOST_AUTO_TEST_SUITE_END()
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(CR3BP_LT_GENERAL_DIR_CONST_F)

auto randAlphaGen = data::random(-PI, PI);
auto randBetaGen = data::random(-PI/2, PI/2);

BOOST_DATA_TEST_CASE(test_continuity, 
	randAlphaGen ^ randBetaGen ^ data::xrange(5) * data::make(tofTypes), 
	alpha, beta, index, tofType){

	(void) index;
	unsigned int lawType = basicLaw | ltlaw::GENERAL | ltlaw::CONST_F;

	SysData_cr3bp_lt sys("earth", "moon", 14);
	std::vector<double> ltParams {1e-2, 1500};
	ControlLaw_cr3bp_lt law(lawType, ltParams);
	std::vector<double> thrustAngles {alpha, beta};
	
	Arcset_cr3bp_lt halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, thrustAngles, 0, emL1Lyap_T, 2, 
		&halfLyapNodeset, &law);

	// Add control continuity constraint for full-rank Jacobian 
	// (check ALL the available partials)
	std::vector<double> conData(law.getNumStates(), 1);
	Constraint ctrlCon(Constraint_tp::CONT_CTRL, 0, conData);
	halfLyapNodeset.addConstraint(ctrlCon);


	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::NO_MSG);
	corrector.setTOFType(tofType);
	corrector.setTol(1e-11);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, 
		corrector, Verbosity_tp::NO_MSG));
	// BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	try{
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
	}catch(const Exception &e){
		printErr("Error: %s\n", e.what());
	}

	for(unsigned int n = 0; n < correctedSet.getNumNodes(); n++){
		BOOST_CHECK_NO_THROW(correctedSet.getNodeByIx(n).getExtraParamVec(PARAMKEY_CTRL));
		std::vector<double> ctrlState = 
			correctedSet.getNodeByIx(n).getExtraParamVec(PARAMKEY_CTRL);
		BOOST_CHECK(ctrlState.size() == law.getNumStates());
	}
}//====================================================

BOOST_DATA_TEST_CASE(test_stateConstraint, 
	randAlphaGen ^ randBetaGen ^ data::xrange(5) * data::make(tofTypes), 
	alpha, beta, index, tofType){

	(void) index;
	unsigned int lawType = basicLaw | ltlaw::GENERAL | ltlaw::CONST_F;

	SysData_cr3bp_lt sys("earth", "moon", 500);
	std::vector<double> ltParams {1e-2, 1500};
	ControlLaw_cr3bp_lt law(lawType, ltParams);
	std::vector<double> thrustAngles {alpha, beta};

	Arcset_cr3bp_lt halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, thrustAngles, 0, emL1Lyap_T/2, 2, 
		&halfLyapNodeset, &law);

	// Add control continuity constraint for full-rank Jacobian 
	// (check ALL the available partials)
	std::vector<double> conData(law.getNumStates(), 1);
	Constraint ctrlCon(Constraint_tp::CONT_CTRL, 0, conData);

	std::vector<double> stateConData {0.78, 0, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 1, stateConData);

	halfLyapNodeset.addConstraint(ctrlCon);
	halfLyapNodeset.addConstraint(stateCon);

	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::NO_MSG);
	// corrector.setVerbosity(Verbosity_tp::NO_MSG);
	corrector.setTOFType(tofType);
	corrector.setTol(1e-11);
	
	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, 
		corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));
}//====================================================

BOOST_DATA_TEST_CASE(test_HLTConstraint, 
	randAlphaGen ^ randBetaGen ^ data::xrange(5), 
	alpha, beta, index){

	(void) index;
	unsigned int lawType = basicLaw | ltlaw::GENERAL | ltlaw::CONST_F;

	SysData_cr3bp_lt sys("earth", "moon", 500);
	std::vector<double> ltParams {9e-3, 1500};
	ControlLaw_cr3bp_lt law(lawType, ltParams);
	std::vector<double> thrustAngles {alpha, beta};

	Arcset_cr3bp_lt halfLyapArcset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, thrustAngles, 0, emL1Lyap_T/2, 3, 
		&halfLyapArcset, &law);

	MultShootEngine corrector;
	// corrector.setVerbosity(Verbosity_tp::SOME_MSG);
	corrector.setVerbosity(Verbosity_tp::NO_MSG);

	double hltConData[] = {-1.6};
	Constraint hltCon(Constraint_tp::HLT, 2, hltConData, 1);
	halfLyapArcset.addConstraint(hltCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapArcset,
		corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapArcset, &correctedSet));

	std::vector<double> q = correctedSet.getState(hltCon.getID());
	std::vector<double> g = 
		correctedSet.getNode(hltCon.getID()).getExtraParamVec(PARAMKEY_CTRL);
	q.insert(q.end(), g.begin(), g.end());
	double H = DynamicsModel_cr3bp_lt::getHamiltonian(0, &(q[0]), &sys, &law);
	BOOST_CHECK_SMALL(H - hltConData[0], 1e-12);
}//====================================================
BOOST_AUTO_TEST_SUITE_END()
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(CR3BP_LT_GEN_INERT_CONST_F)

auto randPsiGen = data::random(-PI, PI);
auto randBetaGen = data::random(-PI/2, PI/2);

BOOST_DATA_TEST_CASE(test_continuity, 
	randPsiGen ^ randBetaGen ^ data::xrange(5) * data::make(tofTypes), 
	psi, beta, index, tofType){

	(void) index;
	unsigned int lawType = basicLaw | ltlaw::GEN_INERT | ltlaw::CONST_F;

	SysData_cr3bp_lt sys("earth", "moon", 14);
	std::vector<double> ltParams {1.234, 1e-2, 1500};
	ControlLaw_cr3bp_lt law(lawType, ltParams);
	std::vector<double> thrustAngles {psi, beta};
	
	Arcset_cr3bp_lt halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, thrustAngles, 0, 0.25*emL1Lyap_T, 2, 
		&halfLyapNodeset, &law);

	// Add control continuity constraint for full-rank Jacobian 
	// (check ALL the available partials)
	std::vector<double> conData(law.getNumStates(), 1);
	Constraint ctrlCon(Constraint_tp::CONT_CTRL, 0, conData);
	halfLyapNodeset.addConstraint(ctrlCon);


	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::NO_MSG);
	// corrector.setVerbosity(Verbosity_tp::NO_MSG);
	corrector.setTOFType(tofType);
	corrector.setTol(1e-11);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, 
		corrector, Verbosity_tp::NO_MSG, false));
	// BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));

	try{
		corrector.multShoot(&halfLyapNodeset, &correctedSet);
	}catch(const Exception &e){
		printErr("Error: %s\n", e.what());
	}

	for(unsigned int n = 0; n < correctedSet.getNumNodes(); n++){
		BOOST_CHECK_NO_THROW(correctedSet.getNodeByIx(n).getExtraParamVec(PARAMKEY_CTRL));
		std::vector<double> ctrlState = 
			correctedSet.getNodeByIx(n).getExtraParamVec(PARAMKEY_CTRL);
		BOOST_CHECK(ctrlState.size() == law.getNumStates());
	}
}//====================================================

BOOST_DATA_TEST_CASE(test_stateConstraint, 
	randPsiGen ^ randBetaGen ^ data::xrange(5) * data::make(tofTypes), 
	psi, beta, index, tofType){

	(void) index;
	unsigned int lawType = basicLaw | ltlaw::GEN_INERT | ltlaw::CONST_F;

	SysData_cr3bp_lt sys("earth", "moon", 500);
	std::vector<double> ltParams {1.234, 1e-2, 1500};
	ControlLaw_cr3bp_lt law(lawType, ltParams);
	std::vector<double> thrustAngles {psi, beta};

	Arcset_cr3bp_lt halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, thrustAngles, 0, emL1Lyap_T/2, 2, 
		&halfLyapNodeset, &law);

	// Add control continuity constraint for full-rank Jacobian 
	// (check ALL the available partials)
	std::vector<double> conData(law.getNumStates(), 1);
	Constraint ctrlCon(Constraint_tp::CONT_CTRL, 0, conData);

	std::vector<double> stateConData {0.78, 0, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 1, stateConData);

	halfLyapNodeset.addConstraint(ctrlCon);
	halfLyapNodeset.addConstraint(stateCon);

	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::NO_MSG);
	// corrector.setVerbosity(Verbosity_tp::NO_MSG);
	corrector.setTOFType(tofType);
	corrector.setTol(1e-11);
	
	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, 
		corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapNodeset, &correctedSet));
}//====================================================

BOOST_DATA_TEST_CASE(test_HLTConstraint, 
	randPsiGen ^ randBetaGen ^ data::xrange(5), 
	psi, beta, index){

	(void) index;
	unsigned int lawType = basicLaw | ltlaw::GEN_INERT | ltlaw::CONST_F;

	SysData_cr3bp_lt sys("earth", "moon", 500);
	std::vector<double> ltParams {1.234, 9e-3, 1500};
	ControlLaw_cr3bp_lt law(lawType, ltParams);
	std::vector<double> thrustAngles {psi, beta};

	Arcset_cr3bp_lt halfLyapArcset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, thrustAngles, 0, emL1Lyap_T/2, 3, 
		&halfLyapArcset, &law);

	MultShootEngine corrector;
	// corrector.setVerbosity(Verbosity_tp::SOME_MSG);
	corrector.setVerbosity(Verbosity_tp::NO_MSG);

	double hltConData[] = {-1.6};
	Constraint hltCon(Constraint_tp::HLT, 2, hltConData, 1);
	halfLyapArcset.addConstraint(hltCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapArcset,
		corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halfLyapArcset, &correctedSet));

	std::vector<double> q = correctedSet.getState(hltCon.getID());
	std::vector<double> g = 
		correctedSet.getNode(hltCon.getID()).getExtraParamVec(PARAMKEY_CTRL);
	q.insert(q.end(), g.begin(), g.end());
	double H = DynamicsModel_cr3bp_lt::getHamiltonian(0, &(q[0]), &sys, &law);
	BOOST_CHECK_SMALL(H - hltConData[0], 1e-12);
}//====================================================
BOOST_AUTO_TEST_SUITE_END()
//-----------------------------------------------------------------------------

//************************************************************
//* CR3BP Earth-Moon Constraints, Variable Thrust Laws
//************************************************************
BOOST_AUTO_TEST_SUITE(CR3BP_LT_VAR_F)
// The vector varF_lawTypes does not include general pointing,
// so params = {fmax, Isp} and ctrl = {g} where f = fmax/2*(sin(g) + 1)

BOOST_DATA_TEST_CASE(test_continuity, 
	data::make(law_bases) * data::make(law_varF) * data::make(tofTypes), 
	baseType, fType, tofType){

	// Construct the ID
	unsigned int lawType = basicLaw | baseType | fType;

	// Set the parameters and states
	std::vector<double> ltParams {1e-1, 3000};	// {fmax, Isp}
	std::vector<double> ctrlState {};
	switch(fType){
		case ltlaw::VAR_F_BND:
			ctrlState.push_back(-0.9273);	// yields f = 1e-2
			break;
		case ltlaw::VAR_F_UBND:
			ctrlState.push_back(0.3162);	// yields f = 1e-2
			break;
		default: break;
	}
	SysData_cr3bp_lt sys("earth", "moon", 14);

	ControlLaw_cr3bp_lt law(lawType, ltParams);

	Arcset_cr3bp_lt halfLyapArcset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL3Lyap_ic, ctrlState, 0, emL3Lyap_T, 2, 
		&halfLyapArcset, &law);

	// Remove final ctrl state; no constraints on final ctrl state, 
	// results in column of zeros
	Constraint rmEndCtrl(Constraint_tp::RM_CTRL, 
		halfLyapArcset.getNodeByIx(-1).getID(), nullptr, 0);
	halfLyapArcset.addConstraint(rmEndCtrl);

	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::NO_MSG);
	corrector.setTOFType(tofType);
	corrector.setTol(1e-11);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapArcset, 
		corrector, Verbosity_tp::NO_MSG, false));

	// waitForUser();

	try{
		corrector.multShoot(&halfLyapArcset, &correctedSet);
	}catch(const Exception &e){
		printErr("Error (%s, %s): %s\n", 
			ControlLaw_cr3bp_lt::typeToString(lawType).c_str(),
			MSTOF_tp_cStr(tofType), e.what());
		// char filename[128];
		// sprintf(filename, "%s_%s_input.mat", 
		// 	ControlLaw_cr3bp_lt::typeToString(lawType).c_str(),
		// 	MSTOF_tp_cStr(tofType));
		// halfLyapArcset.saveToMat(filename);

		// sprintf(filename, "%s_%s_corrected.mat", 
		// 	ControlLaw_cr3bp_lt::typeToString(lawType).c_str(),
		// 	MSTOF_tp_cStr(tofType));
		// correctedSet.saveToMat(filename);
	}
}//====================================================

BOOST_DATA_TEST_CASE(test_fullContinuity, 
	data::make(law_bases) * data::make(law_varF) * data::make(tofTypes), 
	baseType, fType, tofType){

	// Construct the ID
	unsigned int lawType = basicLaw | baseType | fType;

	// Set the parameters and states
	std::vector<double> ltParams {1e-1, 3000};	// {fmax, Isp}
	std::vector<double> ctrlState {};
	switch(fType){
		case ltlaw::VAR_F_BND:
			ctrlState.push_back(-0.9273);	// yields f = 1e-2
			break;
		case ltlaw::VAR_F_UBND:
			ctrlState.push_back(0.3162);	// yields f = 1e-2
			break;
		default: break;
	}

	SysData_cr3bp_lt sys("earth", "moon", 14);

	ControlLaw_cr3bp_lt law(lawType, ltParams);

	Arcset_cr3bp_lt halfLyapArcset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL3Lyap_ic, ctrlState, 0, emL3Lyap_T, 3, 
		&halfLyapArcset, &law);

	// Constrain the control to be continuous between the final segment and 
	// final node
	std::vector<double> ctrlContData(1, law.getNumStates());
	Constraint ctrlCont(Constraint_tp::CONT_CTRL, 1, ctrlContData);
	halfLyapArcset.addConstraint(ctrlCont);

	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::NO_MSG);
	corrector.setTOFType(tofType);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapArcset, 
		corrector, Verbosity_tp::NO_MSG, false));

	try{
		corrector.multShoot(&halfLyapArcset, &correctedSet);
	}catch(const Exception &e){
		printErr("Error (%s, %s): %s\n", 
			ControlLaw_cr3bp_lt::typeToString(lawType).c_str(),
			MSTOF_tp_cStr(tofType), e.what());
		char filename[128];
		snprintf(filename, 128, "%s_%s_input.mat", 
			ControlLaw_cr3bp_lt::typeToString(lawType).c_str(),
			MSTOF_tp_cStr(tofType));
		halfLyapArcset.saveToMat(filename);

		snprintf(filename, 128, "%s_%s_corrected.mat", 
			ControlLaw_cr3bp_lt::typeToString(lawType).c_str(),
			MSTOF_tp_cStr(tofType));
		correctedSet.saveToMat(filename);
	}
}//====================================================

BOOST_DATA_TEST_CASE(test_stateConstraint, 
	data::make(law_bases) * data::make(law_varF) * data::make(tofTypes), 
	baseType, fType, tofType){

	// Construct the ID
	unsigned int lawType = basicLaw | baseType | fType;

	// Set the parameters and states
	std::vector<double> ltParams {1e-1, 1500};	// {fmax, Isp}
	std::vector<double> ctrlState {};
	switch(fType){
		case ltlaw::VAR_F_BND:
			ctrlState.push_back(-0.9614);	// yields f = 9e-3
			break;
		case ltlaw::VAR_F_UBND:
			ctrlState.push_back(0.3);	// yields f = 9e-3
			break;
		default: break;
	}

	SysData_cr3bp_lt sys("earth", "moon", 500);
	ControlLaw_cr3bp_lt law(lawType, ltParams);

	Arcset_cr3bp_lt halfLyapArcset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, ctrlState, 0, emL1Lyap_T/2, 3, 
		&halfLyapArcset, &law);

	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::NO_MSG);
	corrector.setTOFType(tofType);
	corrector.setFullFinalProp(false);
	// corrector.setVerbosity(Verbosity_tp::NO_MSG);

	double stateConData[] = {0.78, 0, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 2, stateConData, 7);
	halfLyapArcset.addConstraint(stateCon);

	// Remove final ctrl state; no constraints on final ctrl state, results in 
	// column of zeros
	Constraint rmEndCtrl(Constraint_tp::RM_CTRL, 
		halfLyapArcset.getNodeByIx(-1).getID(), nullptr, 0);
	halfLyapArcset.addConstraint(rmEndCtrl);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapArcset, 
		corrector, Verbosity_tp::NO_MSG));
	BOOST_REQUIRE_NO_THROW(corrector.multShoot(&halfLyapArcset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_DATA_TEST_CASE(test_endSegCon, 
	data::make(law_bases) * data::make(law_varF) * data::make(tofTypes), 
	baseType, fType, tofType){

	// Construct the ID
	unsigned int lawType = basicLaw | baseType | fType;

	// Set the parameters and states
	std::vector<double> ltParams {1e-1, 1500};	// {fmax, Isp}
	std::vector<double> ctrlState {};
	switch(fType){
		case ltlaw::VAR_F_BND:
			ctrlState.push_back(-0.9614);	// yields f = 9e-3
			break;
		case ltlaw::VAR_F_UBND:
			ctrlState.push_back(0.3);	// yields f = 9e-3
			break;
		default: break;
	}

	SysData_cr3bp_lt sys("earth", "moon", 500);
	ControlLaw_cr3bp_lt law(lawType, ltParams);

	Arcset_cr3bp_lt halfLyapArcset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, ctrlState, 0, emL1Lyap_T/2, 3, 
		&halfLyapArcset, &law);
	halfLyapArcset.deleteNode(2);

	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::NO_MSG);
	corrector.setTOFType(tofType);
	corrector.setFullFinalProp(false);
	// corrector.setVerbosity(Verbosity_tp::NO_MSG);

	double stateConData[] = {0.78, 0, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::ENDSEG_STATE, 1, stateConData, 7);
	halfLyapArcset.addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapArcset, 
		corrector, Verbosity_tp::NO_MSG, true));
	BOOST_REQUIRE_NO_THROW(corrector.multShoot(&halfLyapArcset, &correctedSet));

	std::vector<double> fullState = 
		correctedSet.getSegByIx(stateCon.getID()).getStateByRow(-1);
	std::vector<double> finalState(fullState.begin(), fullState.begin()+7);
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_DATA_TEST_CASE(test_rmState, 
	data::make(law_bases) * data::make(law_varF) * data::make(tofTypes), 
	baseType, fType, tofType){

	// Construct the ID
	unsigned int lawType = basicLaw | baseType | fType;

	// Set the parameters and states
	std::vector<double> ltParams {1e-1, 3000};	// {fmax, Isp}
	std::vector<double> ctrlState {};
	switch(fType){
		case ltlaw::VAR_F_BND:
			ctrlState.push_back(-0.9273);	// yields f = 1e-2
			break;
		case ltlaw::VAR_F_UBND:
			ctrlState.push_back(0.3162);	// yields f = 1e-2
			break;
		default: break;
	}

	SysData_cr3bp_lt sys("earth", "moon", 14);
	ControlLaw_cr3bp_lt law(lawType, ltParams);

	Arcset_cr3bp_lt halfLyapArcset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, ctrlState, 0, emL1Lyap_T/2, 3, 
		&halfLyapArcset, &law);

	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::NO_MSG);
	corrector.setTOFType(tofType);
	corrector.setTol(1e-11);
	corrector.setFullFinalProp(false);
	// corrector.setVerbosity(Verbosity_tp::NO_MSG);

	double stateConData[] = {NAN, 0, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 2, stateConData, 7);
	halfLyapArcset.addConstraint(stateCon);

	Constraint rmState(Constraint_tp::RM_STATE, 0, nullptr, 0);
	halfLyapArcset.addConstraint(rmState);

	// Remove final ctrl state; no constraints on final ctrl state, results in 
	// column of zeros
	Constraint rmEndCtrl(Constraint_tp::RM_CTRL, 
		halfLyapArcset.getNodeByIx(-1).getID(), nullptr, 0);
	halfLyapArcset.addConstraint(rmEndCtrl);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&halfLyapArcset, 
		corrector, Verbosity_tp::NO_MSG));
	BOOST_REQUIRE_NO_THROW(corrector.multShoot(&halfLyapArcset, &correctedSet));

	std::vector<double> finalState = correctedSet.getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));

	// correctedSet.saveToMat("data/lt_correctedSet.mat");
}//====================================================

BOOST_AUTO_TEST_SUITE_END()
//------------------------------------------------------------------------------



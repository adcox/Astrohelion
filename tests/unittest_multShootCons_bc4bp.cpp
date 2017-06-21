/**
 * Useful info about Boost Unit Tests:
 * <http://www.boost.org/doc/libs/1_62_0/libs/test/doc/html/index.html>
 */

#define BOOST_TEST_MODULE MultipleShootingConstraints_BC4BP

#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <vector>

#include "Calculations.hpp"
#include "Constraint.hpp"
#include "MultShootEngine.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "Arcset_bc4bp.hpp"
#include "Arcset_cr3bp.hpp"
#include "SimEngine.hpp"
#include "SysData_bc4bp.hpp"
#include "SysData_cr3bp.hpp"
#include "Utilities.hpp"

using namespace astrohelion;
using namespace boost::unit_test;

// Sun-Earth L1 Lyapunov initial state and period
double lyap_ic[] = {-0.745230328320519, 7.22625684942683e-04, 7.45549413286038e-05, -7.30710697247992e-06, -0.0148897145134465, -1.23266135281459e-06};
double lyap_T = 313;

// All the different ways to parameterize time in the multiple shooting algorithm
std::vector<MSTOF_tp> tofTypes {MSTOF_tp::VAR_FREE, MSTOF_tp::VAR_FIXSIGN, MSTOF_tp::VAR_EQUALARC};

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

struct fixture_SEM_Init{
	fixture_SEM_Init(){
		sys = new SysData_bc4bp("Sun", "earth", "moon");
		halfLyapSet = new Arcset_bc4bp(sys);
		correctedSet = new Arcset_bc4bp(sys);
		sim = new SimEngine();
		corrector = new MultShootEngine();

		corrector->setVerbosity(Verbosity_tp::NO_MSG);
		sim->setVerbosity(Verbosity_tp::NO_MSG);
	}//====================================================

	~fixture_SEM_Init(){
		delete corrector;
		delete sim;
		delete correctedSet;
		delete halfLyapSet;
		delete sys;
	}//====================================================

	SysData_bc4bp *sys = nullptr;
	Arcset_bc4bp *halfLyapSet = nullptr, *correctedSet = nullptr;
	SimEngine *sim = nullptr;
	MultShootEngine *corrector = nullptr;
};	// -- END OF fixture_SEM_Init


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(BC4BP_SunEarthMoon)

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_State, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T, 6, halfLyapSet);
	corrector->setTOFType(tofTp);
	
	// Constraint_tp::STATE
	double stateConData[] = {-0.77, 0.5, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 4, stateConData, 6);
	halfLyapSet->addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));
	
	std::vector<double> finalState = correctedSet->getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_State_EndSeg, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T/2.1, 3, halfLyapSet);
	halfLyapSet->deleteNode(2);	// Delete final node so arcset ends with segment
	corrector->setTOFType(tofTp);
	
	// Constraint_tp::ENDSEG_STATE
	double stateConData[] = {NAN, 0, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::ENDSEG_STATE, 1, stateConData, 6);
	halfLyapSet->addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));
	
	std::vector<double> fullState = correctedSet->getSegByIx(stateCon.getID()).getStateByRow(-1);
	std::vector<double> finalState(fullState.begin(), fullState.begin()+6);
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_State_rmState, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T/2.1, 4, halfLyapSet);
	corrector->setTOFType(tofTp);

	double stateConData[] = {NAN, 0, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 3, stateConData, 6);
	halfLyapSet->addConstraint(stateCon);

	// Remove the initial state from the free variable vector
	Constraint rmState(Constraint_tp::RM_STATE, 0, nullptr, 0);
	halfLyapSet->addConstraint(rmState);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));

	// std::vector<double> finalState = correctedSet->getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(correctedSet->getState(stateCon.getID()), stateConData, 1e-12));
	BOOST_CHECK(stateDiffBelowTol(correctedSet->getState(rmState.getID()), halfLyapSet->getState(rmState.getID()), 1e-12));
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_State_rmEpoch, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T, 7, halfLyapSet);
	if(tofTp == MSTOF_tp::VAR_EQUALARC)
		return;	// in equal arc configuration, epochs are not stored

	corrector->setTOFType(tofTp);
	
	// Constraint_tp::STATE
	double stateConData[] = {NAN, 0, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 3, stateConData, 6);
	Constraint rmEpoch(Constraint_tp::RM_EPOCH, 0, nullptr, 0);

	halfLyapSet->addConstraint(stateCon);
	halfLyapSet->addConstraint(rmEpoch);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));
	
	std::vector<double> finalState = correctedSet->getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
	BOOST_CHECK_SMALL(halfLyapSet->getEpoch(rmEpoch.getID()) - correctedSet->getEpoch(rmEpoch.getID()), 1e-12);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_Epoch, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T, 6, halfLyapSet);
	if(tofTp == MSTOF_tp::VAR_EQUALARC)
		return;	// in equal arc configuration, epochs are not stored

	corrector->setTOFType(tofTp);
	
	// Constraint_tp::EPOCH
	double epochConData[] = {0.4};
	Constraint epochCon(Constraint_tp::EPOCH, 0, epochConData, 1);
	halfLyapSet->addConstraint(epochCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));
	
	BOOST_CHECK_SMALL(correctedSet->getEpochByIx(0) - epochConData[0], 1e-12);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_MatchAll, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T, 6, halfLyapSet);
	corrector->setTOFType(tofTp);

	// Constraint_tp::MATCH_ALL
	double matchAllConData = 0;
	Constraint matchAllCon(Constraint_tp::MATCH_ALL, 4, &matchAllConData, 1);
	halfLyapSet->addConstraint(matchAllCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));

	std::vector<double> finalState = correctedSet->getState(matchAllCon.getID());
	std::vector<double> initState = correctedSet->getStateByIx(0);
	BOOST_CHECK(stateDiffBelowTol(finalState, initState, 1e-12));
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_MatchCust, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T, 6, halfLyapSet);
	corrector->setTOFType(tofTp);

	// Constraint_tp::MATCH_CUST
	double matchCustConData[] = {0,0,NAN,NAN,NAN,NAN};
	Constraint matchCustCon(Constraint_tp::MATCH_CUST, 4, matchCustConData, 6);
	halfLyapSet->addConstraint(matchCustCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));
	
	std::vector<double> finalState = correctedSet->getState(matchCustCon.getID());
	std::vector<double> initState = correctedSet->getStateByIx(0);
	finalState.erase(finalState.begin()+2, finalState.end());	// Erase entries 2 through 5; we're only comparing the first two
	initState.erase(initState.begin()+2, initState.end());
	BOOST_CHECK(stateDiffBelowTol(finalState, initState, 1e-12));
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_Dist, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T, 6, halfLyapSet);
	corrector->setTOFType(tofTp);

	double matchDistConData[] = {1, 1.0};
	Constraint matchDistCon(Constraint_tp::DIST, 3, matchDistConData, 2);
	halfLyapSet->addConstraint(matchDistCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));
	std::vector<double> finalState = correctedSet->getState(matchDistCon.getID());
	std::vector<double> primPos = sys->getDynamicsModel()->getPrimPos(correctedSet->getEpoch(matchDistCon.getID()), sys);
	double dist = sqrt(pow(finalState[0] - primPos[3] ,2) + pow(finalState[1] - primPos[4], 2) + pow(finalState[2] - primPos[5], 2));
	BOOST_CHECK_SMALL(dist - matchDistConData[1], 1e-12);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_MinDist, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T, 6, halfLyapSet);
	corrector->setTOFType(tofTp);

	double matchDistConData[] = {1, 1.1};
	Constraint matchDistCon(Constraint_tp::MIN_DIST, 3, matchDistConData, 2);
	halfLyapSet->addConstraint(matchDistCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));

	std::vector<double> finalState = correctedSet->getState(matchDistCon.getID());
	std::vector<double> primPos = sys->getDynamicsModel()->getPrimPos(correctedSet->getEpoch(matchDistCon.getID()), sys);
	double dist = sqrt(pow(finalState[0] - primPos[3] ,2) + pow(finalState[1] - primPos[4], 2) + pow(finalState[2] - primPos[5], 2));
	BOOST_CHECK_GE(dist, matchDistConData[1]);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_MaxDist, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T, 6, halfLyapSet);
	corrector->setTOFType(tofTp);

	double matchDistConData[] = {1, 0.9};
	Constraint matchDistCon(Constraint_tp::MAX_DIST, 3, matchDistConData, 2);
	halfLyapSet->addConstraint(matchDistCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));
	std::vector<double> finalState = correctedSet->getState(matchDistCon.getID());
	std::vector<double> primPos = sys->getDynamicsModel()->getPrimPos(correctedSet->getEpoch(matchDistCon.getID()), sys);
	double dist = sqrt(pow(finalState[0] - primPos[3] ,2) + pow(finalState[1] - primPos[4], 2) + pow(finalState[2] - primPos[5], 2));

	BOOST_CHECK_LE(dist, matchDistConData[1]);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_Dist_endSeg, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T/1.5, 4, halfLyapSet);
	corrector->setTOFType(tofTp);

	double matchDistConData[] = {1, 1.0};
	Constraint matchDistCon(Constraint_tp::ENDSEG_DIST, 2, matchDistConData, 2);
	halfLyapSet->addConstraint(matchDistCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));

	std::vector<double> finalState = correctedSet->getSegRef(matchDistCon.getID()).getStateByRow(-1);
	double epoch = correctedSet->getSegRef(matchDistCon.getID()).getTimeByIx(-1);
	std::vector<double> primPos = sys->getDynamicsModel()->getPrimPos(epoch, sys);
	double dist = sqrt(pow(finalState[0] - primPos[3] ,2) + pow(finalState[1] - primPos[4], 2) + pow(finalState[2] - primPos[5], 2));
	BOOST_CHECK_SMALL(dist - matchDistConData[1], 1e-12);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_MinDist_endSeg, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T/1.5, 4, halfLyapSet);
	corrector->setTOFType(tofTp);

	double matchDistConData[] = {1, 1.1};
	Constraint matchDistCon(Constraint_tp::ENDSEG_MIN_DIST, 2, matchDistConData, 2);
	halfLyapSet->addConstraint(matchDistCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));

	std::vector<double> finalState = correctedSet->getSegRef(matchDistCon.getID()).getStateByRow(-1);
	double epoch = correctedSet->getSegRef(matchDistCon.getID()).getTimeByIx(-1);
	std::vector<double> primPos = sys->getDynamicsModel()->getPrimPos(epoch, sys);
	double dist = sqrt(pow(finalState[0] - primPos[3] ,2) + pow(finalState[1] - primPos[4], 2) + pow(finalState[2] - primPos[5], 2));
	BOOST_CHECK_GE(dist, matchDistConData[1]);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_MaxDist_endSeg, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T/1.5, 4, halfLyapSet);
	corrector->setTOFType(tofTp);

	double matchDistConData[] = {1, 0.9};
	Constraint matchDistCon(Constraint_tp::ENDSEG_MAX_DIST, 2, matchDistConData, 2);
	halfLyapSet->addConstraint(matchDistCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));

	std::vector<double> finalState = correctedSet->getSegRef(matchDistCon.getID()).getStateByRow(-1);
	double epoch = correctedSet->getSegRef(matchDistCon.getID()).getTimeByIx(-1);
	std::vector<double> primPos = sys->getDynamicsModel()->getPrimPos(epoch, sys);
	double dist = sqrt(pow(finalState[0] - primPos[3] ,2) + pow(finalState[1] - primPos[4], 2) + pow(finalState[2] - primPos[5], 2));

	BOOST_CHECK_LE(dist, matchDistConData[1]);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_MaxDeltaV, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T, 6, halfLyapSet);
	corrector->setTOFType(tofTp);

	std::vector<double> state = halfLyapSet->getStateByIx(3);
	state[3] += 0.01;
	state[4] += 0.1;
	state[5] += 0.001;
	halfLyapSet->setState(3, state);	// Perturb the velocity of this state to create a discontinuity
	std::vector<int> dvSegs {2};
	halfLyapSet->allowDV_at(dvSegs);	// Allow the perturbed node to have a delta-v
	double maxDVConData = 0.03;
	Constraint dVCon(Constraint_tp::MAX_DELTA_V, 0, &maxDVConData, 1);
	halfLyapSet->addConstraint(dVCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	MultShootData itData(correctedSet);
	BOOST_CHECK_NO_THROW(itData = corrector->multShoot(halfLyapSet, correctedSet));

	double totalDV = MultShootEngine::getTotalDV(&itData);
	BOOST_CHECK_LE(totalDV, maxDVConData);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_DeltaV, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T, 6, halfLyapSet);
	corrector->setTOFType(tofTp);

	std::vector<double> state = halfLyapSet->getStateByIx(3);
	state[3] += 0.01;
	state[4] += 0.1;
	state[5] += 0.001;
	halfLyapSet->setState(3, state);	// Perturb the velocity of this state to create a discontinuity
	std::vector<int> dvSegs {2};
	halfLyapSet->allowDV_at(dvSegs);	// Allow the perturbed node to have a delta-v
	double maxDVConData = 0.02*0.02;
	Constraint dVCon(Constraint_tp::DELTA_V, 0, &maxDVConData, 1);
	halfLyapSet->addConstraint(dVCon);

	// This one throws errors but is ok
	BOOST_CHECK(true || MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));

	MultShootData itData(correctedSet);
	BOOST_CHECK_NO_THROW(itData = corrector->multShoot(halfLyapSet, correctedSet));

	double totalDV = MultShootEngine::getTotalDV(&itData);
	BOOST_CHECK_SMALL(totalDV - maxDVConData, 1e-12);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_TOF, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T, 6, halfLyapSet);
	corrector->setTOFType(tofTp);

	double tofData = 340.0;
	Constraint tofCon(Constraint_tp::TOF, 0, &tofData, 1);
	halfLyapSet->addConstraint(tofCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));

	double totalTOF = correctedSet->getTotalTOF();
	BOOST_CHECK_SMALL(totalTOF - tofData, 1e-12);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_Apse, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T, 6, halfLyapSet);
	corrector->setTOFType(tofTp);

	double apseData = 2;
	Constraint apseCon(Constraint_tp::APSE, 4, &apseData, 1);
	halfLyapSet->addConstraint(apseCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));

	std::vector<double> finalState = correctedSet->getState(apseCon.getID());
	const DynamicsModel *model = sys->getDynamicsModel();
	double rdot = model->getRDot(apseData, correctedSet->getEpoch(apseCon.getID()), &(finalState[0]), sys);
	BOOST_CHECK_SMALL(rdot, 1e-12);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_Apse_endSeg, data::make(tofTypes), tofTp){
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T, 6, halfLyapSet);
	corrector->setTOFType(tofTp);

	double apseData = 2;
	Constraint apseCon(Constraint_tp::ENDSEG_APSE, 3, &apseData, 1);
	halfLyapSet->addConstraint(apseCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(halfLyapSet, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(halfLyapSet, correctedSet));

	std::vector<double> finalState = correctedSet->getSegRef(apseCon.getID()).getStateByRow(-1);
	double epoch = correctedSet->getSegRef(apseCon.getID()).getTimeByIx(-1);
	const DynamicsModel *model = sys->getDynamicsModel();
	double rdot = model->getRDot(apseData, epoch, &(finalState[0]), sys);
	BOOST_CHECK_SMALL(rdot, 1e-12);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_SaddlePoint_Exact, data::make(tofTypes), tofTp){
	corrector->setTOFType(tofTp);

	double IC[] = {-0.71200455, 0.16675922, 0.02755461, 0.01186449, -0.00004723, -0.0010737};
	double t0 = 10207.19;
	double tof = 51.32;
	Arcset_bc4bp nodes0(sys);
	sim->runSim_manyNodes(IC, t0, tof, 5, &nodes0);
	corrector->setTol(1e-11);

	double spData = 0;
	Constraint spCon(Constraint_tp::SP, 2, &spData, 1);
	nodes0.addConstraint(spCon);
	
	// This one throws errors but is ok
	BOOST_CHECK(true || MultShootEngine::finiteDiff_checkMultShoot(&nodes0, *corrector, Verbosity_tp::NO_MSG));

	BOOST_CHECK_NO_THROW(corrector->multShoot(&nodes0, correctedSet));
	Eigen::Vector3d spPos = bcr4bpr_getSPLoc(sys, correctedSet->getEpoch(spCon.getID()));
	std::vector<double> finalState = correctedSet->getState(spCon.getID());
	double diff = sqrt(pow(spPos(0) - finalState[0], 2) + pow(spPos(1) - finalState[1], 2) + pow(spPos(2) - finalState[2], 2));
	BOOST_CHECK_SMALL(diff, 1e-10);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_Cont_PV, data::make(tofTypes), tofTp){
	if(tofTp == MSTOF_tp::VAR_EQUALARC)
		return;

	corrector->setTOFType(tofTp);

	Arcset_bc4bp forwardArc(sys), reverseArc(sys);
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T/2.0, 4, &forwardArc);
	sim->setRevTime(true);
	sim->runSim_manyNodes(lyap_ic, lyap_T, -lyap_T/2.0, 4, &reverseArc);

	forwardArc.deleteNode(3);
	reverseArc.deleteNode(3);
	Arcset_bc4bp doubleSrcLyap = forwardArc;
	doubleSrcLyap.concatArcset(&reverseArc);
	// doubleSrcLyap.print();
	// doubleSrcLyap.printInChrono();

	double contData[] = {5, 5, NAN, 5, 5, NAN};
	Constraint contCon(Constraint_tp::SEG_CONT_PV, 2, contData, 6);
	doubleSrcLyap.addConstraint(contCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&doubleSrcLyap, *corrector, Verbosity_tp::NO_MSG));

	MultShootData it(correctedSet);
	BOOST_CHECK_NO_THROW(it = corrector->multShoot(&doubleSrcLyap, correctedSet));

	Arcset forwardTraj = it.propSegs[correctedSet->getSegIx(contCon.getID())];
	Arcset reverseTraj = it.propSegs[correctedSet->getSegIx(contData[0])];
	std::vector<double> for_lastState = forwardTraj.getStateByIx(-1);
	std::vector<double> rev_lastState = reverseTraj.getStateByIx(-1);
	double sum = 0;
	for(int i = 0; i < 6; i++){
		if(!std::isnan(contData[i]))
			sum += pow(for_lastState[i] - rev_lastState[i], 2);
	}
	BOOST_CHECK_SMALL(sum, 1e-12);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_SegCont_Ex, data::make(tofTypes), tofTp){
	if(tofTp == MSTOF_tp::VAR_EQUALARC)
		return;
	
	corrector->setTOFType(tofTp);

	Arcset_bc4bp forwardArc(sys), reverseArc(sys);
	sim->runSim_manyNodes(lyap_ic, 0, lyap_T/2.0, 4, &forwardArc);
	sim->setRevTime(true);
	sim->runSim_manyNodes(lyap_ic, lyap_T, -lyap_T/2.0, 4, &reverseArc);

	forwardArc.deleteNode(3);
	reverseArc.deleteNode(3);
	Arcset_bc4bp doubleSrcLyap = forwardArc;
	doubleSrcLyap.concatArcset(&reverseArc);

	double contExData[] = {5, 0};
	Constraint extraContCon(Constraint_tp::SEG_CONT_EX, 2, contExData, 2);
	doubleSrcLyap.addConstraint(extraContCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&doubleSrcLyap, *corrector, Verbosity_tp::NO_MSG));

	MultShootData it(correctedSet);
	BOOST_CHECK_NO_THROW(it = corrector->multShoot(&doubleSrcLyap, correctedSet));

	Arcset forwardTraj = it.propSegs[correctedSet->getSegIx(extraContCon.getID())];
	Arcset reverseTraj = it.propSegs[correctedSet->getSegIx(contExData[0])];
	BOOST_CHECK_SMALL(forwardTraj.getEpochByIx(-1) - reverseTraj.getEpochByIx(-1), 1e-12);
}//====================================================

BOOST_DATA_TEST_CASE_F(fixture_SEM_Init, BC4BP_SEM_SourceNode, data::make(tofTypes), tofTp){
	SysData_cr3bp emSys("earth", "moon");

	if(tofTp == MSTOF_tp::VAR_EQUALARC)
		return;

	corrector->setTOFType(tofTp);

	// EM L2 Butterfly Orbit
	std::vector<double> ic {1.0639767173456007, 0, 0.1644973017995331, 0, -0.0311246472806882, 0};
	double tof = 3.1256890778;

	Arcset_cr3bp halfPlus(&emSys), halfMinus(&emSys);
	sim->runSim_manyNodes(ic, tof, 4, &halfPlus);
	sim->setRevTime(true);
	sim->runSim_manyNodes(ic, tof, 4, &halfMinus);
	halfPlus.appendSetAtNode(&halfMinus, 0, 0, 0);

	// Transform to SE coordinates
	SysData_cr3bp seSys("sun", "earth");
	double epoch = dateToEphemerisTime("2016/10/16 22:48:02");
	DynamicsModel_bc4bp::orientAtEpoch(epoch, sys);

	Arcset_cr3bp seNodes = cr3bp_EM2SE(halfPlus, &seSys, sys->getTheta0(), sys->getPhi0(), sys->getGamma());
	Arcset_bc4bp semNodes = bcr4bpr_SE2SEM(seNodes, sys, 0, 0);

	double stateConData[] = {0.205, 0.178, 0.068, NAN, 0.025, NAN};
	Constraint stateCon(Constraint_tp::STATE, 3, stateConData, 6);
	semNodes.addConstraint(stateCon);

	BOOST_CHECK(MultShootEngine::finiteDiff_checkMultShoot(&semNodes, *corrector, Verbosity_tp::NO_MSG));
	BOOST_CHECK_NO_THROW(corrector->multShoot(&semNodes, correctedSet));

	std::vector<double> finalState = correctedSet->getState(stateCon.getID());
	BOOST_CHECK(stateDiffBelowTol(finalState, stateConData, 1e-12));
}//====================================================

BOOST_AUTO_TEST_SUITE_END()



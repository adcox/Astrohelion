#define BOOST_TEST_MODULE LambertArcEngine

#include <boost/test/unit_test.hpp>
/**
 *	Test functions in the Lambert Engine
 *	
 *	TODO:
 *	* Check 2A Lambert Solver case
 *  * Write 2B Lambert Solver test case (532 HW 11.1)
 *  * Develop 1B, 1H, and 2H test cases for Lambert Solver
 */

#include "AsciiOutput.hpp"
#include "LambertArcEngine.hpp"
#include "Node.hpp"
#include "SimEngine.hpp"
#include "SysData_2bp.hpp"
#include "Arcset_2bp.hpp"

#include <cmath>
#include <iostream>
#include <vector>

using namespace astrohelion;

BOOST_AUTO_TEST_CASE(TYPE_1A){
	SysData_2bp sys("earth");
	std::vector<double> r1 {9567.2055, 0, 0};
	std::vector<double> r2 {10824.0574176883, -54120.2870884419, -26513.4176199238};
	double tof = 8*3600;

	LambertArcEngine engine;
	engine.setVerbosity(Verbosity_tp::NO_MSG);
	Arcset_2bp type1_traj = engine.getLambertArc(&sys, r1, r2, tof, 1);

	Node node = type1_traj.getNodeByIx(0);
	Node node2 = type1_traj.getNodeByIx(1);

	// node.print();
	BOOST_CHECK(std::abs(node.getExtraParam(PARAMKEY_SMA) - 32773.3694) < 1e-3);
	BOOST_CHECK(std::abs(node.getExtraParam(PARAMKEY_ECC) - 0.8685) < 1e-3);
	BOOST_CHECK(std::abs(node.getExtraParam("ta") - 100.4889*PI/180) < 1e-3);
	BOOST_CHECK(std::abs(node.getExtraParam("fpa") - 45.4072*PI/180) < 1e-3);
	BOOST_CHECK(std::abs(node.getExtraParam("speed") - 8.4359) < 1e-3);
	BOOST_CHECK(std::abs(node2.getExtraParam("ta") - 180.3068*PI/180) < 1e-3);
	BOOST_CHECK(std::abs(node2.getExtraParam("fpa") - -2.0249*PI/180) < 1e-3);
	BOOST_CHECK(std::abs(node2.getExtraParam("speed") - 0.9260) < 1e-3);

	Arcset_2bp type1_traj_sim(&sys);
	SimEngine sim;
	sim.setMakeDefaultEvents(false);
	sim.runSim(node.getState(), node2.getEpoch() - node.getEpoch(), &type1_traj_sim);

	std::vector<double> finalState = type1_traj_sim.getStateByIx(-1);
	double stateDiff[] = {	finalState[0] - node2.getState()[0],
							finalState[1] - node2.getState()[1],
							finalState[2] - node2.getState()[2],
							finalState[3] - node2.getState()[3],
							finalState[4] - node2.getState()[4],
							finalState[5] - node2.getState()[5]};
	// printf("State Diff = [%.2e, %.2e, %.2e, %.2e, %.2e, %.2e]\n", stateDiff[0], stateDiff[1], stateDiff[2],
	// 	stateDiff[3], stateDiff[4], stateDiff[5]);

	BOOST_CHECK(std::abs(stateDiff[0]) < 1e-3);
	BOOST_CHECK(std::abs(stateDiff[1]) < 1e-3);
	BOOST_CHECK(std::abs(stateDiff[2]) < 1e-3);
	BOOST_CHECK(std::abs(stateDiff[3]) < 1e-3);
	BOOST_CHECK(std::abs(stateDiff[4]) < 1e-3);
	BOOST_CHECK(std::abs(stateDiff[5]) < 1e-3);
}//====================================================

BOOST_AUTO_TEST_CASE(TYPE_2A){
	SysData_2bp sys("earth");
	std::vector<double> r1 {9567.2055, 0, 0};
	std::vector<double> r2 {10824.0574176883, -54120.2870884419, -26513.4176199238};
	double tof = 8*3600;

	LambertArcEngine engine;
	engine.setVerbosity(Verbosity_tp::NO_MSG);

	Arcset_2bp type2_traj = engine.getLambertArc(&sys, r1, r2, tof, 2);

	Node node = type2_traj.getNodeByIx(0);
	Node node2 = type2_traj.getNodeByIx(1);

	BOOST_CHECK(std::abs(node.getExtraParam(PARAMKEY_SMA) - 32794.2790) < 1e-3);
	BOOST_CHECK(std::abs(node.getExtraParam(PARAMKEY_ECC) - 0.8718) < 1e-3);
	BOOST_CHECK(std::abs(node.getExtraParam("ta") - 258.2606*PI/180) < 1e-3);
	BOOST_CHECK(std::abs(node.getExtraParam("fpa") - -46.0570*PI/180) < 1e-3);
	BOOST_CHECK(std::abs(node.getExtraParam("speed") - 8.4363) < 1e-3);
	BOOST_CHECK(std::abs(node2.getExtraParam("ta") - 178.4427*PI/180) < 1e-3);
	BOOST_CHECK(std::abs(node2.getExtraParam("fpa") - 10.4440*PI/180) < 1e-3);
	BOOST_CHECK(std::abs(node2.getExtraParam("speed") - 0.9301) < 1e-3);

	Arcset_2bp type2_traj_sim(&sys);
	SimEngine sim;
	sim.setMakeDefaultEvents(false);
	// sim.setVerbosity(Verbosity_tp::ALL_MSG);
	sim.runSim(node.getState(), node2.getEpoch() - node.getEpoch(), &type2_traj_sim);

	std::vector<double> finalState = type2_traj_sim.getStateByIx(-1);
	double stateDiff[] = {	finalState[0] - node2.getState()[0],
							finalState[1] - node2.getState()[1],
							finalState[2] - node2.getState()[2],
							finalState[3] - node2.getState()[3],
							finalState[4] - node2.getState()[4],
							finalState[5] - node2.getState()[5]};
	printf("State Diff = [%.2e, %.2e, %.2e, %.2e, %.2e, %.2e]\n", stateDiff[0], stateDiff[1], stateDiff[2],
		stateDiff[3], stateDiff[4], stateDiff[5]);

	BOOST_CHECK(std::abs(stateDiff[0]) < 1e-3);
	BOOST_CHECK(std::abs(stateDiff[1]) < 1e-3);
	BOOST_CHECK(std::abs(stateDiff[2]) < 1e-3);
	BOOST_CHECK(std::abs(stateDiff[3]) < 1e-3);
	BOOST_CHECK(std::abs(stateDiff[4]) < 1e-3);
	BOOST_CHECK(std::abs(stateDiff[5]) < 1e-3);
}//====================================================
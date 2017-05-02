#define BOOST_TEST_MODULE LinMotion

#include <boost/test/unit_test.hpp>

#include <iostream>

#include "AsciiOutput.hpp"
#include "Calculations.hpp"
#include "SysData_cr3bp.hpp"
#include "Arcset_cr3bp.hpp"
#include "LinMotionEngine.hpp"
#include "Utilities.hpp"

using namespace std;
using namespace astrohelion;

// static const char* PASS = BOLDGREEN "PASS" RESET;
// static const char* FAIL = BOLDRED "FAIL" RESET;

LinMotionEngine engine;
double tol = 1e-12;

//************************************************************
//* Linkable Tests
//************************************************************

BOOST_AUTO_TEST_SUITE(Linear_Motion)

BOOST_AUTO_TEST_CASE(Linear_Nonzero_Time){
	double r[] = {0.1,0.1,0};

	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp traj = engine.getCR3BPLinear(1, r, LinMotion_tp::ELLIP, &sys);

	BOOST_CHECK(traj.getTimeByIx(-1) > 0);
}//====================================================

/*
 * For these tests to work, t_step should be set to 0.001, and the body data must match what I have in MATLAB
 */
BOOST_AUTO_TEST_CASE(L1_Elliptical_Linear){
	// printf("Testing Elliptical Linearization near L1:\n");
	double r[] = {0.1,0.1,0};
	double q0[] = {0.936915132364302, 0.1, 0, 0.065088146132374, -0.837227319488528, 0};
	double q9[] = {0.937478813500055, 0.092443439352275};

	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp traj = engine.getCR3BPLinear(1, r, LinMotion_tp::ELLIP, &sys);
	
	for(int i = 0; i < 6; i++){
		BOOST_CHECK(abs(traj.getStateByIx(0)[i] - q0[i]) < tol);
		// cout << "  state(0, " << i << "): " << ( abs(traj.getStateByIx(0)[i] - q0[i]) < tol ? PASS : FAIL) << endl;
	}

	BOOST_CHECK(abs(traj.getStateByIx(9)[0] - q9[0]) < tol);
	BOOST_CHECK(abs(traj.getStateByIx(9)[1] - q9[1]) < tol);
	// cout << "  state(9, " << 0 << "): " << ( abs(traj.getStateByIx(9)[0] - q9[0]) < tol ? PASS : FAIL) << endl;
	// cout << "  state(9, " << 1 << "): " << ( abs(traj.getStateByIx(9)[1] - q9[1]) < tol ? PASS : FAIL) << endl;
}//====================================================

BOOST_AUTO_TEST_CASE(L2_Hyperbolic_Linear){
	// printf("Testing Hyperbolic Linearization near L2:\n");
	double r[] = {0.1,0.1,0};
	double q0[] = {1.25568216029234, 0.1, 0, -0.342515005674306, -0.136048780251457, 0};
	double q9[] = {1.25261820440111, 0.098794357035758};

	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp traj = engine.getCR3BPLinear(2, r, LinMotion_tp::HYP, &sys);
	
	for(int i = 0; i < 6; i++){
		BOOST_CHECK(abs(traj.getStateByIx(0)[i] - q0[i]) < tol);
		// cout << "  state(0, " << i << "): " << ( abs(traj.getStateByIx(0)[i] - q0[i]) < tol ? PASS : FAIL) << endl;
	}

	BOOST_CHECK(abs(traj.getStateByIx(9)[0] - q9[0]) < tol);
	BOOST_CHECK(abs(traj.getStateByIx(9)[1] - q9[1]) < tol);
	// cout << "  state(9, " << 0 << "): " << ( abs(traj.getStateByIx(9)[0] - q9[0]) < tol ? PASS : FAIL) << endl;
	// cout << "  state(9, " << 1 << "): " << ( abs(traj.getStateByIx(9)[1] - q9[1]) < tol ? PASS : FAIL) << endl;
}//====================================================

BOOST_AUTO_TEST_CASE(L4_SPO_Linear){
	// printf("Testing SPO Linearization near L4:\n");
	double r[] = {0.1,0.1,0};
	double q0[] = {0.58784941573006, 0.966025403784439, 0, 0.221427092899258, -0.146427092899258, 0};
	double q9[] = {0.589838545236831, 0.964703886338579};

	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp traj = engine.getCR3BPLinear(4, r, LinMotion_tp::SPO, &sys);

	for(int i = 0; i < 6; i++){
		BOOST_CHECK(abs(traj.getStateByIx(0)[i] - q0[i]) < tol);
		// cout << "  state(0, " << i << "): " << ( abs(traj.getStateByIx(0)[i] - q0[i]) < tol ? PASS : FAIL) << endl;
	}

	BOOST_CHECK(abs(traj.getStateByIx(9)[0] - q9[0]) < tol);
	BOOST_CHECK(abs(traj.getStateByIx(9)[1] - q9[1]) < tol);
	// cout << "  state(9, " << 0 << "): " << ( abs(traj.getStateByIx(9)[0] - q9[0]) < tol ? PASS : FAIL) << endl;
	// cout << "  state(9, " << 1 << "): " << ( abs(traj.getStateByIx(9)[1] - q9[1]) < tol ? PASS : FAIL) << endl;
}//====================================================

BOOST_AUTO_TEST_CASE(L5_LPO_Linear){
	// printf("Testing LPO Linearization near L5:\n");
	double r[] = {0.1,0.1,0};
	double q0[] = {0.58784941573006, -0.766025403784439, 0, 0.0535729071007422, 0.0214270928992578, 0};
	double q9[] = {0.58833121115652, -0.765832920338464};

	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp traj = engine.getCR3BPLinear(5, r, LinMotion_tp::LPO, &sys);
	
	for(int i = 0; i < 6; i++){
		BOOST_CHECK(abs(traj.getStateByIx(0)[i] - q0[i]) < tol);
		// cout << "  state(0, " << i << "): " << ( abs(traj.getStateByIx(0)[i] - q0[i]) < tol ? PASS : FAIL) << endl;
	}

	BOOST_CHECK(abs(traj.getStateByIx(9)[0] - q9[0]) < tol);
	BOOST_CHECK(abs(traj.getStateByIx(9)[1] - q9[1]) < tol);
	// cout << "  state(9, " << 0 << "): " << ( abs(traj.getStateByIx(9)[0] - q9[0]) < tol ? PASS : FAIL) << endl;
	// cout << "  state(9, " << 1 << "): " << ( abs(traj.getStateByIx(9)[1] - q9[1]) < tol ? PASS : FAIL) << endl;
}//====================================================

BOOST_AUTO_TEST_CASE(L5_MPO_Linear){
	// printf("Testing MPO Linearization near L5:\n");
	double r[] = {0.1,0.1,0};
	double q0[] = {0.58784941573006, -0.766025403784439, 0, 0.0770913737661744, 0.000873498086544126, 0};
	double q9[] = {0.58854120640494, -0.766019806926301};

	SysData_cr3bp sys("earth", "moon");
	Arcset_cr3bp traj = engine.getCR3BPLinear(5, r, LinMotion_tp::MPO, &sys);
	
	for(int i = 0; i < 6; i++){
		BOOST_CHECK(abs(traj.getStateByIx(0)[i] - q0[i]) < tol);
		// cout << "  state(0, " << i << "): " << ( abs(traj.getStateByIx(0)[i] - q0[i]) < tol ? PASS : FAIL) << endl;
	}

	BOOST_CHECK(abs(traj.getStateByIx(9)[0] - q9[0]) < tol);
	BOOST_CHECK(abs(traj.getStateByIx(9)[1] - q9[1]) < tol);
	// cout << "  state(9, " << 0 << "): " << ( abs(traj.getStateByIx(9)[0] - q9[0]) < tol ? PASS : FAIL) << endl;
	// cout << "  state(9, " << 1 << "): " << ( abs(traj.getStateByIx(9)[1] - q9[1]) < tol ? PASS : FAIL) << endl;
}//====================================================

BOOST_AUTO_TEST_CASE(L4_Convergent_Linear){
	// printf("Testing Convergent Linearization near L4:\n");

	double r[] = {0.1,0.1,0};
	double q0[] = {0.491460111646493, 0.966025403784439, 0, 0.132761907938009, -0.180873270208478, 0};
	double q9[] = {0.492647479254839, 0.964400034052942};

	SysData_cr3bp sys("pluto", "charon");
	Arcset_cr3bp traj = engine.getCR3BPLinear(4, r, LinMotion_tp::CONVERGE, &sys);
	
	for(int i = 0; i < 6; i++){
		BOOST_CHECK(abs(traj.getStateByIx(0)[i] - q0[i]) < tol);
		// cout << "  state(0, " << i << "): " << ( abs(traj.getStateByIx(0)[i] - q0[i]) < tol ? PASS : FAIL) << endl;
	}

	BOOST_CHECK(abs(traj.getStateByIx(9)[0] - q9[0]) < tol);
	BOOST_CHECK(abs(traj.getStateByIx(9)[1] - q9[1]) < tol);
	// cout << "  state(9, " << 0 << "): " << ( abs(traj.getStateByIx(9)[0] - q9[0]) < tol ? PASS : FAIL) << endl;
	// cout << "  state(9, " << 1 << "): " << ( abs(traj.getStateByIx(9)[1] - q9[1]) < tol ? PASS : FAIL) << endl;
}//====================================================

BOOST_AUTO_TEST_SUITE_END()
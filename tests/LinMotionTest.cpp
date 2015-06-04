#include <iostream>

#include "adtk_ascii_output.hpp"
#include "adtk_calculations.hpp"
#include "adtk_cr3bp_sys_data.hpp"
#include "adtk_cr3bp_traj.hpp"
#include "adtk_linear_motion_engine.hpp"
#include "adtk_utilities.hpp"

using namespace std;

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;
adtk_linear_motion_engine engine;
double tol = 1e-12;

/*
 * For these tests to work, t_step should be set to 0.001, and the body data must match what I have in MATLAB
 */
void testEllip(){
	printf("Testing Elliptical Linearization near L1:\n");
	double r[] = {0.1,0.1,0};
	double q0[] = {0.936915132364302, 0.1, 0, 0.065088146132374, -0.837227319488528, 0};
	double q9[] = {0.937478813500055, 0.092443439352275};

	adtk_cr3bp_traj traj = engine.getCR3BPLinear(1, r, adtk_linear_motion_engine::ELLIP,
		"earth", "moon");
	
	for(int i = 0; i < 6; i++){
		cout << "  state(0, " << i << "): " << ( abs(traj.getState(0)[i] - q0[i]) < tol ? PASS : FAIL) << endl;
	}

	cout << "  state(9, " << 0 << "): " << ( abs(traj.getState(9)[0] - q9[0]) < tol ? PASS : FAIL) << endl;
	cout << "  state(9, " << 1 << "): " << ( abs(traj.getState(9)[1] - q9[1]) < tol ? PASS : FAIL) << endl;
}

void testHyp(){
	printf("Testing Hyperbolic Linearization near L2:\n");
	double r[] = {0.1,0.1,0};
	double q0[] = {1.25568216029234, 0.1, 0, -0.342515005674306, -0.136048780251457, 0};
	double q9[] = {1.25261820440111, 0.098794357035758};

	adtk_cr3bp_traj traj = engine.getCR3BPLinear(2, r, adtk_linear_motion_engine::HYP,
		"earth", "moon");
	
	for(int i = 0; i < 6; i++){
		cout << "  state(0, " << i << "): " << ( abs(traj.getState(0)[i] - q0[i]) < tol ? PASS : FAIL) << endl;
	}

	cout << "  state(9, " << 0 << "): " << ( abs(traj.getState(9)[0] - q9[0]) < tol ? PASS : FAIL) << endl;
	cout << "  state(9, " << 1 << "): " << ( abs(traj.getState(9)[1] - q9[1]) < tol ? PASS : FAIL) << endl;
}

void testSPO(){
	printf("Testing SPO Linearization near L4:\n");
	double r[] = {0.1,0.1,0};
	double q0[] = {0.58784941573006, 0.966025403784439, 0, 0.221427092899258, -0.146427092899258, 0};
	double q9[] = {0.589838545236831, 0.964703886338579};

	adtk_cr3bp_traj traj = engine.getCR3BPLinear(4, r, adtk_linear_motion_engine::SPO,
		"earth", "moon");

	for(int i = 0; i < 6; i++){
		cout << "  state(0, " << i << "): " << ( abs(traj.getState(0)[i] - q0[i]) < tol ? PASS : FAIL) << endl;
	}

	cout << "  state(9, " << 0 << "): " << ( abs(traj.getState(9)[0] - q9[0]) < tol ? PASS : FAIL) << endl;
	cout << "  state(9, " << 1 << "): " << ( abs(traj.getState(9)[1] - q9[1]) < tol ? PASS : FAIL) << endl;
}

void testLPO(){
	printf("Testing LPO Linearization near L5:\n");
	double r[] = {0.1,0.1,0};
	double q0[] = {0.58784941573006, -0.766025403784439, 0, 0.0535729071007422, 0.0214270928992578, 0};
	double q9[] = {0.58833121115652, -0.765832920338464};

	adtk_cr3bp_traj traj = engine.getCR3BPLinear(5, r, adtk_linear_motion_engine::LPO,
		"earth", "moon");
	
	for(int i = 0; i < 6; i++){
		cout << "  state(0, " << i << "): " << ( abs(traj.getState(0)[i] - q0[i]) < tol ? PASS : FAIL) << endl;
	}

	cout << "  state(9, " << 0 << "): " << ( abs(traj.getState(9)[0] - q9[0]) < tol ? PASS : FAIL) << endl;
	cout << "  state(9, " << 1 << "): " << ( abs(traj.getState(9)[1] - q9[1]) < tol ? PASS : FAIL) << endl;
}

void testMPO(){
	printf("Testing MPO Linearization near L5:\n");
	double r[] = {0.1,0.1,0};
	double q0[] = {0.58784941573006, -0.766025403784439, 0, 0.0770913737661744, 0.000873498086544126, 0};
	double q9[] = {0.58854120640494, -0.766019806926301};

	adtk_cr3bp_traj traj = engine.getCR3BPLinear(5, r, adtk_linear_motion_engine::MPO,
		"earth", "moon");
	
	for(int i = 0; i < 6; i++){
		cout << "  state(0, " << i << "): " << ( abs(traj.getState(0)[i] - q0[i]) < tol ? PASS : FAIL) << endl;
	}

	cout << "  state(9, " << 0 << "): " << ( abs(traj.getState(9)[0] - q9[0]) < tol ? PASS : FAIL) << endl;
	cout << "  state(9, " << 1 << "): " << ( abs(traj.getState(9)[1] - q9[1]) < tol ? PASS : FAIL) << endl;
}

void testConverge(){
	printf("Testing Convergent Linearization near L4:\n");

	double r[] = {0.1,0.1,0};
	double q0[] = {0.491460111646493, 0.966025403784439, 0, 0.132761907938009, -0.180873270208478, 0};
	double q9[] = {0.492647479254839, 0.964400034052942};

	adtk_cr3bp_traj traj = engine.getCR3BPLinear(4, r, adtk_linear_motion_engine::CONVERGE,
		"pluto", "charon");
	
	for(int i = 0; i < 6; i++){
		cout << "  state(0, " << i << "): " << ( abs(traj.getState(0)[i] - q0[i]) < tol ? PASS : FAIL) << endl;
	}

	cout << "  state(9, " << 0 << "): " << ( abs(traj.getState(9)[0] - q9[0]) < tol ? PASS : FAIL) << endl;
	cout << "  state(9, " << 1 << "): " << ( abs(traj.getState(9)[1] - q9[1]) < tol ? PASS : FAIL) << endl;
}

int main(void){

	testEllip();
	testHyp();
	testSPO();
	testLPO();
	testMPO();
	testConverge();
	return 0;
}
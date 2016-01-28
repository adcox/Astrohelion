
#include "tpat_ascii_output.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_bcr4bp.hpp"
#include "tpat_traj_cr3bp.hpp"

#include <iostream>

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

using namespace std;

void testCR3BPTraj(){
	tpat_sys_data_cr3bp emData("earth", "moon");
	tpat_simulation_engine sim(&emData);
	double ic[] = {0.887415132364297, 0, 0, 0, -0.332866299501083, 0};	// EM L1
	double T = 3.02796323553149;	// EM L1 Period
	sim.runSim(ic, T);
	tpat_traj_cr3bp crTraj = sim.getCR3BP_Traj();
	crTraj.saveToMat("data/crTraj.mat");

	tpat_traj_cr3bp crTemp(&emData);
	crTemp.readFromMat("data/crTraj.mat");

	printf("Testing Save/Read functions on CR3BP Trajectory\n");
	cout << "Same Final State: " << (crTraj.getState(-1) == crTemp.getState(-1) ? PASS : FAIL) << endl;
	cout << "Same Final Accel: " << (crTraj.getAccel(-1) == crTemp.getAccel(-1) ? PASS : FAIL) << endl;
	cout << "Same Final Time: " << (crTraj.getTime(-1) == crTemp.getTime(-1) ? PASS : FAIL) << endl;
	cout << "Same Final STM: " << (crTraj.getSTM(-1) == crTemp.getSTM(-1) ? PASS : FAIL) << endl;
	cout << "Same Final Jacobi: " << (crTraj.getJacobi(-1) == crTemp.getJacobi(-1) ? PASS : FAIL) << endl;
}

void testBC4BPTraj(){
	tpat_sys_data_bcr4bpr semData("sun", "earth", "moon");
	tpat_simulation_engine sim(&semData);
	double ic[] = {-0.745230328320519, 7.22625684942683e-04, 7.45549413286038e-05, -7.30710697247992e-06, -0.0148897145134465, -1.23266135281459e-06};
	double T = 313;	// SE L1 Period
	sim.runSim(ic, T);
	tpat_traj_bcr4bp bcTraj = sim.getBCR4BPR_Traj();
	bcTraj.saveToMat("data/bcTraj.mat");

	tpat_traj_bcr4bp bcTemp(&semData);
	bcTemp.readFromMat("data/bcTraj.mat");

	printf("Testing Save/Read functions on CR3BP Trajectory\n");
	cout << "Same Final State: " << (bcTraj.getState(-1) == bcTemp.getState(-1) ? PASS : FAIL) << endl;
	cout << "Same Final Accel: " << (bcTraj.getAccel(-1) == bcTemp.getAccel(-1) ? PASS : FAIL) << endl;
	cout << "Same Final Time: " << (bcTraj.getTime(-1) == bcTemp.getTime(-1) ? PASS : FAIL) << endl;
	cout << "Same Final STM: " << (bcTraj.getSTM(-1) == bcTemp.getSTM(-1) ? PASS : FAIL) << endl;
	cout << "Same Final dqdT: " << (bcTraj.get_dqdT(-1) == bcTemp.get_dqdT(-1) ? PASS : FAIL) << endl;
}


int main(){
	
	testCR3BPTraj();
	testBC4BPTraj();

	return EXIT_SUCCESS;
}
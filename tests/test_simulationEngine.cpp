/**
 *	Test the simulation engine
 */

#include "tpat_ascii_output.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_traj_bcr4bp.hpp"
#include "tpat_constants.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_utilities.hpp"

#include <iostream>

using namespace std;

void test_cr3bp_sim(){
	tpat_sys_data_cr3bp sys("earth", "moon");
	tpat_simulation_engine simEngine;
	tpat_traj_cr3bp traj(&sys);

	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	simEngine.setVerbose(NO_MSG);
	simEngine.setAbsTol(1e-14);
	simEngine.setRelTol(1e-16);
	simEngine.runSim(ic, 2.77, &traj);
	
	

	cout << "Trajectory contains " << traj.getNumNodes() << " points" << endl;
	
	MatrixXRd lastSTM = traj.getSTMByIx(-1);
	std::cout << lastSTM << std::endl;
	// lastSTM.print("%12.4f");

	traj.saveToMat("data/HaloTest.mat");
}//====================================================

void test_bcr4bpr_sim(){
	// Do a simulation in the BCR4BP
	tpat_sys_data_bcr4bpr bcSys("sun", "earth", "moon");
	tpat_simulation_engine bcEngine;
	tpat_traj_bcr4bp bcTraj(&bcSys);

	double haloCross177_IC[] = {0.0359099212, 		0, 					285.85225655914e-05, 
								0.0405130514527453, 0.114452371350509, 	0.0310985198400586};
	double t0 = 2.57;
	
	bcEngine.setVarStepSize(false);
	bcEngine.setNumSteps(500);
	bcEngine.setVerbose(ALL_MSG);
	bcEngine.runSim(haloCross177_IC, t0, 2*PI, &bcTraj);

	cout << "Trajectory contains " << bcTraj.getNumNodes() << " points" << endl;
	bcTraj.saveToMat("data/bcHaloManifoldProp.mat");
}//====================================================

void test_cr3bp_events(){
	tpat_sys_data_cr3bp sys("earth", "moon");
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};	// L1 Halo

	tpat_simulation_engine engine;
	tpat_traj_cr3bp traj(&sys);
	engine.setVerbose(ALL_MSG);
	engine.addEvent(tpat_event(&sys, tpat_event::XZ_PLANE, 0, true));
	engine.setRevTime(true);
	engine.runSim(ic, 2.7, &traj);

	traj.saveToMat("data/HaloHalfTest.mat");
}//====================================================

void test_bcr4bpr_events(){
	tpat_sys_data_bcr4bpr sys("sun", "earth", "moon");
	double ic[] = {0.0359099212, 		0, 	285.85225655914e-05, 
					0.0405130514527453, 0, 	0.0310985198400586};
	double t0 = 2.57;

	tpat_simulation_engine engine;
	tpat_traj_bcr4bp traj(&sys);
	engine.setVerbose(ALL_MSG);
	engine.addEvent(tpat_event(&sys, tpat_event::XY_PLANE, 0, true));
	engine.runSim(ic, t0, 2*PI, &traj);
	traj.saveToMat("data/BC_HaloManifold.mat");
}//====================================================

int main(void){
	printColor(RED, "*************************\n* Test CR3BP Sim        *\n*************************\n");
	test_cr3bp_sim();
	
	printColor(RED, "*************************\n* Test BCR4BPR Sim      *\n*************************\n");
	test_bcr4bpr_sim();

	printColor(RED, "*************************\n* Test CR3BP Events     *\n*************************\n");
	test_cr3bp_events();
	
	printColor(RED, "*************************\n* Test BCR4BPR Events   *\n*************************\n");
	test_bcr4bpr_events();

	return 0;
}//====================================================
/**
 *	Test the simulation engine
 */

#include "tpat_ascii_output.hpp"
#include "tpat_bcr4bpr_sys_data.hpp"
#include "tpat_bcr4bpr_traj.hpp"
#include "tpat_constants.hpp"
#include "tpat_cr3bp_sys_data.hpp"
#include "tpat_cr3bp_traj.hpp"
#include "tpat_matrix.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_utilities.hpp"

#include <iostream>

using namespace std;

void test_cr3bp_sim(){
	tpat_simulation_engine simEngine;

	tpat_cr3bp_sys_data sys("earth", "moon");

	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	simEngine.setSysData(&sys);
	simEngine.setVerbose(false);
	simEngine.setAbsTol(1e-14);
	simEngine.setRelTol(1e-16);
	simEngine.runSim(ic, 2.77);
	tpat_cr3bp_traj traj = simEngine.getCR3BPTraj();

	cout << "Trajectory contains " << traj.getLength() << " points" << endl;
	
	tpat_matrix lastSTM = traj.getSTM(traj.getLength()-1);
	lastSTM.print("%12.4f");

	traj.saveToMat("HaloTest.mat");
}

void test_bcr4bpr_sim(){
	// Do a simulation in the BCR4BP
	tpat_bcr4bpr_sys_data bcSys("sun", "earth", "moon");
	tpat_simulation_engine bcEngine( &bcSys );

	double haloCross177_IC[] = {100.0359099212, 	0, 					285.85225655914e-05, 
								0.0405130514527453, 0.114452371350509, 	0.0310985198400586};
	double t0 = 2.5740750742085;
	
	bcEngine.setVarStepSize(false);
	bcEngine.setNumSteps(500);
	bcEngine.setVerbose(false);
	bcEngine.runSim(haloCross177_IC, t0, 2*PI);

	tpat_bcr4bpr_traj bcTraj = bcEngine.getBCR4BPRTraj();
	cout << "Trajectory contains " << bcTraj.getLength() << " points" << endl;
	bcTraj.saveToMat("bcHaloManifoldProp.mat");
}

void test_cr3bp_events(){
	tpat_cr3bp_sys_data sys("earth", "moon");
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};	// L1 Halo

	tpat_simulation_engine engine(&sys);
	engine.setVerbose(true);
	engine.addEvent(tpat_event::XZ_PLANE, 0, true);
	engine.setRevTime(true);
	engine.runSim(ic, 2.7);

	tpat_cr3bp_traj traj = engine.getCR3BPTraj();
	traj.saveToMat("HaloHalfTest.mat");
}

void test_bcr4bpr_events(){
	tpat_bcr4bpr_sys_data sys("sun", "earth", "moon");
	double ic[] = {100.0359099212, 	0, 					285.85225655914e-05, 
					0.0405130514527453, 0, 	0.0310985198400586};
	double t0 = 2.57;

	tpat_simulation_engine engine(&sys);
	engine.setVerbose(true);
	engine.addEvent(tpat_event::XY_PLANE, 0, true);
	engine.runSim(ic, t0, 2*PI);

	tpat_bcr4bpr_traj traj = engine.getBCR4BPRTraj();
	traj.saveToMat("BC_HaloManifold.mat");
}

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
}
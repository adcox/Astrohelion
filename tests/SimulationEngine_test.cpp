/**
 *	Test the simulation engine
 */


#include "adtk_bcr4bpr_sys_data.hpp"
#include "adtk_bcr4bpr_traj.hpp"
#include "adtk_constants.hpp"
#include "adtk_cr3bp_sys_data.hpp"
#include "adtk_cr3bp_traj.hpp"
#include "adtk_matrix.hpp"
#include "adtk_simulation_engine.hpp"
#include "adtk_utilities.hpp"

#include <iostream>

using namespace std;

void test_cr3bp_sim(){
	adtk_simulation_engine simEngine;

	adtk_cr3bp_sys_data sys("earth", "moon");

	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	simEngine.setSysData(&sys);
	simEngine.setVerbose(true);
	simEngine.setAbsTol(1e-14);
	simEngine.setRelTol(1e-16);
	simEngine.runSim(ic, 2.77);
	adtk_cr3bp_traj traj = simEngine.getCR3BPTraj();

	cout << "Trajectory contains " << traj.getLength() << " points" << endl;
	
	adtk_matrix lastSTM = traj.getSTM(traj.getLength()-1);
	lastSTM.print("%12.4f");

	traj.saveToMat("HaloTest.mat");
}

void test_bcr4bpr_sim(){
	// Do a simulation in the BCR4BP
	adtk_bcr4bpr_sys_data bcSys("sun", "earth", "moon");
	adtk_simulation_engine bcEngine( &bcSys );

	double haloCross177_IC[] = {100.0359099212, 	0, 					285.85225655914e-05, 
								0.0405130514527453, 0.114452371350509, 	0.0310985198400586};
	double t0 = 2.5740750742085;
	
	bcEngine.setVarStepSize(false);
	bcEngine.setNumSteps(500);
	bcEngine.setVerbose(true);
	bcEngine.runSim(haloCross177_IC, t0, 2*PI);

	adtk_bcr4bpr_traj bcTraj = bcEngine.getBCR4BPRTraj();
	cout << "Trajectory contains " << bcTraj.getLength() << " points" << endl;
	bcTraj.saveToMat("bcHaloManifoldProp.mat");
}

void test_cr3bp_events(){
	adtk_cr3bp_sys_data sys("earth", "moon");
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};	// L1 Halo

	adtk_simulation_engine engine(&sys);
	engine.setVerbose(true);
	engine.addEvent(adtk_event::XZ_PLANE, 0, true);
	engine.runSim(ic, 2.77);

	adtk_cr3bp_traj traj = engine.getCR3BPTraj();
	traj.saveToMat("HaloHalfTest.mat");
}

void test_bcr4bpr_events(){
	adtk_bcr4bpr_sys_data sys("sun", "earth", "moon");
	double ic[] = {100.0359099212, 	0, 					285.85225655914e-05, 
					0.0405130514527453, 0.1, 	0.0310985198400586};
	double t0 = 2.57;

	adtk_simulation_engine engine(&sys);
	engine.setVerbose(true);
	engine.addEvent(adtk_event::XY_PLANE, 0, true);
	engine.runSim(ic, t0, 6*PI);

	adtk_bcr4bpr_traj traj = engine.getBCR4BPRTraj();
	traj.saveToMat("BC_HaloManifold.mat");
}
int main(void){
	
	test_cr3bp_sim();
	test_bcr4bpr_sim();

	test_cr3bp_events();
	test_bcr4bpr_events();

	return 0;
}
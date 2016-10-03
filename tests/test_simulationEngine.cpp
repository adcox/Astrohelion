/**
 *	Test the simulation engine
 */

#include "AsciiOutput.hpp"
#include "SysData_bc4bp.hpp"
#include "Traj_bc4bp.hpp"
#include "Common.hpp"
#include "SysData_cr3bp.hpp"
#include "Traj_cr3bp.hpp"
#include "SimEngine.hpp"
#include "Utilities.hpp"

#include <iostream>

using namespace std;
using namespace astrohelion;

void test_cr3bp_sim(){
	SysData_cr3bp sys("earth", "moon");
	SimEngine simEngine;
	Traj_cr3bp traj(&sys);

	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	simEngine.setVerbosity(Verbosity_tp::NO_MSG);
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
	SysData_bc4bp bcSys("sun", "earth", "moon");
	SimEngine bcEngine;
	Traj_bc4bp bcTraj(&bcSys);

	double haloCross177_IC[] = {0.0359099212, 		0, 					285.85225655914e-05, 
								0.0405130514527453, 0.114452371350509, 	0.0310985198400586};
	double t0 = 2.57;
	
	bcEngine.setVarStepSize(false);
	bcEngine.setNumSteps(500);
	bcEngine.setVerbosity(Verbosity_tp::ALL_MSG);
	bcEngine.runSim(haloCross177_IC, t0, 2*PI, &bcTraj);

	cout << "Trajectory contains " << bcTraj.getNumNodes() << " points" << endl;
	bcTraj.saveToMat("data/bcHaloManifoldProp.mat");
}//====================================================

void test_cr3bp_events(){
	SysData_cr3bp sys("earth", "moon");
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};	// L1 Halo

	SimEngine engine;
	Traj_cr3bp traj(&sys);
	engine.setVerbosity(Verbosity_tp::ALL_MSG);
	engine.addEvent(Event(Event_tp::XZ_PLANE, 0, true));
	engine.setRevTime(true);
	engine.runSim(ic, 2.7, &traj);

	traj.saveToMat("data/HaloHalfTest.mat");
}//====================================================

void test_bcr4bpr_events(){
	SysData_bc4bp sys("sun", "earth", "moon");
	double ic[] = {0.0359099212, 		0, 	285.85225655914e-05, 
					0.0405130514527453, 0, 	0.0310985198400586};
	double t0 = 2.57;

	SimEngine engine;
	Traj_bc4bp traj(&sys);
	engine.setVerbosity(Verbosity_tp::ALL_MSG);
	engine.addEvent(Event(Event_tp::XY_PLANE, 0, true));
	engine.runSim(ic, t0, 2*PI, &traj);
	traj.saveToMat("data/BC_HaloManifold.mat");
}//====================================================

int main(void){
	astrohelion::printColor(RED, "*************************\n* Test CR3BP Sim        *\n*************************\n");
	test_cr3bp_sim();
	
	astrohelion::printColor(RED, "*************************\n* Test BCR4BPR Sim      *\n*************************\n");
	test_bcr4bpr_sim();

	astrohelion::printColor(RED, "*************************\n* Test CR3BP Events     *\n*************************\n");
	test_cr3bp_events();
	
	astrohelion::printColor(RED, "*************************\n* Test BCR4BPR Events   *\n*************************\n");
	test_bcr4bpr_events();

	return 0;
}//====================================================
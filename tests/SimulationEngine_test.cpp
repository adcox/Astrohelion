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

#include <iostream>

using namespace std;

int main(void){
	adtk_simulation_engine simEngine;

	adtk_cr3bp_sys_data sys("earth", "moon");

	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	simEngine.setSysData(&sys);
	simEngine.setVerbose(true);
	simEngine.setAbsTol(1e-14);
	simEngine.setRelTol(1e-16);
	simEngine.runSim(ic, 2.77);
	cout << "Extracting trajectory" << endl;
	adtk_cr3bp_traj traj = simEngine.getCR3BPTraj();

	cout << "Trajectory contains " << traj.getLength() << " points" << endl;
	cout << "Final STM is:" << endl;
	
	adtk_matrix lastSTM = traj.getSTM(traj.getLength()-1);
	lastSTM.print("%12.4f");

	traj.saveToMat("HaloTest.mat");


	// Do a simulation in the BCR4BP
	adtk_bcr4bpr_sys_data bcSys("sun", "earth", "moon");
	cout << "System type: " << bcSys.getTypeStr() << endl;
	
	adtk_simulation_engine bcEngine( &bcSys );

	double haloCross177_IC[] = {100.0359099212, 	0, 					285.85225655914e-05, 
								0.0405130514527453, 0.114452371350509, 	0.0310985198400586};
	double t0 = 2.5740750742085;
	// double t0 = 0;

	bcEngine.runSim(haloCross177_IC, t0, 2*PI);

	adtk_bcr4bpr_traj bcTraj = simEngine.getBCR4BPRTraj();
	bcTraj.saveToMat("bcHaloManifoldProp.mat");

	return 0;
}
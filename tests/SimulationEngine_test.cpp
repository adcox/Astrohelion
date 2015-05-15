/**
 *	Test the simulation engine
 */

#include "adtk_simulation_engine.hpp"
#include "adtk_cr3bp_sys_data.hpp"

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

	return 0;
}
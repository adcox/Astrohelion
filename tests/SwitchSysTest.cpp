/**
 *	Test system conversion functions: SE2EM and EM2SE
 */

#include "adtk_bcr4bpr_traj.hpp"
#include "adtk_bcr4bpr_sys_data.hpp"
#include "adtk_calculations.hpp"
#include "adtk_correction_engine.hpp"
#include "adtk_cr3bp_nodeset.hpp"
#include "adtk_cr3bp_sys_data.hpp"
#include "adtk_cr3bp_traj.hpp"
#include "adtk_simulation_engine.hpp"

#include <iostream>

int main(void){
	// Define system data objects
	adtk_cr3bp_sys_data emSys("earth", "moon");
	adtk_cr3bp_sys_data seSys("sun", "earth");

	// Try converting Earth-Moon to Sun-Earth
	double haloIC[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	adtk_simulation_engine engine(&emSys);
	engine.runSim(haloIC, 2.77);
	adtk_cr3bp_traj emHalo = engine.getCR3BPTraj();
	emHalo.saveToMat("EM_Halo.mat");

	adtk_cr3bp_nodeset emNodes(haloIC, emSys, 2.77, 10, adtk_nodeset::TIME);
	emNodes.saveToMat("EM_Nodes.mat");

	adtk_cr3bp_traj emHalo_inSE = cr3bp_EM2SE(emHalo, 0.1, 0.2, 0.05);
	emHalo_inSE.saveToMat("EM_Halo_inSE.mat");

	adtk_cr3bp_nodeset emNodes_inSE = cr3bp_EM2SE(emNodes, 0, 0.1, 0.2, 0.05);
	emNodes_inSE.saveToMat("EM_Nodes_inSE.mat");

	// Try converting Sun-Earth to Earth-Moon
	engine.reset();
	engine.setSysData(&seSys);
	engine.runSim(haloIC, 2.77);
	adtk_cr3bp_traj seTraj = engine.getCR3BPTraj();
	seTraj.saveToMat("SE_Traj.mat");

	adtk_cr3bp_nodeset seNodes(haloIC, seSys, 2.77, 10, adtk_nodeset::ARCLENGTH);
	seNodes.saveToMat("SE_Nodes.mat");

	adtk_cr3bp_traj seTraj_inEM = cr3bp_SE2EM(seTraj, 0.1, 0.2, 0.05);
	seTraj_inEM.saveToMat("SE_Traj_inEM.mat");

	adtk_cr3bp_nodeset seNodes_inEM = cr3bp_SE2EM(seNodes, 0, 0.1, 0.2, 0.05);
	seNodes_inEM.saveToMat("SE_Nodes_inEM.mat");

	// Try converting Sun-Earth to Sun-Earth-Moon
	adtk_bcr4bpr_sys_data bcSys("sun", "earth", "moon");
	adtk_bcr4bpr_traj bcTraj = bcr4bpr_SE2SEM(seTraj, bcSys, 7.08);
	bcTraj.saveToMat("SEM_Traj.mat");

	return 0;
}
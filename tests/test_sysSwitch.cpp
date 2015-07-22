/**
 *	Test system conversion functions: SE2EM and EM2SE
 */

#include "tpat_nodeset_bcr4bpr.hpp"
#include "tpat_traj_bcr4bpr.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_calculations.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"

#include <iostream>

int main(void){
	// Define system data objects
	tpat_sys_data_cr3bp emSys("earth", "moon");
	tpat_sys_data_cr3bp seSys("sun", "earth");

	// Try converting Earth-Moon to Sun-Earth
	double haloIC[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	tpat_simulation_engine engine(emSys);
	engine.runSim(haloIC, 2.77);
	tpat_traj_cr3bp emHalo = engine.getCR3BP_Traj();
	emHalo.saveToMat("EM_Halo.mat");

	tpat_nodeset_cr3bp emNodes(haloIC, emSys, 2.77, 10, tpat_nodeset::DISTRO_TIME);
	emNodes.saveToMat("EM_Nodes.mat");

	tpat_traj_cr3bp emHalo_inSE = cr3bp_EM2SE(emHalo, 0.1, 0.2, 0.05);
	emHalo_inSE.saveToMat("EM_Halo_inSE.mat");

	tpat_nodeset_cr3bp emNodes_inSE = cr3bp_EM2SE(emNodes, 0, 0.1, 0.2, 0.05);
	emNodes_inSE.saveToMat("EM_Nodes_inSE.mat");

	// Try converting Sun-Earth to Earth-Moon
	engine.reset();
	engine.setSysData(seSys);
	engine.runSim(haloIC, 2.77);
	tpat_traj_cr3bp seTraj = engine.getCR3BP_Traj();
	seTraj.saveToMat("SE_Traj.mat");

	tpat_nodeset_cr3bp seNodes(haloIC, seSys, 2.77, 10, tpat_nodeset::DISTRO_ARCLENGTH);
	seNodes.saveToMat("SE_Nodes.mat");

	tpat_traj_cr3bp seTraj_inEM = cr3bp_SE2EM(seTraj, 0.1, 0.2, 0.05);
	seTraj_inEM.saveToMat("SE_Traj_inEM.mat");

	tpat_nodeset_cr3bp seNodes_inEM = cr3bp_SE2EM(seNodes, 0, 0.1, 0.2, 0.05);
	seNodes_inEM.saveToMat("SE_Nodes_inEM.mat");

	// Try converting Sun-Earth to Sun-Earth-Moon
	tpat_sys_data_bcr4bpr bcSys("sun", "earth", "moon");
	tpat_traj_bcr4bpr bcTraj = bcr4bpr_SE2SEM(seTraj, bcSys, 7.08);
	bcTraj.saveToMat("SEM_Traj.mat");

	tpat_nodeset_bcr4bpr bcNodes = bcr4bpr_SE2SEM(seNodes, bcSys, 7.08);
	bcNodes.saveToMat("SEM_Nodes.mat");
	return 0;
}
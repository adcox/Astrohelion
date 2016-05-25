/**
 *	Test system conversion functions: SE2EM and EM2SE
 */

#include "tpat_nodeset_bc4bp.hpp"
#include "tpat_traj_bc4bp.hpp"
#include "tpat_sys_data_bc4bp.hpp"
#include "tpat_calculations.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_sim_engine.hpp"

#include <iostream>

int main(void){
	// Define system data objects
	TPAT_Sys_Data_CR3BP emSys("earth", "moon");
	TPAT_Sys_Data_CR3BP seSys("sun", "earth");

	// Try converting Earth-Moon to Sun-Earth
	double haloIC[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	TPAT_Sim_Engine engine;
	TPAT_Traj_CR3BP emHalo(&emSys);
	engine.runSim(haloIC, 2.77, &emHalo);
	emHalo.saveToMat("data/EM_Halo.mat");

	TPAT_Nodeset_CR3BP emNodes(haloIC, &emSys, 2.77, 10, TPAT_Nodeset::DISTRO_TIME);
	emNodes.saveToMat("data/EM_Nodes.mat");

	TPAT_Traj_CR3BP emHalo_inSE = cr3bp_EM2SE(emHalo, &seSys, 0.1, 0.2, 0.05);
	emHalo_inSE.saveToMat("data/EM_Halo_inSE.mat");

	TPAT_Nodeset_CR3BP emNodes_inSE = cr3bp_EM2SE(emNodes, &seSys, 0, 0.1, 0.2, 0.05);
	emNodes_inSE.saveToMat("data/EM_Nodes_inSE.mat");

	// Try converting Sun-Earth to Earth-Moon
	engine.reset();
	TPAT_Traj_CR3BP seTraj(&seSys);
	engine.runSim(haloIC, 2.77, &seTraj);
	seTraj.saveToMat("data/SE_Traj.mat");

	TPAT_Nodeset_CR3BP seNodes(haloIC, &seSys, 2.77, 10, TPAT_Nodeset::DISTRO_ARCLENGTH);
	seNodes.saveToMat("data/SE_Nodes.mat");

	TPAT_Traj_CR3BP seTraj_inEM = cr3bp_SE2EM(seTraj, &emSys, 0.1, 0.2, 0.05);
	seTraj_inEM.saveToMat("data/SE_Traj_inEM.mat");

	TPAT_Nodeset_CR3BP seNodes_inEM = cr3bp_SE2EM(seNodes, &emSys, 0, 0.1, 0.2, 0.05);
	seNodes_inEM.saveToMat("data/SE_Nodes_inEM.mat");

	// Try converting Sun-Earth to Sun-Earth-Moon
	TPAT_Sys_Data_BC4BP bcSys("sun", "earth", "moon");
	TPAT_Traj_BC4BP bcTraj = bcr4bpr_SE2SEM(seTraj, &bcSys, 7.08);
	bcTraj.saveToMat("data/SEM_Traj.mat");

	TPAT_Nodeset_BC4BP bcNodes = bcr4bpr_SE2SEM(seNodes, &bcSys, 7.08);
	bcNodes.saveToMat("data/SEM_Nodes.mat");
	return 0;
}
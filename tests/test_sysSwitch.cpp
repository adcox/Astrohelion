/**
 *	Test system conversion functions: SE2EM and EM2SE
 */

#include "Nodeset_bc4bp.hpp"
#include "Traj_bc4bp.hpp"
#include "SysData_bc4bp.hpp"
#include "Calculations.hpp"
#include "CorrectionEngine.hpp"
#include "Nodeset_cr3bp.hpp"
#include "SysData_cr3bp.hpp"
#include "Traj_cr3bp.hpp"
#include "SimEngine.hpp"

#include <iostream>

using namespace astrohelion;

int main(void){
	// Define system data objects
	SysData_cr3bp emSys("earth", "moon");
	SysData_cr3bp seSys("sun", "earth");

	// Try converting Earth-Moon to Sun-Earth
	double haloIC[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	SimEngine engine;
	Traj_cr3bp emHalo(&emSys);
	engine.runSim(haloIC, 2.77, &emHalo);
	emHalo.saveToMat("data/EM_Halo.mat");

	Nodeset_cr3bp emNodes(haloIC, &emSys, 2.77, 10, Nodeset::DISTRO_TIME);
	emNodes.saveToMat("data/EM_Nodes.mat");

	Traj_cr3bp emHalo_inSE = cr3bp_EM2SE(emHalo, &seSys, 0.1, 0.2, 0.05);
	emHalo_inSE.saveToMat("data/EM_Halo_inSE.mat");

	Nodeset_cr3bp emNodes_inSE = cr3bp_EM2SE(emNodes, &seSys, 0.1, 0.2, 0.05);
	emNodes_inSE.saveToMat("data/EM_Nodes_inSE.mat");


	double haloIC_SE[] = {1.00211887846215, 0.000211363695548905, 0.000101518083014945, -0.00327841590034548, 0.0327449242150319, 0.000347950888417211};
	double halo_final_SE[] = {0.998125006225457, 0.000978215985711854, 0.000310194273165479, -0.0151681647612735, -0.0292266231714366, 0.000365438430426524}
	
	// Try converting Sun-Earth to Earth-Moon
	engine.reset();
	Traj_cr3bp seTraj(&seSys);
	engine.runSim(haloIC, 2.77, &seTraj);
	seTraj.saveToMat("data/SE_Traj.mat");

	Nodeset_cr3bp seNodes(haloIC, &seSys, 2.77, 10, Nodeset::DISTRO_ARCLENGTH);
	seNodes.saveToMat("data/SE_Nodes.mat");

	Traj_cr3bp seTraj_inEM = cr3bp_SE2EM(seTraj, &emSys, 0.1, 0.2, 0.05);
	seTraj_inEM.saveToMat("data/SE_Traj_inEM.mat");

	Nodeset_cr3bp seNodes_inEM = cr3bp_SE2EM(seNodes, &emSys, 0.1, 0.2, 0.05);
	seNodes_inEM.saveToMat("data/SE_Nodes_inEM.mat");

	// Try converting Sun-Earth to Sun-Earth-Moon
	SysData_bc4bp bcSys("sun", "earth", "moon");
	Traj_bc4bp bcTraj = bcr4bpr_SE2SEM(seTraj, &bcSys, 0, 7.08);
	bcTraj.saveToMat("data/SEM_Traj.mat");

	Nodeset_bc4bp bcNodes = bcr4bpr_SE2SEM(seNodes, &bcSys, 0, 7.08);
	bcNodes.saveToMat("data/SEM_Nodes.mat");


	return 0;
}
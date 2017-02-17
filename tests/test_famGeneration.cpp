
#include "Fam_cr3bp.hpp"
#include "SysData_cr3bp.hpp"
#include "Traj_cr3bp.hpp"
#include "FamGenerator.hpp"
#include "Utilities.hpp"

#include <cmath>
#include <iostream>
#include <vector>

using namespace astrohelion;

int main(int argc, char *argv[]){

	(void) argc;
	(void) argv;

	// Create system data and generator family generator
	// SysData_cr3bp sysData("sun", "earth");
	SysData_cr3bp sysData("../share/SE_DPO_fromHill.mat");
	// SysData_cr3bp sysData("earth", "moon");
	Fam_cr3bp DPO(sysData);

	FamGenerator gen;
	gen.setNumOrbits(4);
	gen.setContType(Continuation_tp::PSEUDO_ARC);
	gen.setNumNodes(11);
	gen.setTol(1e-11);

	// gen.cr3bp_generateDPO(&DPO);
	Traj_cr3bp dpoTraj("../share/SE_DPO_fromHill.mat", &sysData);
	std::vector<int> dir {0,0,0,0,1,0};	// Step in + ydot first
	
	// Use pseudo arclength with the precomputed DPO
	gen.cr3bp_pacFromTraj(dpoTraj, Mirror_tp::MIRROR_XZ, dir, &DPO);

	// Run the other direction too
	dir[4] *= -1;
	gen.cr3bp_pacFromTraj(dpoTraj, Mirror_tp::MIRROR_XZ, dir, &DPO);
		
	DPO.sortEigs();
	DPO.setName("Sun-Earth Distant Prograde Orbit");
	DPO.saveToMat("../share/families/SE_DPO_PAC.mat");


	// printf("Generating Sun-Earth L1 Northern Halo Family : PAC\n");
	// gen.setContType(Continuation_tp::PSEUDO_ARC);
	// gen.setTol(6e-12);
	// gen.setNumNodes(12);
	
	// Fam_cr3bp L1_Halo = gen.cr3bp_generateHalo("../share/families_natParam_checked/SE_L1_Lyap.mat", -1e-4);
	// L1_Halo.sortEigs();
	// L1_Halo.setName("Sun-Earth L1 Northern Halo");
	// L1_Halo.saveToMat("../share/families/SE_L1_NHalo.mat");

	return EXIT_SUCCESS;
}
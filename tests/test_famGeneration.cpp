
#include "Fam_cr3bp.hpp"
#include "SysData_cr3bp.hpp"
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
	SysData_cr3bp sysData("earth", "moon");
	FamGenerator gen;
	gen.setContType(Continuation_tp::PSEUDO_ARC);

	Fam_cr3bp DPO = gen.cr3bp_generateDPO(&sysData);
	DPO.sortEigs();
	DPO.setName("Earth-Moon Distant Prograde Orbit");
	DPO.saveToMat("../share/families/EM_DPO.mat");

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
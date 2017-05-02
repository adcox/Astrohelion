
#include "Fam_cr3bp.hpp"
#include "SysData_cr3bp.hpp"
#include "Arcset_cr3bp.hpp"
#include "FamGenerator.hpp"
#include "Utilities.hpp"

#include <cmath>
#include <iostream>
#include <vector>

using namespace astrohelion;

void createEMLyap(){
	SysData_cr3bp sys("earth", "moon");
	Fam_cr3bp fam(sys);

	// Natural Parameter Continuation
	FamGenerator gen;
	gen.setStep_simple(0.0005);
	gen.setStep_fitted_1(0.005);
	gen.setStep_fitted_2(0.005);
	gen.setContType(Continuation_tp::NAT_PARAM);
	// gen.setNumOrbits(5);

	// Run the generator
	gen.cr3bp_generateLyap(1, 0.001, &fam);
	fam.sortEigs();
	fam.setName("Earth-Moon L1 Lyapunov");
	fam.saveToMat("../share/families/EM_L1_Lyap.mat");
}

void createSEDPO(){
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
	Arcset_cr3bp dpoTraj("../share/SE_DPO_fromHill.mat", &sysData);
	std::vector<int> dir {0,0,0,0,1,0};	// Step in + ydot first
	
	// Use pseudo arclength with the precomputed DPO
	gen.cr3bp_pacFromArcset(dpoTraj, Mirror_tp::MIRROR_XZ, dir, &DPO);

	// Run the other direction too
	dir[4] *= -1;
	gen.cr3bp_pacFromArcset(dpoTraj, Mirror_tp::MIRROR_XZ, dir, &DPO);
		
	DPO.sortEigs();
	DPO.setName("Sun-Earth Distant Prograde Orbit");
	DPO.saveToMat("../share/families/SE_DPO_PAC.mat");
}


int main(int argc, char *argv[]){

	(void) argc;
	(void) argv;

	createEMLyap();	
	// createSEDPO();


	return EXIT_SUCCESS;
}
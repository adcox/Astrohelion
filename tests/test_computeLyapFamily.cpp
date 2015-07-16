
#include "tpat_cr3bp_family.hpp"
#include "tpat_cr3bp_sys_data.hpp"
#include "tpat_family_generator.hpp"

#include <cmath>
#include <iostream>
#include <vector>

/**
 *	Compute the three Lyapunov families in the Earth-Moon system
 */
int main(void){
	tpat_cr3bp_sys_data sysData("earth", "Moon");
	tpat_family_generator gen;
	gen.setStep_fitted(0.0005);
	tpat_cr3bp_family L1_Lyap = gen.cr3bp_generateLyap(sysData, 1, 0.001);
	L1_Lyap.setName("Earth-Moon L1 Lyapunov");
	L1_Lyap.saveToMat("../share/families/EM_L1_Lyap.mat");

	gen.setStep_fitted(0.001);
	tpat_cr3bp_family L2_Lyap = gen.cr3bp_generateLyap(sysData, 2, 0.01);
	L2_Lyap.setName("Earth-Moon L2 Lyapunov");
	L2_Lyap.saveToMat("../share/families/EM_L2_Lyap.mat");

	gen.setNumNodes(5);
	gen.setStep_fitted(0.005);
	tpat_cr3bp_family L3_Lyap = gen.cr3bp_generateLyap(sysData, 3, 0.01);
	L3_Lyap.setName("Earth-Moon L3 Lyapunov");
	L3_Lyap.saveToMat("../share/families/EM_L3_Lyap.mat");
	
}
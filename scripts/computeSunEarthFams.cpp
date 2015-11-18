
#include "tpat_family_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_family_generator.hpp"
#include "tpat_utilities.hpp"

#include <cmath>
#include <iostream>
#include <vector>

int main(int argc, char *argv[]){
	int probNum = 0;
	if(argc > 1)
		probNum = atoi(argv[1]);
	else{
		printf("Please enter an integer to specify which Sun-Earth family to generate\n");
		return EXIT_SUCCESS;
	}
	// Create system data and generator family generator
	tpat_sys_data_cr3bp sysData("sun", "earth");
	tpat_family_generator gen;

	switch(probNum){
		case 0:
		{
			printf("Generating Sun-Earth L1 Lyapunov Family : Nat Param\n");
			gen.setContType(tpat_family_generator::NAT_PARAM);
			gen.setStep_simple(0.00001);
			gen.setStep_fitted_1(0.00007);
			gen.setStep_fitted_2(0.00007);
			gen.setNumOrbits(1500);

			// Run the generator
			tpat_family_cr3bp L1_Lyap = gen.cr3bp_generateLyap(sysData, 1, 0.001);
			L1_Lyap.sortEigs();
			L1_Lyap.setName("Sun-Earth L1 Lyapunov");
			L1_Lyap.saveToMat("../share/families/SE_L1_Lyap.mat");
			break;
		}
		case 100:
		{
			printf("Generating Sun-Earth L1 Lyapunov Family : PAC\n");
			gen.setContType(tpat_family_generator::PSEUDO_ARC);
			printErr("Not yet implemented...\n");
			break;
		}
		case 1:
		{
			printf("Generating Sun-Earth L2 Lyapunov Family : Nat Param\n");
			gen.setContType(tpat_family_generator::NAT_PARAM);
			printErr("Not yet implemented...\n");
			break;
		}
		case 101:
		{
			printf("Generating Sun-Earth L2 Lyapunov Family : PAC\n");
			gen.setContType(tpat_family_generator::PSEUDO_ARC);
			printErr("Not yet implemented...\n");
			break;
		}
		case 2:
		{
			printf("Generating Sun-Earth L3 Lyapunov Family : Nat Param\n");
			gen.setContType(tpat_family_generator::NAT_PARAM);
			printErr("Not yet implemented...\n");
			break;
		}
		case 102:
		{
			printf("Generating Sun-Earth L3 Lyapunov Family : PAC\n");
			gen.setContType(tpat_family_generator::PSEUDO_ARC);
			printErr("Not yet implemented...\n");
			break;
		}
		case 3:
		{
			printf("Generating Sun-Earth L1 Northern Halo Family : Nat Param\n");
			gen.setContType(tpat_family_generator::NAT_PARAM);
			gen.setStep_fitted_1(0.0001);
			gen.setStep_fitted_2(0.0001);
			gen.setStep_simple(0.0001);
			gen.setNumNodes(5);

			tpat_family_cr3bp L1_Halo = gen.cr3bp_generateHalo("../share/families_natParam_checked/SE_L1_Lyap.mat", -1e-4);
			L1_Halo.sortEigs();
			L1_Halo.setName("Sun-Earth L1 Northern Halo");
			L1_Halo.saveToMat("../share/families/SE_L1_NHalo.mat");
			break;
		}
		case 103:
		{
			printf("Generating Sun-Earth L1 Northern Halo Family : PAC\n");
			gen.setContType(tpat_family_generator::PSEUDO_ARC);
			gen.setTol(6e-12);

			tpat_family_cr3bp L1_Halo = gen.cr3bp_generateHalo("../share/families_natParam_checked/SE_L1_Lyap.mat", -1e-4);
			L1_Halo.sortEigs();
			L1_Halo.setName("Sun-Earth L1 Northern Halo");
			L1_Halo.saveToMat("../share/families/SE_L1_NHalo.mat");
			break;
		}
		case 4:
		{
			printf("Generating Sun-Earth L2 Northern Halo Family : Nat Param\n");
			gen.setContType(tpat_family_generator::NAT_PARAM);
			printErr("Not yet implemented...\n");
			break;
		}
		case 104:
		{
			printf("Generating Sun-Earth L2 Northern Halo Family : PAC\n");
			gen.setContType(tpat_family_generator::PSEUDO_ARC);
			printErr("Not yet implemented...\n");
			break;
		}
		case 5:
		{
			printf("Generating Sun-Earth L3 Northern Halo Family : Nat Param\n");
			gen.setContType(tpat_family_generator::NAT_PARAM);
			printErr("Not yet implemented...\n");
			break;
		}
		case 105:
		{
			printf("Generating Sun-Earth L3 Northern Halo Family : PAC\n");
			gen.setContType(tpat_family_generator::PSEUDO_ARC);
			printErr("Not yet implemented...\n");
			break;
		}
		case 6:
		{
			printf("Generating Sun-Earth L1 Northern Axial Family : Nat Param\n");
			gen.setContType(tpat_family_generator::NAT_PARAM);
			gen.setStep_fitted_1(0.0003);
			gen.setStep_fitted_2(0.0003);
			gen.setNumNodes(5);
			gen.setNumOrbits(320);	// family starts repeating after this many...

			tpat_family_cr3bp L1_Axial = gen.cr3bp_generateAxial("../share/families_natParam_checked/SE_L1_Lyap.mat", 1e-4);
			L1_Axial.sortEigs();
			L1_Axial.setName("Sun-Earth L1 Northern Axial");
			L1_Axial.saveToMat("../share/families/SE_L1_NAxial.mat");
			break;
		}
		case 7:
		{
			printf("Generating Sun-Earth L2 Northern Axial Family : Nat Param\n");
			gen.setContType(tpat_family_generator::NAT_PARAM);
			printErr("Not yet implemented...\n");
			break;	
		}
		case 8:
		{
			printf("Generating Sun-Earth L3 Northern Axial Family : Nat Param\n");
			gen.setContType(tpat_family_generator::NAT_PARAM);
			printErr("Not yet implemented...\n");
			break;	
		}
		case 9:
		{
			printf("Generating Sun-Earth L1 Vertical Family : Nat Param\n");
			gen.setContType(tpat_family_generator::NAT_PARAM);
			printErr("Not yet implemented...\n");
			break;
		}
		case 109:
		{
			printf("Generating Sun-Earth L1 Vertical Family : PAC\n");
			gen.setContType(tpat_family_generator::PSEUDO_ARC);
			gen.setNumOrbits(1000);

			// tpat_family_cr3bp L1_Vert = gen.cr3bp_generateVertical("../share/families_natParam_checked/SE_L1_NAxial.mat", 0.0001);
			// L1_Vert.sortEigs();
			// L1_Vert.setName("Sun-Earth L1 Vertical");
			// L1_Vert.saveToMat("../share/families/SE_L1_Vert_PAC.mat");

			// Generate other side of family
			// gen.setNumOrbits(3);
			tpat_family_cr3bp L1_Vert_small = gen.cr3bp_generateVertical("../share/families_natParam_checked/SE_L1_NAxial.mat", -0.0001);
			L1_Vert_small.sortEigs();
			L1_Vert_small.setName("Sun-Earth L1 Vertical");
			L1_Vert_small.saveToMat("../share/families/SE_L1_Vert_Small_PAC.mat");
			break;
		}
		default:
			printf("The problem number %d has not been implemented\n", probNum);
			break;
	}

	return EXIT_SUCCESS;
}
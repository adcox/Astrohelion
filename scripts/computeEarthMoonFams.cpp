
#include "tpat_family_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_family_generator.hpp"
#include "tpat_utilities.hpp"

#include <cmath>
#include <iostream>
#include <vector>

/**
 *	Compute the three Lyapunov families in the Earth-Moon system
 */
int main(int argc, char *argv[]){
	int probNum = 0;
	if(argc > 1)
		probNum = atoi(argv[1]);
	else{
		printf("Please enter an integer to specify which Earth-Moon family to generate\n");
		return EXIT_SUCCESS;
	}

	int p = 0;
	int q = 0;
	if(argc == 4){
		p = atoi(argv[2]);
		q = atoi(argv[3]);
	}

	// Create system data and generator family generator
	tpat_sys_data_cr3bp sysData("earth", "Moon");
	tpat_family_generator gen;

	switch(probNum){
		case 0:
		{
			printf("Generating Earth-Moon L1 Lyapunov Family : Nat Param\n");
			
			// Natural Parameter Continuation
			gen.setStep_simple(0.0005);
			gen.setStep_fitted_1(0.005);
			gen.setStep_fitted_2(0.005);
			gen.setContType(tpat_family_generator::NAT_PARAM);
			// gen.setNumOrbits(5);

			// Run the generator
			tpat_family_cr3bp L1_Lyap = gen.cr3bp_generateLyap(sysData, 1, 0.001);
			L1_Lyap.sortEigs();
			L1_Lyap.setName("Earth-Moon L1 Lyapunov");
			L1_Lyap.saveToMat("../share/families/EM_L1_Lyap.mat");
			break;
		}
		case 100:
		{
			printf("Generating Earth-Moon L1 Lyapunov Family : PAC\n");
			
			// Pseudo Arclength (comment out for Natural Parameter)
			gen.setNumOrbits(1500);
			gen.setContType(tpat_family_generator::PSEUDO_ARC);

			// Run the generator
			tpat_family_cr3bp L1_Lyap = gen.cr3bp_generateLyap(sysData, 1, 0.001);
			L1_Lyap.sortEigs();
			L1_Lyap.setName("Earth-Moon L1 Lyapunov");
			L1_Lyap.saveToMat("../share/families/EM_L1_Lyap_PAC.mat");
			break;
		}
		case 1:
		{
			printf("Generating Earth-Moon L2 Lyapunov Family\n");

			// Natural Parameter Continuation
			gen.setStep_fitted_1(0.005);
			gen.setStep_fitted_2(0.005);
			gen.setContType(tpat_family_generator::NAT_PARAM);

			tpat_family_cr3bp L2_Lyap = gen.cr3bp_generateLyap(sysData, 2, 0.01);
			L2_Lyap.sortEigs();
			L2_Lyap.setName("Earth-Moon L2 Lyapunov");
			L2_Lyap.saveToMat("../share/families/EM_L2_Lyap.mat");
			break;
		}
		case 101:
		{
			printf("Generating Earth-Moon L2 Lyapunov Family : PAC\n");
			
			// Pseudo Arclength (comment out for Natural Parameter)
			gen.setNumOrbits(2000);
			gen.setNumNodes(5);
			gen.setContType(tpat_family_generator::PSEUDO_ARC);

			// Run the generator
			tpat_family_cr3bp L2_Lyap = gen.cr3bp_generateLyap(sysData, 2, 0.01);
			L2_Lyap.sortEigs();
			L2_Lyap.setName("Earth-Moon L2 Lyapunov");
			L2_Lyap.saveToMat("../share/families/EM_L2_Lyap_PAC.mat");
			break;
		}
		case 2:
		{
			printf("Generating Earth-Moon L3 Lyapunov Family\n");

			// Natural Parameter Continuation
			gen.setNumNodes(5);
			gen.setStep_fitted_1(0.05);
			gen.setStep_fitted_2(0.05);
			gen.setNumOrbits(800);
			gen.setContType(tpat_family_generator::NAT_PARAM);

			tpat_family_cr3bp L3_Lyap = gen.cr3bp_generateLyap(sysData, 3, 0.01);
			L3_Lyap.sortEigs();
			L3_Lyap.setName("Earth-Moon L3 Lyapunov");
			L3_Lyap.saveToMat("../share/families/EM_L3_Lyap.mat");
			break;
		}
		case 102:
		{
			printf("Generating Earth-Moon L3 Lyapunov Family : PAC\n");
			
			// Pseudo Arclength (comment out for Natural Parameter)
			gen.setNumOrbits(3000);
			gen.setContType(tpat_family_generator::PSEUDO_ARC);

			// Run the generator
			tpat_family_cr3bp L3_Lyap = gen.cr3bp_generateLyap(sysData, 3, 0.01);
			L3_Lyap.sortEigs();
			L3_Lyap.setName("Earth-Moon L3 Lyapunov");
			L3_Lyap.saveToMat("../share/families/EM_L3_Lyap_PAC.mat");
			break;
		}
		case 3:
		{
			printf("Generating Earth-Moon L1 Northern Halo Family\n");
			gen.setStep_fitted_1(0.001);
			gen.setStep_fitted_2(0.001);
			
			tpat_family_cr3bp L1_Halo = gen.cr3bp_generateHalo("../share/families/EM_L1_Lyap.mat", -5.2e-5);
			L1_Halo.sortEigs();
			L1_Halo.setName("Earth-Moon L1 Northern Halo");
			L1_Halo.saveToMat("../share/families/EM_L1_NHalo.mat");
			break;
		}
		case 103:
		{
			printf("Generating Earth-Moon L1 Northern Halo Family : PAC\n");
			
			// Pseudo Arclength (comment out for Natural Parameter)
			gen.setNumOrbits(3000);
			gen.setContType(tpat_family_generator::PSEUDO_ARC);
			gen.setTol(6e-12);

			// Run the generator
			tpat_family_cr3bp L1_Halo = gen.cr3bp_generateHalo("../share/families_natParam_checked/EM_L1_Lyap.mat", -5.2e-5);
			L1_Halo.sortEigs();
			L1_Halo.setName("Earth-Moon L1 Northern Halo");
			L1_Halo.saveToMat("../share/families/EM_L1_NHalo_PAC.mat");
			break;
		}
		case 4:
		{
			printf("Generating Earth-Moon L2 Northern Halo Family\n");
			gen.setStep_fitted_1(0.001);
			gen.setStep_fitted_2(0.001);
			tpat_family_cr3bp L2_Halo = gen.cr3bp_generateHalo("../share/families/EM_L2_Lyap.mat", 5.2e-5);
			L2_Halo.sortEigs();
			L2_Halo.setName("Earth-Moon L2 Northern Halo");
			L2_Halo.saveToMat("../share/families/EM_L2_NHalo.mat");
			break;
		}
		case 104:
		{
			printf("Generating Earth-Moon L2 Northern Halo Family : PAC\n");
			
			// Pseudo Arclength (comment out for Natural Parameter)
			gen.setNumOrbits(3000);
			gen.setContType(tpat_family_generator::PSEUDO_ARC);
			gen.setTol(6e-12);

			// Run the generator
			tpat_family_cr3bp L2_Halo = gen.cr3bp_generateHalo("../share/families_natParam_checked/EM_L2_Lyap.mat", 5.2e-5);
			L2_Halo.sortEigs();
			L2_Halo.setName("Earth-Moon L2 Northern Halo");
			L2_Halo.saveToMat("../share/families/EM_L2_NHalo_PAC.mat");
			break;
		}
		case 5:
		{
			printf("Generating Earth-Moon L3 Northern Halo Family\n");
			gen.setStep_fitted_1(0.001);
			gen.setStep_fitted_2(0.001);
			tpat_family_cr3bp L3_Halo = gen.cr3bp_generateHalo("../share/families/EM_L3_Lyap.mat", -5.2e-5);
			L3_Halo.sortEigs();
			L3_Halo.setName("Earth-Moon L3 Northern Halo");
			L3_Halo.saveToMat("../share/families/EM_L3_NHalo.mat");
			break;
		}
		case 6:
		{
			printf("Generating Earth-Moon L1 Northern Axial Family\n");
			gen.setNumNodes(5);
			gen.setNumOrbits(285);	// family starts repeating after this many...
			tpat_family_cr3bp L1_Axial = gen.cr3bp_generateAxial("../share/families/EM_L1_Lyap.mat", 1e-4);
			L1_Axial.sortEigs();
			L1_Axial.setName("Earth-Moon L1 Northern Axial");
			L1_Axial.saveToMat("../share/families/EM_L1_NAxial.mat");
			break;
		}
		case 7:
		{
			printf("Generating Earth-Moon L2 Northern Axial Family\n");
			gen.setNumNodes(5);
			gen.setSlopeThresh(1.3);
			gen.setNumOrbits(316);	// family starts repeating after this many...
			tpat_family_cr3bp L2_Axial = gen.cr3bp_generateAxial("../share/families/EM_L2_Lyap.mat", 1e-4);
			L2_Axial.sortEigs();
			L2_Axial.setName("Earth-Moon L2 Northern Axial");
			L2_Axial.saveToMat("../share/families/EM_L2_NAxial.mat");
			break;	
		}
		case 8:
		{
			printf("Generating Earth-Moon L3 Northern Axial Family\n");
			gen.setNumNodes(5);
			gen.setStep_fitted_1(0.01);
			gen.setStep_fitted_2(0.01);
			gen.setSlopeThresh(3);
			gen.setNumOrbits(1000);
			// gen.setNumOrbits(316);	// family starts repeating after this many...
			tpat_family_cr3bp L3_Axial = gen.cr3bp_generateAxial("../share/families/EM_L3_Lyap.mat", 1e-2);
			L3_Axial.sortEigs();
			L3_Axial.setName("Earth-Moon L3 Northern Axial");
			L3_Axial.saveToMat("../share/families/EM_L3_NAxial.mat");
			break;	
		}
		case 9:
		{
			printf("Generating Earth-Moon L1 Vertical Family via Natural Parameter\n");
			gen.setNumNodes(3);
			gen.setStep_fitted_1(0.001);
			gen.setStep_fitted_2(0.001);
			gen.setContType(tpat_family_generator::NAT_PARAM);

			tpat_family_cr3bp L1_Vert = gen.cr3bp_generateVertical("../share/families_natParam_checked/EM_L1_NAxial.mat", -0.0005);
			L1_Vert.sortEigs();
			L1_Vert.setName("Earth-Moon L1 Vertical");
			L1_Vert.saveToMat("../share/families/EM_L1_Vert.mat");
			break;
		}
		case 109:
		{
			printf("Generating Earth-Moon L1 Vertical Family via Pseudo Arc-Length\n");
			gen.setContType(tpat_family_generator::PSEUDO_ARC);
			gen.setNumOrbits(1000);

			tpat_family_cr3bp L1_Vert = gen.cr3bp_generateVertical("../share/families_natParam_checked/EM_L1_NAxial.mat", 0.001);
			L1_Vert.sortEigs();
			L1_Vert.setName("Earth-Moon L1 Vertical");
			L1_Vert.saveToMat("../share/families/EM_L1_Vert_PAC.mat");
			break;
		}
		case 10:
		{
			printf("Generating Earth-Moon DRO Family via Natural Parameter\n");
			gen.setNumNodes(4);
			gen.setStep_simple(0.001);
			gen.setStep_fitted_1(0.01);
			gen.setStep_fitted_2(0.01);
			gen.setMaxStepSize(0.01);
			gen.setNumOrbits(1000);
			gen.setContType(tpat_family_generator::NAT_PARAM);

			tpat_family_cr3bp DRO = gen.cr3bp_generateDRO(&sysData);
			DRO.sortEigs();
			DRO.setName("Earth-Moon DRO");
			DRO.saveToMat("../share/families/EM_DRO.mat");
			break;
		}
		case 110:
		{
			printf("Generating Earth-Moon DRO Family via Pseudo Arclength\n");
			printWarn("Not Yet Implemented");
			break;	
		}
		case 11:
		{
			printf("Generating Earth-Moon LPO Family via Natural Parameter\n");
			gen.setNumNodes(4);
			gen.setStep_simple(0.001);
			gen.setStep_fitted_1(0.01);
			gen.setStep_fitted_2(0.01);
			gen.setMaxStepSize(0.01);
			gen.setNumOrbits(1000);
			gen.setContType(tpat_family_generator::NAT_PARAM);

			tpat_family_cr3bp LPO = gen.cr3bp_generateLPO(&sysData);
			LPO.sortEigs();
			LPO.setName("Earth-Moon LPO");
			LPO.saveToMat("../share/families/EM_LPO.mat");
			break;
		}
		case 111:
		{
			printf("Generating Earth-Moon LPO Family via Pseudo Arclength\n");
			printWarn("Not Yet Implemented");
			break;	
		}
		case 12:
		{
			printf("Generating Earth-Moon DPO Family via Natural Parameter\n");
			printWarn("Not Yet Implemented");
			break;
		}
		case 112:
		{
			printf("Generating Earth-Moon DPO Family via Pseudo Arclength\n");
			printWarn("Not Yet Implemented");
			break;	
		}
		case 20:
		{
			printf("Generating Earth-Moon Resonant Orbit via Natural Parameter\n");
			gen.setContType(tpat_family_generator::NAT_PARAM);
			gen.setNumOrbits(1000);

			tpat_family_cr3bp res = gen.cr3bp_generateRes(&sysData, p, q);
			res.sortMembers();
			res.sortEigs();

			char famName[32], filename[128];
			sprintf(famName, "Earth-Moon %d:%d", p, q);
			sprintf(filename, "../share/families/EM_Res_%d_%d.mat", p, q);
			res.setName(famName);
			res.saveToMat(filename);
			break;
		}
		case 201:
		{
			printf("Generating Earth-Moon Resonant Orbit via Pseudo Arclength\n");
			gen.setContType(tpat_family_generator::PSEUDO_ARC);
			gen.setNumOrbits(1000);
			
			tpat_family_cr3bp res = gen.cr3bp_generateRes(&sysData, p, q);
			res.sortMembers();
			res.sortEigs();

			char famName[32], filename[128];
			sprintf(famName, "Earth-Moon %d:%d", p, q);
			sprintf(filename, "../share/families/EM_Res_%d_%d_PAC.mat", p, q);
			res.setName(famName);
			res.saveToMat(filename);
			break;
		}
		case 999:
		{
			printf("Generating Earth-Moon L2 Northern Butterfly Family\n");
			// Natural Parameter Continuation
			gen.setStep_simple(0.0001);
			gen.setStep_fitted_1(0.0005);
			gen.setStep_fitted_2(0.0005);
			gen.setNumNodes(4);

			// PAC
			gen.setContType(tpat_family_generator::PSEUDO_ARC);
			gen.setNumOrbits(5000);
			gen.setNumNodes(5);

			tpat_family_cr3bp L2_NButterfly = gen.cr3bp_generateButterfly(&sysData, 2);
			L2_NButterfly.sortEigs();
			L2_NButterfly.setName("Earth-Moon L2 Northern Butterfly");
			L2_NButterfly.saveToMat("../share/families/EM_L2_NButterfly.mat");
			break;
		}
		default:
			printf("The problem number %d has not been implemented\n", probNum);
			break;
	}

	return EXIT_SUCCESS;
}
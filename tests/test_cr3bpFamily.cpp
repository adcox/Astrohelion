/**
 *	Test out CR3BP Families of Orbits
 */
#include "tpat_ascii_output.hpp"
#include "tpat_family_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"

#include <iostream>

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

int main(void){
	// Load the family
	// tpat_family_cr3bp fam("../share/families/EM_L1_Lyap.mat");
	tpat_family_cr3bp fam("../share/families/EM_L2_NButterfly_combined.mat");
	fam.sortEigs();
	std::vector<int> bifs = fam.findBifurcations();
	if(bifs.size() > 0){
		printf("Found bifurcations at:\n");
		for(size_t i = 0; i < bifs.size(); i++){
			printf("  Ix = %04d\n", bifs[i]);
		}
	}
	fam.saveToMat("data/LoadedButterflyFam.mat");	// Check to see if data was re-loaded correctly

	printf("Checing Match State: X\n");
	double matchX = 0.9;
	std::vector<tpat_family_member_cr3bp> matches = fam.getMemberByStateVar(matchX, 0);
	printf("  Found %zu Potential members\n", matches.size());
	for(size_t i = 0; i < matches.size(); i++){
		printf("   %03zu: x = %f ", i, matches[i].getIC()[0]);
		std::cout << (std::abs(matches[i].getIC()[0] - matchX) < 1e-9 ? PASS : FAIL) << std::endl;
	}

	printf("Checking Match Jacobi:\n");
	double matchJC = 3;
	matches.clear();
	matches = fam.getMemberByJacobi(matchJC);
	printf("  Found %zu Potential members\n", matches.size());
	for(size_t i = 0; i < matches.size(); i++){
		printf("   %03zu: JC = %f ", i, matches[i].getJacobi());
		std::cout << (std::abs(matches[i].getJacobi() - matchJC) < 1e-9 ? PASS : FAIL) << std::endl;
	}

	printf("Checking Match TOF:\n");
	double matchTOF = 4.5;
	matches.clear();
	matches = fam.getMemberByTOF(matchTOF);
	printf("  Found %zu Potential members\n", matches.size());
	for(size_t i = 0; i < matches.size(); i++){
		printf("   %03zu: TOF = %f ", i, matches[i].getTOF());
		std::cout << (std::abs(matches[i].getTOF() - matchTOF) < 1e-9 ? PASS : FAIL) << std::endl;
	}
}
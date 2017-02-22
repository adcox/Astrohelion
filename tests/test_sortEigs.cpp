/**
 * Test the eigenvalue sorting algorithm on a family file
 */

#include "Fam_cr3bp.hpp"

#include <iostream>

using namespace astrohelion;

int main(){

	// Fam_cr3bp fam("../share/families_natParam_checked/EM_L1_NAxial.mat");
	// Fam_cr3bp fam("../share/families/EM_L1_Lyap.mat");
	Fam_cr3bp fam("../../Astrohelion_scripts/share/families/SE_DPO_PAC.mat");

	// fam.setSortType(FamSort_tp::SORT_TOF);
	// fam.sortMembers();
 	fam.sortEigs();

 	fam.saveToMat("../../Astrohelion_scripts/share/EigSorted_TestFam.mat");

 	return EXIT_SUCCESS;
 }
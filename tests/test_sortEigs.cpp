/**
 * Test the eigenvalue sorting algorithm on a family file
 */

#include "tpat_family_cr3bp.hpp"

#include <iostream>


int main(){

	// tpat_family_cr3bp fam("../share/families_natParam_checked/EM_L1_NAxial.mat");
	tpat_family_cr3bp fam("../share/families/EM_L1_Lyap.mat");

	// fam.sortMembers();
 	fam.sortEigs();

 	fam.saveToMat("../share/EigSorted_TestFam_Axial.mat");

 	return EXIT_SUCCESS;
 }
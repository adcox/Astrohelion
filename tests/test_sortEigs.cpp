/**
 * Test the eigenvalue sorting algorithm on a family file
 */

#include "tpat_fam_cr3bp.hpp"

#include <iostream>


int main(){

	// TPAT_Fam_CR3BP fam("../share/families_natParam_checked/EM_L1_NAxial.mat");
	TPAT_Fam_CR3BP fam("../share/families/EM_L1_Lyap.mat");

	// fam.sortMembers();
 	fam.sortEigs();

 	fam.saveToMat("../share/EigSorted_TestFam_Axial.mat");

 	return EXIT_SUCCESS;
 }
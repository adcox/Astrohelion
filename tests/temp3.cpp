/**
 * @file ltpo_manifolds.cpp
 * This script computes stable and unstable manifolds of a low-thrust periodic
 * orbit.
 * 
 * Usage:
 * 
 * 	ltpo_manifolds <fam_file> <orb_ix> <numMan> <tof> 
 * 		Compute the manifolds of an orbit from the family specified by fam_file
 * 		at index orb_ix (orb_ix = 0 is the first family member). To specify the 
 * 		number of manifolds, append a positive integer as numMan and a positive 
 * 		double as tof
 * 	
 * 
 * @author Andrew Cox
 * @version June 14, 2018
 */
#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;
using ltlaw = ControlLaw_cr3bp_lt;

double stepDist = 50;	// kilometers

// Function declarations
void freeMem(std::vector<ControlLaw*>&);

// Function definitions
void freeMem(std::vector<ControlLaw*>& laws){
	for(ControlLaw* pLaw : laws){
		delete pLaw;
		pLaw = nullptr;
	}
}//====================================================


int main(){
	
	SysData_cr3bp sys("../../data/families/cr3bp_earth-moon/L1_Lyap.mat");
	Family_PO_cr3bp fam(&sys);
	std::vector<ControlLaw*> laws;
	fam.readFromMat("../../data/families/cr3bp_earth-moon/L1_Lyap.mat", laws);

	printf("Checing Match State: X\n");
	double matchX = 0.9;
	std::vector<Arcset_periodic> matches = fam.getMemberByState(matchX, 0);

	matches[0].print();
	
	freeMem(laws);
	return EXIT_SUCCESS;
}
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
	
	char filename[] = "../../data/families/cr3bp-lt_earth-moon/"
		"E2_Lyap_f7.0e-02_alph000.00_law2112.mat";
	SysData_cr3bp_lt sys(filename);
	Family_PO_cr3bp_lt fam(&sys);
	std::vector<ControlLaw*> laws;
	fam.readFromMat(filename, laws);

	double H = -1.55;
	std::vector<Arcset_periodic> matches = fam.getMemberByH_lt(H);
	for(unsigned int i = 0; i < matches.size(); i++){
		Arcset_cr3bp_lt temp(matches[i]);
		if(std::abs(temp.getHltByIx(0) - H) > 1e-9)
			printErr("Hlt is not close to desired value");
	}
	
	std::vector<double> elem = matches[0].getSTMElementsByIx(0);
	std::cout << elem.size() << std::endl;
	std::cout << std::sqrt(elem.size()) << std::endl;
	
	freeMem(laws);
	return EXIT_SUCCESS;
}
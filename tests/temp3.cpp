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

#include <cstdlib>
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

void mexErrMsgTxt(const char* msg){
	printErr(msg);
}


int main(){
	SimEngine sim;
	std::vector<double> q0 {0.5, 0, 0, 0, 1.5, 0, 1};
	std::vector<double> ctrl0 {0, 0};

	SysData_cr3bp_lt sys("earth", "moon", 1);
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::CONST_MF_GENERAL,
		std::vector<double> {0.07});

	Arcset_cr3bp_lt arc(&sys);

	Event evt(Event_tp::YZ_PLANE, 0, true);

	sim.addEvent(evt);
	sim.setSimpleInt(true);
	// sim.setVerbosity(Verbosity_tp::DEBUG);
	sim.runSim(q0, ctrl0, 0, 6*PI, &arc, &law);

	arc.print();

	return EXIT_SUCCESS;
}
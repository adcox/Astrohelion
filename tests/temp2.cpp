/**
 * @file ltpo_manifolds.cpp
 * This script computes stable and unstable manifolds of a low-thrust periodic
 * orbit at a specified low-thrust Hamiltonian value
 * 
 * Usage:
 * 
 * 	ltpo_manifolds <fam_file> <Hlt> <numMan> <tof> 
 * 		Compute the manifolds of an orbit from the family specified by fam_file
 * 		with low-thrust Hamiltonian equal to Hlt. To specify the 
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


int main(int argc, char** argv){
	// Input 1: famFile
	// Input 2: low-thrust Hamiltonian
	// Input 2: numMan
	// Input 3: tof (nondim)

	if(argc < 3){
		printErr("Too few arguments. See script documentation.\n");
		return EXIT_SUCCESS;
	}

	const char* famFile = argv[1];
	double Hlt = std::atof(argv[2]);
	int numMan = 50;
	int tof = 4*PI;

	if(argc > 3)
		numMan = std::atoi(argv[3]);

	if(argc > 4)
		tof = std::atof(argv[4]);

	// Error handling

	if(tof < 0){
		printErr("tof = %f is invalid; must be positive\n");
		return EXIT_SUCCESS;
	}

	SysData_cr3bp_lt ltSys(famFile);
	Family_PO_cr3bp_lt fam(&ltSys);

	std::vector<ControlLaw*> loadedLaws {};
	fam.readFromMat(famFile, loadedLaws);

	std::vector<Arcset_periodic> matches = fam.getMemberByH_lt(Hlt);
	if(matches.empty()){
		printErr("Could not find an orbit with Hlt = %f\n", Hlt);
		freeMem(loadedLaws);
		return EXIT_SUCCESS;
	}
	Arcset_cr3bp_lt famArc = static_cast<Arcset_cr3bp_lt>(matches[0]);

	assert(famArc.getCtrlLawByIx(0) == loadedLaws[0]);
	ControlLaw_cr3bp_lt *pLaw = static_cast<ControlLaw_cr3bp_lt*>(loadedLaws[0]);

	char outFile[128];
	sprintf(outFile, "orb_H%.4f.mat", Hlt);
	famArc.saveToMat(outFile, Save_tp::SAVE_CURVE);

	std::vector<double> ctrl0 = famArc.getExtraParamVecByIx(0, PARAMKEY_CTRL);

	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	Arcset_cr3bp_lt singleSegArc(&ltSys);

	sim.runSim_manyNodes(famArc.getStateByIx(0), ctrl0, 0, famArc.getTotalTOF(), 
		numMan, &singleSegArc, famArc.getCtrlLawByIx(0));
	singleSegArc.setSTMs_sequence();

	// Compute the manifold initial conditions (no current support for 
	// propagating with low-thrust in the ManifoldEngine)
	ManifoldEngine man;
	man.setStepOffDist(stepDist);
	std::vector<Arcset_cr3bp_lt> ics_u, ics_s;
	try{
		ics_u = man.computeSetFromLTPeriodic(Manifold_tp::MAN_U, &singleSegArc,
			pLaw, numMan, 0.1);
		ics_s = man.computeSetFromLTPeriodic(Manifold_tp::MAN_S, &singleSegArc, 
			pLaw, numMan, 0.1);
	}catch(const Exception &e){
		freeMem(loadedLaws);
		printErr("%s\n", e.what());
		return EXIT_SUCCESS;
	}

	// Construct two families of orbits: one for the unstable manifolds, another
	// for the stable manifolds
	Family_PO_cr3bp_lt manifolds_u(&ltSys), manifolds_s(&ltSys);

	// Propagate the unstable manifold ICs in forward time
	Arcset_cr3bp_lt temp(&ltSys);
	for(unsigned int m = 0; m < ics_u.size(); m++){
		printf("Unstable Manifold %03u\n", m);
		temp.reset();
		sim.runSim(ics_u[m].getStateByIx(0),
			ics_u[m].getExtraParamVecByIx(0, PARAMKEY_CTRL), 0, tof, &temp, pLaw);

			manifolds_u.addMember(temp);
	}

	// Propagate the stable manifold ICs in reverse time
	sim.setRevTime(true);
	for(unsigned int m = 0; m < ics_s.size(); m++){
		printf("Stable Manifold %03u\n", m);
		temp.reset();
		sim.runSim(ics_s[m].getStateByIx(0),
			ics_s[m].getExtraParamVecByIx(0, PARAMKEY_CTRL), 0, tof, &temp, pLaw);

			manifolds_s.addMember(temp);
	}

	sprintf(outFile, "orb_H%.4f_unstableManifolds.mat", Hlt);
	manifolds_u.saveToMat(outFile, Save_tp::SAVE_CURVE);
	sprintf(outFile, "orb_H%.4f_stableManifolds.mat", Hlt);
	manifolds_s.saveToMat(outFile, Save_tp::SAVE_CURVE);

	freeMem(loadedLaws);
	return EXIT_SUCCESS;
}
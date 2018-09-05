/**
 * @file natpo_manifolds_atHnat.cpp
 * This script computes stable and unstable manifolds of a natural periodic
 * orbit at a specified natural Hamiltonian value
 * 
 * Usage:
 * 
 * 	natpo_manifolds_atHnat <fam_file> <Hnat> <numMan> <tof> 
 * 		Compute the manifolds of an orbit from the family specified by fam_file
 * 		with natural Hamiltonian equal to Hnat. To specify the 
 * 		number of manifolds, append a positive integer as numMan and a positive 
 * 		double as tof
 * 	
 * 
 * @author Andrew Cox
 * @version Sep 5, 2018
 */
#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

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
	// Input 2: natural-thrust Hamiltonian
	// Input 3: numMan
	// Input 4: tof (nondim)

	if(argc < 3){
		printErr("Too few arguments. See script documentation.\n");
		return EXIT_SUCCESS;
	}

	const char* famFile = argv[1];
	double Hnat = std::atof(argv[2]);
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

	SysData_cr3bp natSys(famFile);
	Family_PO_cr3bp fam(&natSys);

	std::vector<ControlLaw*> loadedLaws {};
	fam.readFromMat(famFile, loadedLaws);

	std::vector<Arcset_periodic> matches = fam.getMemberByJacobi(-2*Hnat);
	if(matches.empty()){
		printErr("Could not find an orbit with Hnat = %f\n", Hnat);
		freeMem(loadedLaws);
		return EXIT_SUCCESS;
	}
	Arcset_cr3bp famArc = static_cast<Arcset_cr3bp>(matches[0]);

	char outFile[128];
	sprintf(outFile, "orb_H%.4f.mat", Hnat);
	famArc.saveToMat(outFile, Save_tp::SAVE_CURVE);

	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	Arcset_cr3bp singleSegArc(&natSys);

	sim.runSim_manyNodes(famArc.getStateByIx(0), famArc.getTotalTOF(), 
		numMan, &singleSegArc);
	// singleSegArc.setSTMs_sequence();
	singleSegArc.setSTMs_cumulative();

	// Compute the manifolds
	ManifoldEngine man;
	man.setStepOffDist(stepDist);
	man.setVerbosity(Verbosity_tp::NO_MSG);
	std::vector<Arcset_cr3bp> arcs_u, arcs_s;
	try{
		arcs_u = man.computeSetFromPeriodic(Manifold_tp::MAN_U, &singleSegArc,
			numMan, 0.1);
		arcs_s = man.computeSetFromPeriodic(Manifold_tp::MAN_S, &singleSegArc, 
			numMan, 0.1);
	}catch(const Exception &e){
		freeMem(loadedLaws);
		printErr("%s\n", e.what());
		return EXIT_SUCCESS;
	}

	// Construct two families of orbits: one for the unstable manifolds, another
	// for the stable manifolds
	Family_PO_cr3bp manifolds_u(&natSys), manifolds_s(&natSys);

	// Propagate the unstable manifold ICs in forward time
	// assert(arcs_u.size() == arcs_s.size());
	// for(unsigned int m = 0; m < arcs_u.size(); m++){
	// 	manifolds_u.addMember(arcs_u[m]);
	// 	manifolds_s.addMember(arcs_s[m]);
	// }

	// Propagate the unstable manifold ICs in forward time
	Arcset_cr3bp temp(&natSys);
	for(unsigned int m = 0; m < arcs_u.size(); m++){
		printf("Unstable Manifold %03u\n", m);
		temp.reset();
		sim.runSim(arcs_u[m].getStateByIx(0), tof, &temp);

		manifolds_u.addMember(temp);
	}

	// Propagate the stable manifold ICs in reverse time
	sim.setRevTime(true);
	for(unsigned int m = 0; m < arcs_s.size(); m++){
		printf("Stable Manifold %03u\n", m);
		temp.reset();
		sim.runSim(arcs_s[m].getStateByIx(0), tof, &temp);

		manifolds_s.addMember(temp);
	}

	sprintf(outFile, "orb_H%.4f_unstableManifolds.mat", Hnat);
	manifolds_u.saveToMat(outFile, Save_tp::SAVE_CURVE);
	sprintf(outFile, "orb_H%.4f_stableManifolds.mat", Hnat);
	manifolds_s.saveToMat(outFile, Save_tp::SAVE_CURVE);

	freeMem(loadedLaws);
	return EXIT_SUCCESS;
}
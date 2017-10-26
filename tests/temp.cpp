#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

int main(){
	double C0 = 3.1;
	double M0 = 1000;
	double F = 0, f = 1e-2;
	double Isp = 1500;

	SysData_cr3bp emSys("earth", "moon");
	SysData_cr3bp_lt sys("earth", "moon", M0);
	Family_PO_cr3bp emFam(&emSys);
	std::vector<ControlLaw*> loadedLaws;
	emFam.readFromMat("../../data/families/cr3bp_earth-moon/L1_Lyap.mat", loadedLaws);

	std::vector<Arcset_periodic> matches = emFam.getMemberByJacobi(C0);

	if(matches.size() == 0){
		printErr("Could not find/correct any orbits!\n");
		return 0;
	}

	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::CONST_C_2D_LEFT, F, Isp);
	

	// Generate the natural Lyapunov orbit
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);

	Arcset_cr3bp_lt natPO(&sys);
	std::vector<double> ic = matches[0].getStateByIx(0);
	ic.push_back(1);

	sim.runSim_manyNodes(ic, matches[0].getTotalTOF(), 5, &natPO, &law);

	// Form constraints for periodicity
	double n = static_cast<double>(natPO.getNodeByIx(-1).getID());
	std::vector<double> periodicityData {0, 0, NAN, 0, 0, NAN, NAN};
	Constraint periodicityCon(Constraint_tp::MATCH_CUST, n, periodicityData);

	std::vector<double> initStateData {NAN, 0, 0, NAN, NAN, NAN, 1};
	Constraint fixInitState(Constraint_tp::STATE, 0, initStateData);

	natPO.addConstraint(periodicityCon);
	natPO.addConstraint(fixInitState);

	law.setThrust(f);

	Arcset_cr3bp_lt a1(&sys), a2(&sys);
	Family_PO ltFam(&sys);

	std::vector<Arcset> allArcs {};
	std::vector<int> initDir {1, 0, 0, 0, 0, 0};

	MultShootEngine temp;
	temp.setDoLineSearch(true);
	temp.setMaxIts(200);

	PseudoArcEngine pae;
	pae.setTol(9e-12);
	pae.setStepCountIncrease(10);
	// pae.setNumOrbits(50);
	// pae.setVerbosity(Verbosity_tp::ALL_MSG);
	// pae.setNullspaceCol(1);
	pae.pac(&natPO, &a1, &a2, allArcs, initDir, &temp);

	for(unsigned int i = 0; i < allArcs.size(); i++){
		ltFam.addMember(static_cast<Arcset_cr3bp_lt>(allArcs[i]));
	}

	ltFam.saveToMat("../../data/families_toBeChecked/ltFamTest.mat");

	for(unsigned int i = 0; i < loadedLaws.size(); i++){
		delete loadedLaws[i];
		loadedLaws[i] = nullptr;
	}
}
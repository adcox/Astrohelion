#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

int main(){
	// Load the L1 Lyapunov family
	char filename[128];
	sprintf(filename, "../../data/families/cr3bp-lt_earth-moon/"
		"L1_Lyap_f7.0e-02_Hlt-1.542_law105.mat");
	printf("Opening %s\n", filename);
	SysData_cr3bp_lt ltSys(filename);

	Family_PO_cr3bp_lt fam(&ltSys);
	std::vector<ControlLaw*> loadedLaws {};
	fam.readFromMat(filename, loadedLaws);

	for(unsigned int i = 0; i < loadedLaws.size(); i++){
		loadedLaws[i]->print();
	}

	unsigned int ix = std::floor(3*fam.getNumMembers()/8);
	Arcset_periodic arc = fam.getMember(ix);
	Arcset_cr3bp_lt arc_fixedMass(&ltSys);
	reconstructArc(&arc, &arc_fixedMass);

	// arc.print();
	// arc_fixedMass.print();
	arc_fixedMass.saveToMat("LTPO_fixedMass.mat");

	// Copy the fixed mass case
	Arcset_cr3bp_lt arc_varMass = arc_fixedMass;

	// Change control law
	if(loadedLaws.empty()){
		printErr("Did not load any control laws?\n");
		return EXIT_SUCCESS;
	}

	loadedLaws[0]->setType(ControlLaw_cr3bp_lt::CONST_F_GENERAL);

	// Remove control constraints
	std::vector<Constraint> allCons = arc_varMass.getAllConstraints();
	arc_varMass.clearAllConstraints();
	for(unsigned int i = 0; i < allCons.size(); i++){
		if(allCons[i].getType() != Constraint_tp::CONT_CTRL &&
			allCons[i].getType() != Constraint_tp::CTRL){

			arc_varMass.addConstraint(allCons[i]);
		}
	}

	// Reconverge LTPO
	MultShootEngine shooter;
	shooter.setVerbosity(Verbosity_tp::SOME_MSG);
	shooter.setTOFType(MSTOF_tp::VAR_EQUALARC);

	// Be extra careful when converging the first orbit
	// shooter.setMaxIts(200);
	// shooter.setDoLineSearch(true);

	Arcset_cr3bp_lt converged(&ltSys);
	try{
		shooter.multShoot(&arc_varMass, &converged);
	}catch(DivergeException &e){
		printErr("Did not converge:\n%s\n", e.what());
		return EXIT_SUCCESS;
	}
	
	converged.saveToMat("LTPO_varMass.mat");

	// Free memory
	for(ControlLaw *p : loadedLaws){
		delete p;
		p = nullptr;
	}
}
#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;
using ltlaw = ControlLaw_cr3bp_lt;

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

	char famFile[] = "../../data/families/cr3bp_earth-moon/L2_Halo.mat";

	std::vector<ControlLaw*> loadedLaws {};
	SysData_cr3bp sys(famFile);
	Family_PO_cr3bp fam(&sys);

	double T_syn = 29.5406; 	// days
	double T = 2*T_syn/9 * 24 * 3600 / sys.getCharT();	// 9:2 synodic period, nondim

	fam.readFromMat(famFile, loadedLaws);
	std::vector<Arcset_periodic> matches = fam.getMemberByTOF(T);

	if(!matches.empty()){
		Arcset_cr3bp arc = static_cast<Arcset_cr3bp>(matches[0]);

		arc.clearArcConstraints();	// get rid of the pseudo-arclength constraint

		Constraint tofCon(Constraint_tp::TOF_TOTAL, 0, std::vector<double> {T});
		arc.addConstraint(tofCon);

		int id0 = arc.getNodeRefByIx(0).getID();
		std::vector<double> periodicData {NAN, NAN, id0, NAN, id0, NAN};
		Constraint periodicCon(Constraint_tp::MATCH_CUST, 
			arc.getNodeRefByIx(-1).getID(), periodicData);
		// arc.addConstraint(periodicCon);

		arc.print();

		Arcset_cr3bp correctArc(&sys);
		MultShootEngine shooter;
		shooter.setVerbosity(Verbosity_tp::ALL_MSG);
		shooter.setTOFType(MSTOF_tp::VAR_EQUALARC);
		shooter.multShoot(&arc, &correctArc);

		correctArc.saveToMat("NRHO_9-2_Synodic.mat");

		// Arcset_cr3bp fullArc(&sys);
		// reconstructArc(&arc, &fullArc);
		// fullArc.saveToMat("NRHO_9-2_Synodic.mat");
	}else{
		printWarn("Could not find a Halo at T = %.2f\n", T);
	}

	freeMem(loadedLaws);
	return EXIT_SUCCESS;
}//====================================================
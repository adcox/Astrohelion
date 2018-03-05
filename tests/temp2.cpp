#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

int main(){
	// Load the L1 Lyapunov family
	char filename[128];
	sprintf(filename, "../../data/families/cr3bp_earth-moon/L1_Lyap.mat");
	printf("Opening %s\n", filename);
	SysData_cr3bp crSys(filename);

	Family_PO lyapFam(&crSys);
	std::vector<ControlLaw*> loadedLaws {};
	lyapFam.readFromMat(filename, loadedLaws);

	unsigned int ix0 = 5;

	Arcset_periodic lyap = lyapFam.getMember(ix0);

	// Define low-thrust parameters
	double f = 7e-2;
	double Isp = 1500;

	// Create parameters and other data storage objects
	SysData_cr3bp_lt ltSys("earth", "moon", 1);	// Let M0 = 1
	std::vector<double> ltParams {sqrt(f), Isp},
		ctrl0 {0, 0},
		n0conData {NAN, 0, 0, NAN, NAN, 0, 1},
		ctrlConData {1, 1};
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::CONST_MF_GENERAL, ltParams);

	// Create sim engine and storage arcset
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	Arcset_cr3bp_lt ltGuess(&ltSys), ltConverged(&ltSys);

	// Create Multiple shooting engine
	MultShootEngine shooter;
	shooter.setVerbosity(Verbosity_tp::ALL_MSG);
	shooter.setTOFType(MSTOF_tp::VAR_EQUALARC);
	Family_PO fam(&ltSys);

	double alpha = PI/3;
	
	printf("Converging at alpha = %.2f deg\n", alpha*180/PI);

	ctrl0[0] = alpha;
	ltGuess.reset();
	ltConverged.reset();

	// Fill the storage arcset with arcs propagated from the natural Lyapunov nodes
	// with the specified low-thrust parameters
	double m = 1;
	std::vector<double> q0 {};
	for(unsigned int i = 0; i < lyap.getNumSegs(); i++){
		q0 = lyap.getStateByIx(i);
		q0.push_back(m);
		if(i == 0){
			sim.runSim(q0, ctrl0, 0, lyap.getTOF(i), &ltGuess, &law);
			
		}else{
			Arcset_cr3bp_lt temp(&ltSys);
			sim.runSim(q0, ctrl0, 0, lyap.getTOF(i), &temp, &law);
			ltGuess.appendSetAtNode(&temp, ltGuess.getNodeByIx(-1).getID(), 0, 0);
		}

		m = ltGuess.getStateByIx(-1).back();
	}

	// Periodicity Constraint
	int id0 = ltGuess.getNodeByIx(0).getID();
	int idf = ltGuess.getNodeByIx(-1).getID();

	std::vector<double> conData{ idf, idf, NAN, idf, NAN, NAN, NAN};
	Constraint periodicityCon(Constraint_tp::MATCH_CUST, id0, conData);

	ltGuess.addConstraint(periodicityCon);

	// Fix initial state components
	Constraint n0con(Constraint_tp::STATE, ltGuess.getNodeByIx(0).getID(), n0conData);
	ltGuess.addConstraint(n0con);

	// Remove control from design vector; always fixed
	for(unsigned int i = 0; i < ltGuess.getNumNodes(); i++){
		Constraint con(Constraint_tp::RM_CTRL, ltGuess.getNodeByIx(i).getID(), ctrlConData);
		ltGuess.addConstraint(con);
	}

	try{
		shooter.multShoot(&ltGuess, &ltConverged);
	}catch(DivergeException &e){
		printErr("Corrector diverged\n");
	}

	ltConverged.putInChronoOrder();

	// ltGuess.print();
	// ltConverged.print();

	PseudoArcEngine pae;
	std::vector<Arcset> allArcs;
	std::vector<int> initDir {1, 0, 0, 0, 0, 0};
	Arcset_cr3bp_lt a1(&ltSys), a2(&ltSys);

	pae.setNumOrbits(100);
	pae.pac(&ltConverged, &a1, &a2, allArcs, initDir);

	for(auto arc : allArcs)
		fam.addMember(arc);

	fam.sortEigs();

	sprintf(filename, "L1_Lyap_f%.1e_a%05.1f_law%d.mat", f,
		alpha*180/PI, law.getType());
	fam.saveToMat(filename);

	// Free memory
	for(ControlLaw *p : loadedLaws){
		delete p;
		p = nullptr;
	}
}
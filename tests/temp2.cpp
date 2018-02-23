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

	unsigned int ix0 = 50;	// for lawID = 2, f = 1e-2

	Arcset_periodic member = lyapFam.getMember(ix0);

	// Define low-thrust parameters
	double f = 1e-2;
	double Isp = 1500;
	double alpha = PI/3;

	// Create parameters and other data storage objects
	SysData_cr3bp_lt ltSys("earth", "moon", 1);	// Let M0 = 1
	std::vector<double> ltParams {f, Isp}, ctrl0 {alpha, 0};
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::CONST_F_GENERAL, ltParams);
	
	// Create sim engine and storage arcset
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	Arcset_cr3bp_lt ltGuess(&ltSys);

	// Fill the storage arcset with arcs propagated from the natural Lyapunov nodes
	// with the specified low-thrust parameters
	double m = 1;
	std::vector<double> q0 {};
	for(unsigned int i = 0; i < member.getNumSegs(); i++){
		q0 = member.getStateByIx(i);
		q0.push_back(m);
		if(i == 0){
			sim.runSim(q0, ctrl0, 0, member.getTOF(i), &ltGuess, &law);
			
		}else{
			Arcset_cr3bp_lt temp(&ltSys);
			sim.runSim(q0, ctrl0, 0, member.getTOF(i), &temp, &law);
			ltGuess.appendSetAtNode(&temp, ltGuess.getNodeByIx(-1).getID(), 0, 0);
		}

		m = ltGuess.getStateByIx(-1).back();
	}

	ltGuess.print();
	ltGuess.saveToMat("ltGuess.mat");

	// Now, try corrections
	MultShootEngine shooter;
	shooter.setDoLineSearch(true);
	shooter.setMaxIts(200);

	int id0 = ltGuess.getNodeByIx(0).getID();
	int idf = ltGuess.getNodeByIx(-1).getID();

	std::vector<double> conData{ idf, idf, NAN, idf, idf, NAN, NAN};
	Constraint periodicityCon(Constraint_tp::MATCH_CUST, id0, conData);

	ltGuess.addConstraint(periodicityCon);

	// // Enforce control continuity on all segment
	// std::vector<double> ctrlConData {1, 1};
	// for(unsigned int i = 0; i < ltGuess.getNumSegs(); i++){
	// 	Constraint con(Constraint_tp::CONT_CTRL, ltGuess.getSegByIx(i).getID(), ctrlConData);
	// 	ltGuess.addConstraint(con);
	// }

	Arcset_cr3bp_lt ltConverged(&ltSys);
	try{
		shooter.multShoot(&ltGuess, &ltConverged);
	}catch(DivergeException &e){
		printErr("Corrector diverged\n");
	}

	ltConverged.putInChronoOrder();
	ltConverged.saveToMat("ltConverged.mat");

	// double L4pos[3];
	// DynamicsModel_cr3bp::getEquilibPt(&ltSys, 4, 1e-6, L4pos);

	// std::vector<double> yVal{lawID == 1 ? L4pos[1] : -L4pos[1]};
	// Event planeEvt(Event_tp::XZ_PLANE, 0, true, yVal);

	// printf("Stop at y = %f\n", yVal[0]);
	// Arcset_cr3bp_lt temp(&ltSys), newMember(&ltSys);

	// SimEngine sim;
	// sim.setVerbosity(Verbosity_tp::NO_MSG);
	// sim.addEvent(planeEvt);

	// sim.runSim(member.getStateByIx(0), member.getTotalTOF(), &temp, loadedLaws[0]);

	// sim.clearEvents();
	// sim.runSim_manyNodes(temp.getStateByIx(-1), member.getTotalTOF(), 5, &newMember, loadedLaws[0]);

	// std::vector<double> initStateData {NAN, yVal[0], 0, NAN, NAN, 0, 1};
	// std::vector<double> perData {0, 0, NAN, 0, NAN, NAN};

	// Constraint periodicityCon(Constraint_tp::MATCH_CUST, newMember.getNodeByIx(-1).getID(), perData);
	// Constraint fixInitState(Constraint_tp::STATE, newMember.getNodeByIx(0).getID(), initStateData);

	// newMember.addConstraint(fixInitState);
	// newMember.addConstraint(periodicityCon);

	// newMember.saveToMat("data/tempInitGuess.mat");

	// MultShootEngine shooter;
	// Arcset_cr3bp_lt corrected(&ltSys);

	// try{
	// 	shooter.multShoot(&newMember, &corrected);
	// }catch(Exception &e){
	// 	printErr("%s\n", e.what());
	// }

	// corrected.saveToMat("data/tempCorrected.mat");

	// PseudoArcEngine pae;
	// std::vector<Arcset> allArcs;
	// std::vector<int> initDir {1, 0, 0, 1, -1, 0};
	// Arcset_cr3bp_lt a1(&ltSys), a2(&ltSys);

	// pae.setNumOrbits(50);
	// pae.pac(&newMember, &a1, &a2, allArcs, initDir);

	// Family_PO fam(&ltSys);
	// for(Arcset arc : allArcs)
	// 	fam.addMember(static_cast<Arcset_cr3bp_lt>(arc));

	// sprintf(filename, "../../data/families_toBeChecked/cr3bp-lt_L4Lyap_f%.1e_Isp3000_law%u.mat", f, lawID);
	// fam.saveToMat(filename);

	// Free memory
	for(ControlLaw *p : loadedLaws){
		delete p;
		p = nullptr;
	}
}
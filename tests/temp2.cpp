#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

int main() {
	double C0 = 1.6;
	double M0 = 1000;
	double F = 0, f = 1e-2;
	double Isp = 1500;
	double xStep = 5e-3;
	bool bDoNewMethod = true;

	SysData_cr3bp cr3bpSys("earth", "moon");
	SysData_cr3bp_lt ltSys("earth", "moon", M0);

	Family_PO_cr3bp lyapFam(&cr3bpSys);
	std::vector<ControlLaw*> loadedLaws {};
	lyapFam.readFromMat("../../data/families/cr3bp_earth-moon/L3_Lyap.mat", loadedLaws);

	std::vector<Arcset_periodic> matches = lyapFam.getMemberByJacobi(C0);

	if(matches.size() == 0){
		printErr("Could not find/correct any orbits!\n");
		return 0;
	}

	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::CONST_C_2D_RIGHT, F, Isp);

	// Generate the natural Lyapunov orbit
	SimEngine sim;
	sim.setVarStepSize(false);
	sim.setNumSteps(2);
	sim.setMakeDefaultEvents(false);
	sim.setVerbosity(Verbosity_tp::NO_MSG);

	Arcset_cr3bp_lt natPO(&ltSys);
	std::vector<double> ic = matches[0].getStateByIx(0);
	ic.push_back(1);

	for(double tof = 1; tof < 1000; tof += 50){
		natPO.reset();
		try{
			sim.runSim(ic, tof, &natPO, &law);		
		}catch(Exception &e){
			printWarn("Exception: %s\n", e.what());
		}
	}
	// sim.runSim_manyNodes(ic, matches[0].getTotalTOF(), 7, &natPO, &law);

	// // Form constraints for periodicity
	// std::vector<double> perData {0, 0, NAN, 0, 0, NAN};
	// Constraint periodicityCon(Constraint_tp::MATCH_CUST, natPO.getNodeByIx(-1).getID(), perData);

	// Constraint fixJC(Constraint_tp::JC, 0, &C0, 1);

	// std::vector<double> initStateData {NAN, 0, NAN, NAN, NAN, NAN, 1};
	// Constraint fixInitState(Constraint_tp::STATE, 0, initStateData);

	// natPO.addConstraint(periodicityCon);
	// natPO.addConstraint(fixJC);
	// natPO.addConstraint(fixInitState);

	// law.setThrust(f);

	// Arcset_cr3bp_lt corrected(&ltSys), a1(&ltSys), a2(&ltSys);

	// // Correct natural orbit with low-thrust to be periodic
	// MultShootEngine shooter;

	// shooter.setDoLineSearch(true);
	// shooter.setMaxIts(40);
	// shooter.setTol(1e-11);
	// shooter.setIgnoreCrash(true);
	// shooter.setTOFType(MSTOF_tp::VAR_FIXSIGN);
	// shooter.multShoot(&natPO, &corrected);

	// ControlLaw_cr3bp_lt outputLaw(ControlLaw_cr3bp_lt::GENERAL_CONST_F, f, Isp);

	// // Create family to store results in
	// Family_PO ltFam(&ltSys);
	
	// if(bDoNewMethod){
	// 	// Reset constraints; get rid of Jacobi constraint
	// 	corrected.clearAllConstraints();
	// 	corrected.addConstraint(periodicityCon);

	// 	std::vector<double> alwaysFixVals {NAN, 0, NAN, NAN, NAN, NAN, 1};
	// 	std::vector<unsigned int> indVarIx {0, 4};
	// 	std::vector<unsigned int> depVarIx {4};
	// 	std::vector<Arcset> allArcs {};

	// 	MultShootEngine engineTemplate;
	// 	engineTemplate.setTOFType(MSTOF_tp::VAR_EQUALARC);
	// 	engineTemplate.setIgnoreCrash(false);

	// 	NatParamEngine npe;
	// 	npe.setStep_simple(5e-3);
	// 	npe.setTol(5e-12);
	// 	npe.setNumOrbits(25);
	// 	npe.continuePO(&corrected, &a1, &a2, allArcs, alwaysFixVals, indVarIx, depVarIx, &engineTemplate);

	// 	for(unsigned int i = 0; i < allArcs.size(); i++){
	// 		Arcset_cr3bp_lt temp(allArcs[i]);
	// 		// ControlLaw_cr3bp_lt::convertLaws(&temp, &outputLaw);
	// 		ltFam.addMember(temp);
	// 	}
	// }else{

	// 	// Reset constraints; get rid of Jacobi constraint
	// 	corrected.clearAllConstraints();
	// 	corrected.addConstraint(periodicityCon);
	// 	corrected.addConstraint(fixInitState);

	// 	// Original algorithm
	// 	bool doContinuation = true;
		
	// 	Arcset_cr3bp_lt nextArc = corrected;
	// 	ltFam.addMember(corrected);
	// 	unsigned int count = 0;
	// 	while(doContinuation){
	// 		printf("Orbit %04d\n", count);
	// 		ic = corrected.getStateByIx(0);
	// 		ic[0] += xStep;

	// 		// Fix the initial x-coordinate
	// 		initStateData[0] = ic[0];
	// 		fixInitState.setData(initStateData);

	// 		nextArc.clearAllConstraints();
	// 		nextArc.addConstraint(periodicityCon);
	// 		nextArc.addConstraint(fixInitState);

	// 		// Set the initial state
	// 		nextArc.setStateByIx(0, ic);

	// 		corrected.reset();
	// 		try{
	// 			shooter.multShoot(&nextArc, &corrected);
	// 			ltFam.addMember(corrected);

	// 			nextArc = corrected;
	// 			count++;
	// 		}catch(Exception &e){
	// 			doContinuation = false;
	// 		}
	// 	}
	// }

	// ltFam.saveToMat("../../data/families_toBeChecked/ltFam_constF_npc.mat");

	// for(unsigned int i = 0; i < loadedLaws.size(); i++){
	// 	delete loadedLaws[i];
	// 	loadedLaws[i] = nullptr;
	// }
}
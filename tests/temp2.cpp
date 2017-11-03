#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

int main() {
	std::vector<double> emL1Lyap_ic {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
	std::vector<double> emL3Lyap_ic {-0.628097117249632, 0, 0, 0, -0.87229410151913, 0, 1};	// EM L3
	double emL1Lyap_T = 3.02797;	// EM L1 Period
	double emL3Lyap_T = 6.2238;
	double alpha = 1.343, beta = 0.0955296;

	SysData_cr3bp_lt sys("earth", "moon", 500);
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::Law_tp::GENERAL_CONST_F, 9e-3, 1500);
	std::vector<double> thrustAngles {alpha, beta};
	
	Arcset_cr3bp_lt halfLyapNodeset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, thrustAngles, 0, emL1Lyap_T/2, 2, &halfLyapNodeset, &law);

	halfLyapNodeset.saveToMat("data/temp.mat");

	// Add control continuity constraint for full-rank Jacobian (check ALL the available partials)
	std::vector<double> conData(law.getNumStates(), 1);
	Constraint ctrlCon(Constraint_tp::CONT_CTRL, 0, conData);

	std::vector<double> stateConData {0.78, 0, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 1, stateConData);

	halfLyapNodeset.addConstraint(ctrlCon);
	halfLyapNodeset.addConstraint(stateCon);

	MultShootEngine corrector;
	// corrector.setVerbosity(Verbosity_tp::NO_MSG);
	// corrector.setVerbosity(Verbosity_tp::SOME_MSG);
	corrector.setTOFType(MSTOF_tp::VAR_FREE);
	corrector.setTol(1e-10);
	
	MultShootEngine::finiteDiff_checkMultShoot(&halfLyapNodeset, corrector, Verbosity_tp::SOME_MSG);
	corrector.multShoot(&halfLyapNodeset, &correctedSet);

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
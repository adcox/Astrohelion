#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

std::vector<double> emL1Lyap_ic {0.887415132364297, 0, 0, 0, -0.332866299501083, 0, 1};	// EM L1
std::vector<double> emL3Lyap_ic {-0.628097117249632, 0, 0, 0, -0.87229410151913, 0, 1};	// EM L3
double emL1Lyap_T = 3.02797;	// EM L1 Period
double emL3Lyap_T = 6.2238;

int main(int argc, char **argv) {
	(void) argc;
	(void) argv;

	MSTOF_tp tofType = MSTOF_tp::VAR_FIXSIGN;
	ControlLaw_cr3bp_lt::Law_tp lawType = ControlLaw_cr3bp_lt::Law_tp::VAR_F_PRO_VEL;

	SysData_cr3bp_lt sys("earth", "moon", 14);
	std::vector<double> ltParams {1500}, ctrlState{1e-2};
	ControlLaw_cr3bp_lt law(lawType, ltParams);

	printf("varF test_rmState, tofType = %d, lawType = %d\n", to_underlying(tofType), to_underlying(lawType));

	Arcset_cr3bp_lt halfLyapArcset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL1Lyap_ic, ctrlState, 0, emL1Lyap_T/2.0, 3, &halfLyapArcset, &law);

	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::SOME_MSG);
	corrector.setTOFType(tofType);
	corrector.setTol(1e-11);
	corrector.setFullFinalProp(false);
	// corrector.setVerbosity(Verbosity_tp::SOME_MSG);

	double stateConData[] = {NAN, 0, NAN, NAN, NAN, NAN, NAN};
	Constraint stateCon(Constraint_tp::STATE, 2, stateConData, 7);
	halfLyapArcset.addConstraint(stateCon);

	Constraint rmState(Constraint_tp::RM_STATE, 0, nullptr, 0);
	halfLyapArcset.addConstraint(rmState);

	// Remove final ctrl state; no constraints on final ctrl state, results in column of zeros
	Constraint rmEndCtrl(Constraint_tp::RM_CTRL, halfLyapArcset.getNodeByIx(-1).getID(), nullptr, 0);
	halfLyapArcset.addConstraint(rmEndCtrl);

	MultShootEngine::finiteDiff_checkMultShoot(&halfLyapArcset, corrector, Verbosity_tp::SOME_MSG);
	halfLyapArcset.saveToMat("data/lt_initGuess.mat");
	corrector.multShoot(&halfLyapArcset, &correctedSet);

	// std::vector<double> finalState = correctedSet.getState(stateCon.getID());
	// stateDiffBelowTol(finalState, stateConData, 1e-12);

	// correctedSet.saveToMat("data/lt_correctedSet.mat");
}
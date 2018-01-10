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

	MSTOF_tp tofType = MSTOF_tp::VAR_FREE;
	ControlLaw_cr3bp_lt::Law_tp lawType = ControlLaw_cr3bp_lt::Law_tp::VAR_F_ANTI_VEL;

	SysData_cr3bp_lt sys("earth", "moon", 14);
	std::vector<double> ltParams {3000}, ctrlState {0.3};

	printf("varF test_continuity, tofType = %d, lawType = %d\n", to_underlying(tofType), to_underlying(lawType));

	ControlLaw_cr3bp_lt law(lawType, ltParams);

	Arcset_cr3bp_lt halfLyapArcset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL3Lyap_ic, ctrlState, 0, emL3Lyap_T, 2, &halfLyapArcset, &law);

	// Remove final ctrl state; no constraints on final ctrl state, results in column of zeros
	Constraint rmEndCtrl(Constraint_tp::RM_CTRL, halfLyapArcset.getNodeByIx(-1).getID(), nullptr, 0);
	halfLyapArcset.addConstraint(rmEndCtrl);

	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::SOME_MSG);
	corrector.setTOFType(tofType);

	MultShootEngine::finiteDiff_checkMultShoot(&halfLyapArcset, corrector, Verbosity_tp::DEBUG, true);

	try{
		corrector.multShoot(&halfLyapArcset, &correctedSet);
	}catch(Exception &e){
		printErr("Error (%s, %s): %s\n", ControlLaw_cr3bp_lt::lawTypeToString(lawType).c_str(),
			MSTOF_tp_cStr(tofType), e.what());
		char filename[128];
		sprintf(filename, "%s_%s_input.mat", ControlLaw_cr3bp_lt::lawTypeToString(lawType).c_str(),
			MSTOF_tp_cStr(tofType));
		halfLyapArcset.saveToMat(filename);

		sprintf(filename, "%s_%s_corrected.mat", ControlLaw_cr3bp_lt::lawTypeToString(lawType).c_str(),
			MSTOF_tp_cStr(tofType));
		correctedSet.saveToMat(filename);
	}
}
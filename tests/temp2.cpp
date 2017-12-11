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

	unsigned int lawType = ControlLaw_cr3bp_lt::Law_tp::VAR_F_CONST_C_2D_LEFT;
	MSTOF_tp tofType = MSTOF_tp::VAR_FREE;

	SysData_cr3bp_lt sys("earth", "moon", 14);
	std::vector<double> ltParams {3000}, ctrlState {0.3};

	ControlLaw_cr3bp_lt law(lawType, ltParams);

	Arcset_cr3bp_lt halfLyapArcset(&sys), correctedSet(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(emL3Lyap_ic, ctrlState, 0, emL3Lyap_T, 2, &halfLyapArcset, &law);

	MultShootEngine corrector;
	corrector.setVerbosity(Verbosity_tp::NO_MSG);
	corrector.setTOFType(tofType);

	MultShootEngine::finiteDiff_checkMultShoot(&halfLyapArcset, corrector, Verbosity_tp::SOME_MSG, true);

	// try{
	// 	corrector.multShoot(&halfLyapArcset, &correctedSet);
	// }catch(Exception &e){
	// 	printErr("Error (%s, %s): %s\n", ControlLaw_cr3bp_lt::lawTypeToString(lawType).c_str(),
	// 		MSTOF_tp_cStr(tofType), e.what());
	// 	char filename[128];
	// 	sprintf(filename, "%s_%s_input.mat", ControlLaw_cr3bp_lt::lawTypeToString(lawType).c_str(),
	// 		MSTOF_tp_cStr(tofType));
	// 	halfLyapArcset.saveToMat(filename);

	// 	sprintf(filename, "%s_%s_corrected.mat", ControlLaw_cr3bp_lt::lawTypeToString(lawType).c_str(),
	// 		MSTOF_tp_cStr(tofType));
	// 	correctedSet.saveToMat(filename);
	// }

	//--------------------------------------------

	// std::vector<double> params {3000};
	// SysData_cr3bp_lt sys("earth", "moon", 1000);
	// ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::Law_tp::VAR_F_CONST_C_2D_LEFT, params);
	
	// std::vector<double> thrustVals {1e-2, 5e-2, 2e-1, 1e-2};
	// double tof = 1;

	// SimEngine sim;
	// sim.setVerbosity(Verbosity_tp::DEBUG);
	// Arcset_cr3bp_lt arc(&sys), fullArc(&sys);
	// std::vector<double> ctrl_state {thrustVals[0]};

	// for(unsigned int i = 0; i < thrustVals.size(); i++){
	// 	arc.reset();
	// 	sim.runSim(ic, ctrl_state, 0, tof, &arc, &law);

	// 	if(i == 0)
	// 		fullArc = arc;
	// 	else
	// 		fullArc += arc;

	// 	ic = arc.getStateByIx(-1);
	// 	ctrl_state[0] = thrustVals[i];
	// }

	// fullArc.saveToMat("data/multi_thrust_arc.mat");

}
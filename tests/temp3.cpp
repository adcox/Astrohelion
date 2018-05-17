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

int main(){

	const char* famFile = "../../data/families/cr3bp-lt_earth-moon/"
		"L4_Lyap_f7.0e-02_a060.00_law2112.mat";

	SysData_cr3bp_lt ltSys(famFile);
	Family_PO_cr3bp_lt fam(&ltSys);

	int numMan = 200;
	double tof = 18*PI;


	std::vector<ControlLaw*> loadedLaws {};
	fam.readFromMat(famFile, loadedLaws);

	double H_lt = -1.561;
	std::vector<Arcset_periodic> matches = fam.getMemberByH_lt(H_lt);

	if(matches.empty()){
		printErr("Did not find any matches with H_lt = %.4f\n", H_lt*180/PI);
		freeMem(loadedLaws);
		return EXIT_SUCCESS;
	}

	assert(matches[0].getCtrlLawByIx(0) == loadedLaws[0]);
	ControlLaw_cr3bp_lt *pLaw = static_cast<ControlLaw_cr3bp_lt*>(loadedLaws[0]);

	Arcset_cr3bp_lt famArc = static_cast<Arcset_cr3bp_lt>(matches[0]);
	famArc.saveToMat("L4_ltpo.mat", Save_tp::SAVE_CURVE);

	std::vector<double> ctrl0 = famArc.getExtraParamVecByIx(0, PARAMKEY_CTRL);

	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	Arcset_cr3bp_lt singleSegArc(&ltSys);

	sim.runSim_manyNodes(famArc.getStateByIx(0), ctrl0, 0, famArc.getTotalTOF(), 
		numMan, &singleSegArc, famArc.getCtrlLawByIx(0));
	singleSegArc.setSTMs_sequence();


	ManifoldEngine man;
	man.setStepOffDist(50);
	std::vector<Arcset_cr3bp_lt> ics_u = man.computeSetFromLTPeriodic(\
		Manifold_tp::MAN_U, &singleSegArc, pLaw, numMan, 0.1);
	std::vector<Arcset_cr3bp_lt> ics_s = man.computeSetFromLTPeriodic(\
		Manifold_tp::MAN_S, &singleSegArc, pLaw, numMan, 0.1);

	Event evt_y0(Event_tp::XZ_PLANE, 0, true);

	std::vector<double> dist_data {0, 1.5};
	Event evt_rBig(Event_tp::DIST, 0, true, dist_data);

	std::vector<double> dist2_data {0, 0.5};
	Event evt_rSmall(Event_tp::DIST, 0, true, dist2_data);

	sim.addEvent(evt_y0);
	sim.addEvent(evt_rBig);
	// sim.addEvent(evt_rSmall);

	Family_PO_cr3bp_lt manifolds_u(&ltSys), manifolds_s(&ltSys);
	Arcset_cr3bp_lt temp(&ltSys);
	for(unsigned int m = 0; m < ics_u.size(); m++){
		printf("Unstable Manifold %03u\n", m);
		temp.reset();
		sim.runSim(ics_u[m].getStateByIx(0),
			ics_u[m].getExtraParamVecByIx(0, PARAMKEY_CTRL), 0, tof, &temp, pLaw);

		// if(temp.getNodeRefByIx(-1).getTriggerEvent() != Event_tp::SIM_TOF)
			manifolds_u.addMember(temp);
	}

	sim.setRevTime(true);
	for(unsigned int m = 0; m < ics_s.size(); m++){
		printf("Stable Manifold %03u\n", m);
		temp.reset();
		sim.runSim(ics_s[m].getStateByIx(0),
			ics_s[m].getExtraParamVecByIx(0, PARAMKEY_CTRL), 0, tof, &temp, pLaw);

		// if(temp.getNodeRefByIx(-1).getTriggerEvent() != Event_tp::SIM_TOF)
			manifolds_s.addMember(temp);
	}

	manifolds_u.saveToMat("L4_unstableManifolds.mat", Save_tp::SAVE_CURVE);
	manifolds_s.saveToMat("L4_stableManifolds.mat", Save_tp::SAVE_CURVE);

	freeMem(loadedLaws);
	return EXIT_SUCCESS;
}
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
		"L1_Lyap_f7.0e-02_Hlt-1.600_law2112.mat";

	SysData_cr3bp_lt ltSys(famFile);
	Family_PO_cr3bp_lt fam(&ltSys);

	int numMan = 200;
	double tof = 3*PI;


	std::vector<ControlLaw*> loadedLaws {};
	fam.readFromMat(famFile, loadedLaws);

	double alpha = 70*PI/180;
	std::vector<Arcset_periodic> matches = fam.getMemberBy2DThrustAngle(alpha);

	if(matches.empty()){
		printErr("Did not find any matches with alpha = %.2f deg\n", alpha*180/PI);
		freeMem(loadedLaws);
		return EXIT_SUCCESS;
	}

	assert(matches[0].getCtrlLawByIx(0) == loadedLaws[0]);
	ControlLaw_cr3bp_lt *pLaw = static_cast<ControlLaw_cr3bp_lt*>(loadedLaws[0]);

	Arcset_cr3bp_lt famArc = static_cast<Arcset_cr3bp_lt>(matches[0]);
	std::vector<double> ctrl0 = famArc.getExtraParamVecByIx(0, PARAMKEY_CTRL);

	if(std::abs(ctrl0[0] - alpha) > 1e-10){
		printErr("Loaded family member has the wrong angle...\n");
		freeMem(loadedLaws);
		return EXIT_SUCCESS;
	}

	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	Arcset_cr3bp_lt singleSegArc(&ltSys);

	sim.runSim_manyNodes(famArc.getStateByIx(0), ctrl0, 0, famArc.getTotalTOF(), 
		numMan, &singleSegArc, famArc.getCtrlLawByIx(0));
	singleSegArc.setSTMs_sequence();


	ManifoldEngine man;
	std::vector<Arcset_cr3bp_lt> ics_u = man.computeSetFromLTPeriodic(\
		Manifold_tp::MAN_U, &singleSegArc, pLaw, numMan, 0.1);
	std::vector<Arcset_cr3bp_lt> ics_s = man.computeSetFromLTPeriodic(\
		Manifold_tp::MAN_S, &singleSegArc, pLaw, numMan, 0.1);

	std::vector<double> xMoon_data {1-ltSys.getMu()};
	Event evt_xMoon(Event_tp::YZ_PLANE, 0, true, xMoon_data);

	std::vector<double> xInt_data {0.7};
	Event evt_xInterior(Event_tp::YZ_PLANE, 0, true, xInt_data);

	sim.addEvent(evt_xMoon);
	sim.addEvent(evt_xInterior);
	
	Family_PO_cr3bp_lt manifolds_u(&ltSys), manifolds_s(&ltSys);
	Arcset_cr3bp_lt temp(&ltSys);
	for(unsigned int m = 0; m < ics_u.size(); m++){
		printf("Unstable Manifold %03u\n", m);
		temp.reset();
		sim.runSim(ics_u[m].getStateByIx(0),
			ics_u[m].getExtraParamVecByIx(0, PARAMKEY_CTRL), 0, tof, &temp, pLaw);

		if(temp.getNodeRefByIx(-1).getTriggerEvent() == evt_xMoon.getType())
			manifolds_u.addMember(temp);
	}

	sim.setRevTime(true);
	for(unsigned int m = 0; m < ics_s.size(); m++){
		printf("Stable Manifold %03u\n", m);
		temp.reset();
		sim.runSim(ics_s[m].getStateByIx(0),
			ics_s[m].getExtraParamVecByIx(0, PARAMKEY_CTRL), 0, tof, &temp, pLaw);

		if(temp.getNodeRefByIx(-1).getTriggerEvent() == evt_xMoon.getType())
			manifolds_s.addMember(temp);
	}

	manifolds_u.saveToMat("unstableManifolds.mat", Save_tp::SAVE_CURVE);
	manifolds_s.saveToMat("stableManifolds.mat", Save_tp::SAVE_CURVE);

	freeMem(loadedLaws);
	return EXIT_SUCCESS;
}
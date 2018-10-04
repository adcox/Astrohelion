/**
 * @file ltpo_manifolds.cpp
 * This script computes stable and unstable manifolds of a low-thrust periodic
 * orbit.
 * 
 * Usage:
 * 
 * 	ltpo_manifolds <fam_file> <orb_ix> <numMan> <tof> 
 * 		Compute the manifolds of an orbit from the family specified by fam_file
 * 		at index orb_ix (orb_ix = 0 is the first family member). To specify the 
 * 		number of manifolds, append a positive integer as numMan and a positive 
 * 		double as tof
 * 	
 * 
 * @author Andrew Cox
 * @version June 14, 2018
 */
#include "AllIncludes.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;
using ltlaw = ControlLaw_cr3bp_lt;

double stepDist = 50;	// kilometers

// Function declarations
void freeMem(std::vector<ControlLaw*>&);

// Function definitions
void freeMem(std::vector<ControlLaw*>& laws){
	for(ControlLaw* pLaw : laws){
		delete pLaw;
		pLaw = nullptr;
	}
}//====================================================

void mexErrMsgTxt(const char* msg){
	printErr(msg);
}


int main(){
	
	Manifold_tp tp = Manifold_tp::MAN_ALL;
	int numMan = 50;
	double tof = 4*PI;
	std::vector<Event> events {};

	const char* fFile = "../../data/families/cr3bp-lt_earth-moon/E2_f7.0e-02_alph000.00_law2112.mat";
	SysData_cr3bp_lt sys(fFile);
	std::vector<ControlLaw*> laws {};
	Family_PO_cr3bp_lt fam(&sys);
	fam.readFromMat(fFile, laws);
	
	std::vector<Arcset_periodic> matches = fam.getMemberByH_lt(-1.55);
	Arcset_cr3bp_lt arc(matches[0]);

	try{
		// Get the control law, or create an object that does not apply control
		ControlLaw_cr3bp_lt *pLaw = nullptr;
		std::vector<double> ctrl0 {};
		if(!laws.empty()){
			pLaw = static_cast<ControlLaw_cr3bp_lt*>(arc.getCtrlLawByIx(0));
			ctrl0 = arc.getExtraParamVecByIx(0, PARAMKEY_CTRL);
		}else{
			laws.push_back(new ControlLaw_cr3bp_lt());
			pLaw = static_cast<ControlLaw_cr3bp_lt*>(laws[0]);
		}

		// Space nodes around the arc by running a quick propagation
		SimEngine sim;
		// sim.setMakeDefaultEvents(false);
		sim.setVerbosity(Verbosity_tp::NO_MSG);
		Arcset_cr3bp_lt multiNodeArc(&sys);
		sim.runSim_manyNodes(arc.getStateByIx(0), ctrl0, 0, arc.getTotalTOF(),
			numMan, &multiNodeArc, pLaw);


		// Compute the manifold initial conditions
		multiNodeArc.setSTMs_cumulative();
		ManifoldEngine man;
		man.setStepOffDist(stepDist);
		std::vector<Arcset_cr3bp_lt> ics;
		try{
			ics = man.computeSetFromLTPeriodic(tp, &multiNodeArc, pLaw, numMan, 0.1);
		}catch(const Exception &e){
			freeMem(laws);
			char msg[256];
			snprintf(msg, 256, "%s\n", e.what());
			mexErrMsgTxt(msg);
		}

		// Propagate the manifold initial conditions 
		for(unsigned int e = 0; e < events.size(); e++){
			sim.addEvent(events[e]);
		}

		Arcset_cr3bp_lt temp(&sys);
		std::vector<Arcset_cr3bp_lt> manifolds(ics.size(), temp);
		for(unsigned int i = 0; i < ics.size(); i++){
			astrohelion::printf("Manifold %03u/%03zu\n", i, ics.size());
			
			if(ics[i].getTotalTOF() < 0)
				sim.setRevTime(true);
			else
				sim.setRevTime(false);

			sim.runSim(ics[i].getStateByIx(0), 
				ics[i].getExtraParamVecByIx(0, PARAMKEY_CTRL), 0, tof,
				&(manifolds[i]), pLaw);
		}

		// -------------------------------------------------------------------------
		// Format the outputs
		// -------------------------------------------------------------------------
		// const char* fieldnames[] = {VARNAME_LINKTABLE, VARNAME_NODESTATE,
		// 	VARNAME_NODETIME, VARNAME_NODECTRL, VARNAME_SEGSTATE, VARNAME_SEGTIME,
		// 	VARNAME_SEGCTRL, VARNAME_TOF, VARNAME_STM,
		// 	VARNAME_P1, VARNAME_P2, VARNAME_MU, VARNAME_M0};
		// const unsigned int nFields = 13;

		// plhs[0] = mxCreateStructMatrix(manifolds.size(), 1, nFields, fieldnames);

		// for(unsigned int m = 0; m < manifolds.size(); m++){
		// 	astroHelper::arcset2MexData(&(manifolds[m]), fieldnames, nFields, 
		// 		plhs[0], m);
		// }
	}catch(const std::exception &e){
		freeMem(laws);
		char msg[256];
		snprintf(msg, 256, "%s\n", e.what());
		mexErrMsgTxt(msg);
	}

	return EXIT_SUCCESS;
}
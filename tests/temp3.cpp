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


int main(){
	
	// --------------------------------------------------------------
	// Load the families, control laws, and system data
	// --------------------------------------------------------------
	const char ltpo_L4_famFile[] = "../../data/families/cr3bp-lt_earth-moon/"
		"L4_Lyap_f7.0e-02_a060.00_law2112.mat";
	const char L4_spo_famFile[] = "../../data/families/cr3bp_earth-moon/L4_SPO.mat";
	// const char lyap_famFile[] = "../../data/families/cr3bp_earth-moon/L1_Lyap.mat";
	const char lyap_famFile[] = "../../data/families/cr3bp_earth-moon/DRO.mat";

	std::vector<ControlLaw*> ltpo_laws {}, spo_laws {}, lyap_laws {};

	SysData_cr3bp_lt ltSys(ltpo_L4_famFile);
	Family_PO_cr3bp_lt lt_L4_fam(&ltSys);
	lt_L4_fam.readFromMat(ltpo_L4_famFile, ltpo_laws);
	
	SysData_cr3bp natSys(L4_spo_famFile);
	Family_PO_cr3bp spo_fam(&natSys), lyap_fam(&natSys);

	spo_fam.readFromMat(L4_spo_famFile, spo_laws);
	lyap_fam.readFromMat(lyap_famFile, lyap_laws);

	assert(spo_laws.size() == 0);
	assert(lyap_laws.size() == 0);

	// --------------------------------------------------------------
	// Get a LTPO at a specific low-thrust Hamiltonian value
	// --------------------------------------------------------------
	double Hlt = -1.561;
	std::vector<Arcset_periodic> matches = lt_L4_fam.getMemberByH_lt(Hlt);
	if(matches.size() == 0){
		printErr("Could not find a LTPO with H = %f\n", Hlt);
		freeMem(ltpo_laws);
		return EXIT_SUCCESS;
	}

	Arcset_cr3bp_lt ltpo_frame = matches[0];
	assert(ltpo_frame.getCtrlLawByIx(0) == ltpo_laws[0]);
	ControlLaw_cr3bp_lt *pLaw = static_cast<ControlLaw_cr3bp_lt*>(ltpo_laws[0]);

	// --------------------------------------------------------------
	// Compute unstable manifolds of the LTPO_frame
	// --------------------------------------------------------------
	int numMan = 200;
	double tof = 20*PI;
	std::vector<double> ctrl0 = ltpo_frame.getExtraParamVecByIx(0, PARAMKEY_CTRL);

	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	Arcset_cr3bp_lt ltpo(&ltSys);

	sim.runSim_manyNodes(ltpo_frame.getStateByIx(0), ctrl0, 0, ltpo_frame.getTotalTOF(), 
		numMan, &ltpo, ltpo_frame.getCtrlLawByIx(0));
	ltpo.setSTMs_sequence();

	// Compute the manifold initial conditions (no current support for 
	// propagating with low-thrust in the ManifoldEngine)
	ManifoldEngine man;
	man.setStepOffDist(stepDist);
	std::vector<Arcset_cr3bp_lt> ics_u;
	try{
		ics_u = man.computeSetFromLTPeriodic(Manifold_tp::MAN_U, &ltpo,
			pLaw, numMan, 0.1);
	}catch(const Exception &e){
		freeMem(ltpo_laws);
		printErr("%s\n", e.what());
		return EXIT_SUCCESS;
	}

	// Add event to SimEngine to stop propagation at the y-axis
	Event evt_yAx(Event_tp::XZ_PLANE, 0, true);
	sim.addEvent(evt_yAx);

	// Propagate the unstable manifold ICs in forward time
	Arcset_cr3bp_lt temp(&ltSys), manifold(&ltSys);
	bool foundArc = false;
	for(unsigned int m = 0; m < ics_u.size(); m++){
		printf("Propagating unstable manifold %03u...\n", m);
		temp.reset();
		sim.runSim(ics_u[m].getStateByIx(0),
			ics_u[m].getExtraParamVecByIx(0, PARAMKEY_CTRL), 0, tof, &temp,
			pLaw);

		std::vector<double> qf = temp.getStateByIx(-1);
		if(qf[0] > 0.9){
			// Simulate with variable thrust magnitude
			unsigned int id = ltlaw::GENERAL | ltlaw::VAR_F_BND | ltlaw::CONST_M;
			ctrl0.push_back(0.4*PI);	// max thrust
			pLaw->setType(id);

			sim.runSim_manyNodes(ics_u[m].getStateByIx(0), ctrl0, 0, 
				temp.getTotalTOF(), 15, &manifold, pLaw);
			foundArc = true;
			break;
		}
	}

	if(!foundArc){
		printErr("Could not find a manifold with the specified final state\n");
		freeMem(ltpo_laws);
	}

	// --------------------------------------------------------------
	// Get a Lyap and SPO with an energies similar to the manifold ends
	// --------------------------------------------------------------
	double C_manEnd = manifold.getJacobiByIx(-1);
	double C_manStart = manifold.getJacobiByIx(0);

	if(C_manStart > 2.985){ C_manStart = 2.985; }	// ensure SPO isn't tiny

	matches = spo_fam.getMemberByJacobi(C_manStart);
	if(matches.size() == 0){
		printErr("Could not find an SPO with C = %f\n", C_manStart);
		freeMem(ltpo_laws);
		return EXIT_SUCCESS;
	}
	Arcset_cr3bp spo = matches[0];

	matches = lyap_fam.getMemberByJacobi(C_manEnd);
	if(matches.size() == 0){
		printErr("Could not find a Lyapunov with C = %f\n", C_manEnd);
		freeMem(ltpo_laws);
		return EXIT_SUCCESS;
	}
	Arcset_cr3bp lyap = matches[0];

	spo.saveToMat("transfer_SPO.mat");
	lyap.saveToMat("transfer_Lyap.mat");
	manifold.saveToMat("transfer_Manifold.mat");

	// --------------------------------------------------------------
	// Construct an end-to-end transfer
	// --------------------------------------------------------------
	Arcset_cr3bp_lt spo_lt(&ltSys), lyap_lt(&ltSys), transfer(&ltSys);

	// unsigned int lawID = ltlaw::GENERAL | ltlaw::VAR_F_BND | ltlaw::CONST_M;
	// ControlLaw_cr3bp_lt law_varThrust(lawID, std::vector<double> {7e-2});
	// ctrl0 = std::vector<double> {PI/3.0, 0, -0.25*PI};	// Angle = 0, psi = -pi/2 for f = 0

	unsigned int lawID = ltlaw::NO_CTRL;
	ControlLaw_cr3bp_lt law_noThrust(lawID, std::vector<double> {});
	ctrl0 = std::vector<double> {};

	sim.clearEvents();
	
	std::vector<double> q0 = spo.getStateByIx(0);
	q0.push_back(1);
	sim.runSim_manyNodes(q0, ctrl0, 0, spo.getTotalTOF(), 
		spo.getNumNodes(), &spo_lt, &law_noThrust);

	q0 = lyap.getStateByIx(0);
	q0.push_back(1);
	sim.runSim_manyNodes(q0, ctrl0, 0, lyap.getTotalTOF(),
		lyap.getNumNodes(), &lyap_lt, &law_noThrust);

	transfer = manifold;
	int id0 = transfer.getNodeByIx(0).getID(), 
		idf = transfer.getNodeByIx(-1).getID();
	transfer.appendSetAtNode(&spo_lt, id0, spo_lt.getNodeByIx(-1).getID(), 0);
	transfer.appendSetAtNode(&lyap_lt, idf, lyap_lt.getNodeByIx(0).getID(), 0);

	transfer.sortChrono();
	transfer.saveToMat("transfer_raw.mat");

	// --------------------------------------------------------------
	// Attempt corrections - continuity
	// --------------------------------------------------------------
	MultShootEngine shooter;
	Arcset_cr3bp_lt corrected(&ltSys);
	
	// Constrain second and second to last states on natural orbits
	Constraint fixStart(Constraint_tp::STATE, transfer.getNodeByIx(1).getID(),
		transfer.getStateByIx(1));
	Constraint fixEnd(Constraint_tp::STATE, transfer.getNodeByIx(-2).getID(),
		transfer.getStateByIx(-2));

	transfer.addConstraint(fixStart);
	transfer.addConstraint(fixEnd);

	printf("Correcting transfer...\n");
	shooter.setDoLineSearch(true);
	shooter.setMaxIts(200);
	try{
		shooter.multShoot(&transfer, &corrected);
	}catch(Exception &e){
		printErr("%s\n", e.what());
		freeMem(ltpo_laws);
		return EXIT_SUCCESS;
	}

	corrected.saveToMat("transfer_cont.mat");

	freeMem(ltpo_laws);
	return EXIT_SUCCESS;
}
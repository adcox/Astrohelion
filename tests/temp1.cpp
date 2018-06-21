#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;
using ltlaw = astrohelion::ControlLaw_cr3bp_lt;

void freeMem(std::vector<ControlLaw*>& laws){
	for(auto p : laws){
		delete p;
		p = nullptr;
	}
}//====================================================

int main(int argc, char** argv){
	
	char filepath[] = "../../data/families/cr3bp-lt_earth-moon";
	char filename[] = "L4_Lyap_f7.0e-02_Hlt-1.552_law2112.mat";

	char fullpath[256];
	sprintf(fullpath, "%s/%s", filepath, filename);

	SysData_cr3bp_lt sys(fullpath);
	Family_PO_cr3bp_lt fam(&sys);

	std::vector<ControlLaw*> loadedLaws {};
	fam.readFromMat(fullpath, loadedLaws);

	// Get orbits associated with these discrete angles
	std::vector<double> alpha_vals {4*PI/12, 5*PI/12, 6*PI/12, 7*PI/12, 8*PI/12,
		9*PI/12, 10*PI/12, 11*PI/12,  11.9*PI/12};
	std::vector<Arcset_cr3bp_lt> fullArcs {};

	for(double alpha : alpha_vals){
		std::vector<Arcset_periodic> matches = fam.getMemberBy2DThrustAngle(alpha);

		if(!matches.empty()){
			fullArcs.push_back(static_cast<Arcset_cr3bp_lt>(matches[0]));
			// char fn[64];
			// sprintf(fn, "L4_f7.0e-02_Hlt-1.552_a%06.2f.mat", alpha*180/PI);
			// matches[0].saveToMat(fn);
		}else{
			printErr("Could not load family member with alpha = %.2f deg\n",
				alpha*180/PI);
			freeMem(loadedLaws);
			return EXIT_SUCCESS;
		}
	}

	// Try multiple shooting without any adjustments
	Arcset_cr3bp_lt chain = fullArcs[0];
	for(unsigned int i = 1; i < fullArcs.size(); i++){
		chain.appendSetAtNode(&(fullArcs[i]), chain.getNodeByIx(-1).getID(),
			fullArcs[i].getNodeByIx(0).getID(), 1.5, 
			fullArcs[i].getCtrlLawByIx(-1));
	}

	chain.clearAllConstraints();
	chain.sortChrono();

	// Constrain the first and last states
	std::vector<double> q0 = chain.getStateByIx(0);
	std::vector<double> qf = chain.getStateByIx(-1);
	
	Constraint fix_q0(Constraint_tp::STATE, chain.getNodeRefByIx(0).getID(), q0);
	Constraint fix_qf(Constraint_tp::STATE, chain.getNodeRefByIx(-1).getID(), qf);

	chain.addConstraint(fix_q0);
	chain.addConstraint(fix_qf);

	chain.saveToMat("L4_chain.mat");

	// Adjust the control law to use variable mass
	unsigned int lawID = ltlaw::GENERAL | ltlaw::CSI_VAR_M | ltlaw::CONST_F;
	ltlaw *pLaw = static_cast<ltlaw *>(chain.getCtrlLawByIx(0));
	std::vector<double> params = pLaw->getParams();
	params.push_back(1500);		// Isp
	// *pLaw = ltlaw(lawID, params);

	MultShootEngine shooter;
	shooter.setMaxIts(200);
	// shooter.setDoLineSearch(true);
	// shooter.setTOFType(MSTOF_tp::VAR_FIXSIGN);

	// Converge with fixed mass
	printf("Converging with fixed mass\n");
	Arcset_cr3bp_lt converged(&sys);
	try{
		shooter.multShoot(&chain, &converged);
		converged.saveToMat("L4_chain_converged.mat");
	}catch(const DivergeException &e){
		printErr("Multiple shooting diverged\n");
		freeMem(loadedLaws);
		return EXIT_SUCCESS;
	}catch(const std::exception &e){
		freeMem(loadedLaws);
		throw;
	}

	// Converge again with variable mass
	printf("Converging with variable mass\n");
	Arcset_cr3bp_lt conv_varMass(&sys);
	*pLaw = ltlaw(lawID, params);
	shooter.setTOFType(MSTOF_tp::VAR_FIXSIGN);
	shooter.setDoLineSearch(true);
	try{
		shooter.multShoot(&converged, &conv_varMass);
		conv_varMass.saveToMat("L4_chain_converged_varMass.mat");
	}catch(const std::exception &e){
		freeMem(loadedLaws);
		throw;
	}

	freeMem(loadedLaws);
	return EXIT_SUCCESS;
}
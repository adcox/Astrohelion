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
	std::vector<double> alpha_vals {4*PI/12, 5*PI/12, 6*PI/12, 7*PI/12, 8*PI/12};
	std::vector<Arcset_cr3bp_lt> fullArcs {};

	for(double alpha : alpha_vals){
		std::vector<Arcset_periodic> matches = fam.getMemberBy2DThrustAngle(alpha);

		if(!matches.empty()){
			fullArcs.push_back(static_cast<Arcset_cr3bp_lt>(matches[0]));
			char fn[64];
			sprintf(fn, "L4_f7.0e-02_Hlt-1.552_a%06.2f.mat", alpha*180/PI);
			matches[0].saveToMat(fn);
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
	chain.putInChronoOrder();
	chain.saveToMat("L4_chain.mat");

	MultShootEngine shooter;
	shooter.setMaxIts(200);

	Arcset_cr3bp_lt converged(&sys);
	try{
		shooter.multShoot(&chain, &converged);
		converged.saveToMat("L4_chain_converged.mat");
	}catch(const DivergeException &e){
		printErr("Multiple shooting diverged\n");
		freeMem(loadedLaws);
		return EXIT_SUCCESS;
	}catch(const Exception &e){
		freeMem(loadedLaws);
		throw;
	}catch(const std::exception &e){
		freeMem(loadedLaws);
		throw;
	}

	freeMem(loadedLaws);
	return EXIT_SUCCESS;
}
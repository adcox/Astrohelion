//#include "Calculations.hpp"
#include "Common.hpp"
#include "Fam_cr3bp.hpp"
#include "FamMember_cr3bp.hpp"
#include "ManifoldEngine.hpp"
#include "SimEngine.hpp"
#include "SysData_cr3bp.hpp"
#include "Arcset_cr3bp.hpp"

#include <vector>

using namespace astrohelion;

int main(){

	Fam_cr3bp fam("../share/families/EM_L1_Lyap.mat");
	SysData_cr3bp sysData = fam.getSysData();
	std::vector<FamMember_cr3bp> halos = fam.getMemberByJacobi(3.188);

	printf("Found %zu matches!\n", halos.size());

	SimEngine sim;
	Arcset_cr3bp lyap(&sysData);
	sim.runSim(halos[0].getIC(), halos[0].getTOF(), &lyap);
	lyap.saveToMat("data/manifoldTest_periodicOrbit.mat");

	ManifoldEngine engine;
	std::vector<Arcset_cr3bp> manifolds = engine.computeSetFromPeriodic(Manifold_tp::MAN_S_P, &lyap, 100, 2*PI);
	for(unsigned int i = 0; i < manifolds.size(); i++){
		char filename[128];
		sprintf(filename, "data/man%03zu.mat", i);
		manifolds[i].saveToMat(filename);
	}
	return EXIT_SUCCESS;
}
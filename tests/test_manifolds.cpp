#include "Calculations.hpp"
#include "Common.hpp"
#include "Fam_cr3bp.hpp"
#include "FamMember_cr3bp.hpp"
#include "SimEngine.hpp"
#include "SysData_cr3bp.hpp"
#include "Traj_cr3bp.hpp"

#include <vector>

using namespace astrohelion;

int main(){

	Fam_cr3bp fam("../share/families/EM_L1_Lyap.mat");
	SysData_cr3bp sysData = fam.getSysData();
	std::vector<FamMember_cr3bp> halos = fam.getMemberByJacobi(3.188);

	printf("Found %zu matches!\n", halos.size());

	SimEngine sim;
	Traj_cr3bp aHalo(&sysData);
	sim.runSim(halos[0].getIC(), halos[0].getTOF(), &aHalo);
	aHalo.saveToMat("data/halo_JC312.mat");

	std::vector<Traj_cr3bp> manifolds = getManifolds(Manifold_tp::MAN_S_P, &aHalo, 20, 2*PI);
	for(size_t i = 0; i < manifolds.size(); i++){
		char filename[128];
		sprintf(filename, "data/man%03zu.mat", i);
		manifolds[i].saveToMat(filename);
	}
	return EXIT_SUCCESS;
}
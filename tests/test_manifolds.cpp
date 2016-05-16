#include "tpat_calculations.hpp"
#include "tpat_constants.hpp"
#include "tpat_family_cr3bp.hpp"
#include "tpat_family_member_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"

#include <vector>

int main(){

	tpat_family_cr3bp fam("../share/families/EM_L1_Lyap.mat");
	tpat_sys_data_cr3bp sysData = fam.getSysData();
	std::vector<tpat_family_member_cr3bp> halos = fam.getMemberByJacobi(3.188);

	printf("Found %zu matches!\n", halos.size());

	tpat_simulation_engine sim;
	tpat_traj_cr3bp aHalo(&sysData);
	sim.runSim(halos[0].getIC(), halos[0].getTOF(), &aHalo);
	aHalo.saveToMat("data/halo_JC312.mat");

	std::vector<tpat_traj_cr3bp> manifolds = getManifolds(MAN_S_P, &aHalo, 20, 2*PI);
	for(size_t i = 0; i < manifolds.size(); i++){
		char filename[128];
		sprintf(filename, "data/man%03zu.mat", i);
		manifolds[i].saveToMat(filename);
	}
	return EXIT_SUCCESS;
}
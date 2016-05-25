#include "tpat_calculations.hpp"
#include "tpat_constants.hpp"
#include "tpat_fam_cr3bp.hpp"
#include "tpat_famMember_cr3bp.hpp"
#include "tpat_sim_engine.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"

#include <vector>

int main(){

	TPAT_Fam_CR3BP fam("../share/families/EM_L1_Lyap.mat");
	TPAT_Sys_Data_CR3BP sysData = fam.getSysData();
	std::vector<TPAT_FamMember_CR3BP> halos = fam.getMemberByJacobi(3.188);

	printf("Found %zu matches!\n", halos.size());

	TPAT_Sim_Engine sim;
	TPAT_Traj_CR3BP aHalo(&sysData);
	sim.runSim(halos[0].getIC(), halos[0].getTOF(), &aHalo);
	aHalo.saveToMat("data/halo_JC312.mat");

	std::vector<TPAT_Traj_CR3BP> manifolds = getManifolds(TPAT_Manifold_Tp::Man_S_P, &aHalo, 20, 2*PI);
	for(size_t i = 0; i < manifolds.size(); i++){
		char filename[128];
		sprintf(filename, "data/man%03zu.mat", i);
		manifolds[i].saveToMat(filename);
	}
	return EXIT_SUCCESS;
}
#include "tpat_constants.hpp"
#include "tpat_family_cr3bp.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"

#include <cstdlib>
#include <iostream>

int main(){
	printf("Hello world!\n");

	tpat_family_cr3bp EML1Lyap("../share/families_natParam_checked/EM_L1_Lyap.mat");
	tpat_family_cr3bp EMRes_5_2("../share/families_pac_checked/EM_Res_5_2.mat");
	tpat_family_cr3bp EMRes_4_5("../share/families_pac_checked/EM_Res_4_5.mat");

	tpat_sys_data_cr3bp emSys("earth", "moon");
	tpat_sys_data_bcr4bpr bcSys("sun", "earth", "moon");

	double P_syn = 2*PI*bcSys.getCharT()/(sqrt(bcSys.getMu()/pow(bcSys.getCharLRatio(), 3)) - bcSys.getK())/emSys.getCharT();
	printf("Synodic Period = %.6f, nondim\n", P_syn);

	std::vector<tpat_family_member_cr3bp> L1Matches = EML1Lyap.getMemberByTOF(P_syn);
	L1Matches.at(0).toTraj(&emSys).saveToMat("data/LPF_SynResOrbits/L1Lyap.mat");

	EMRes_5_2.setSortType(tpat_family_cr3bp::SORT_JC);
	EMRes_5_2.sortMembers();
	std::vector<tpat_family_member_cr3bp> Res52_Matches = EMRes_5_2.getMemberByTOF(2*P_syn);
	Res52_Matches.at(0).toTraj(&emSys).saveToMat("data/LPF_SynResOrbits/EMRes_5_2.mat");

	EMRes_4_5.setSortType(tpat_family_cr3bp::SORT_JC);
	EMRes_4_5.sortMembers();
	std::vector<tpat_family_member_cr3bp> Res45_Matches = EMRes_4_5.getMemberByTOF(4*P_syn);
	Res45_Matches.at(0).toTraj(&emSys).saveToMat("data/LPF_SynResOrbits/EMRes_4_5.mat");
}
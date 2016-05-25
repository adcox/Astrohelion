#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_traj_cr3bp_ltvp.hpp"
#include "tpat_sys_data_cr3bp_ltvp.hpp"
#include "tpat_sim_engine.hpp"

int main(void){
	// Try to simulate a circular orbit in the CR3BP
	double ic[] = {0.131231781418776, 0, 0, 0, 2.48142854119997, 0};
	double tof = 0.363031858821079;

	TPAT_Sys_Data_CR3BP sys("earth", "moon");
	TPAT_Sim_Engine sim;
	TPAT_Traj_CR3BP traj(&sys);
	sim.runSim(ic, tof, &traj);
	traj.saveToMat("data/circleOrb.mat");

	double mass = 12;	// kg
	TPAT_Sys_Data_CR3BP_LTVP lowThrustSys("earth", "moon", 0.0012, 2500, 12);
	// TPAT_Sim_Engine sim2(lowThrustSys);
	double ic_lt[] = {0.131231781418776, 0, 0, 0, 2.48142854119997, 0, mass/lowThrustSys.getCharM()};
	TPAT_Traj_CR3BP_LTVP lowThrustTraj(&lowThrustSys);
	sim.setVerbose(TPAT_Verbosity_Tp::ALL_MSG);
	sim.runSim(ic_lt, 5*tof, &lowThrustTraj);
	
	lowThrustTraj.saveToMat("data/circleOrb_lowThrust.mat");
}
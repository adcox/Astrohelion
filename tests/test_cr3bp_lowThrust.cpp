#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_traj_cr3bp_ltvp.hpp"
#include "tpat_sys_data_cr3bp_ltvp.hpp"
#include "tpat_simulation_engine.hpp"

int main(void){
	// Try to simulate a circular orbit in the CR3BP
	double ic[] = {0.131231781418776, 0, 0, 0, 2.48142854119997, 0};
	double tof = 0.363031858821079;

	tpat_sys_data_cr3bp sys("earth", "moon");
	tpat_simulation_engine sim(&sys);
	sim.runSim(ic, tof);
	tpat_traj_cr3bp traj = sim.getCR3BP_Traj();
	traj.saveToMat("data/circleOrb.mat");

	double mass = 12;	// kg
	tpat_sys_data_cr3bp_ltvp lowThrustSys("earth", "moon", 0.0012, 2500, 12);
	// tpat_simulation_engine sim2(lowThrustSys);
	double ic_lt[] = {0.131231781418776, 0, 0, 0, 2.48142854119997, 0, mass/lowThrustSys.getCharM()};
	sim.setSysData(&lowThrustSys);
	sim.setVerbose(ALL_MSG);
	sim.runSim(ic_lt, 5*tof);
	
	tpat_traj_cr3bp_ltvp lowThrustTraj = sim.getCR3BP_LTVP_Traj();
	lowThrustTraj.saveToMat("data/circleOrb_lowThrust.mat");
}
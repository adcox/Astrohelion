#include "SysData_cr3bp.hpp"
#include "Traj_cr3bp.hpp"
#include "Traj_cr3bp_ltvp.hpp"
#include "SysData_cr3bp_ltvp.hpp"
#include "SimEngine.hpp"

using namespace astrohelion;

int main(void){
	// Try to simulate a circular orbit in the CR3BP
	double ic[] = {0.131231781418776, 0, 0, 0, 2.48142854119997, 0};
	double tof = 0.363031858821079;

	SysData_cr3bp sys("earth", "moon");
	SimEngine sim;
	Traj_cr3bp traj(&sys);
	sim.runSim(ic, tof, &traj);
	traj.saveToMat("data/circleOrb.mat");

	double mass = 12;	// kg
	SysData_cr3bp_ltvp lowThrustSys("earth", "moon", 0.0012, 2500, 12);
	// SimEngine sim2(lowThrustSys);
	double ic_lt[] = {0.131231781418776, 0, 0, 0, 2.48142854119997, 0, mass/lowThrustSys.getCharM()};
	Traj_cr3bp_ltvp lowThrustTraj(&lowThrustSys);
	sim.setVerbosity(Verbosity_tp::ALL_MSG);
	sim.runSim(ic_lt, 5*tof, &lowThrustTraj);
	
	lowThrustTraj.saveToMat("data/circleOrb_lowThrust.mat");
}
#include "Exceptions.hpp"
#include "SysData_cr3bp.hpp"
#include "Arcset_cr3bp.hpp"
#include "Arcset_cr3bp_lt.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "SimEngine.hpp"

using namespace astrohelion;

int main(void){
	// Try to simulate a circular orbit in the CR3BP
	double ic[] = {0.131231781418776, 0, 0, 0, 2.48142854119997, 0};
	double tof = 0.363031858821079;

	SysData_cr3bp sys("earth", "moon");
	SimEngine sim;
	Arcset_cr3bp traj(&sys);
	sim.runSim(ic, tof, &traj);
	traj.saveToMat("data/circleOrb.mat");

	double mass = 100;		// kg
	double f = 3.6688;		// (1N, nondim in Earth-Moon)
	double Isp = 5000;		// sec
	SysData_cr3bp_lt lowThrustSys("earth", "moon", mass);
	ControlLaw_cr3bp_lt control(ControlLaw_cr3bp_lt::Law_tp::CONST_FC_2D_RIGHT, f, Isp);

	// SimEngine sim2(lowThrustSys);
	double ic_lt[] = {0.131231781418776, 0, 0, 0, 2.48142854119997, 0, 1};
	Arcset_cr3bp_lt lowThrustTraj(&lowThrustSys);
	sim.setVerbosity(Verbosity_tp::ALL_MSG);
	// sim.setSimpleInt(true);
	try{
		sim.runSim(ic_lt, 5*tof, &lowThrustTraj, &control);
	}catch(DivergeException &e){}

	lowThrustTraj.saveToMat("data/circleOrb_lowThrust.mat");
}
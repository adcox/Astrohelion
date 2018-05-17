#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;
using ltlaw = astrohelion::ControlLaw_cr3bp_lt;

int main(int argc, char** argv){
	
	(void) argc;
	(void) argv;

	SysData_cr3bp_lt sys("earth", "moon", 100);
	double theta0 = 0.123, f = 1e-2, Isp = 1500;
	double t0 = -0.1, tof = PI;
	unsigned int numSegs = 2;

	double seg_tof = tof/static_cast<double>(numSegs);

	// Law for inertially-fixed thrust vector
	unsigned int lawID = ltlaw::GEN_INERT | ltlaw::CONST_F | ltlaw::CSI_VAR_M;
	std::vector<double> params {theta0, f, Isp};
	ControlLaw_cr3bp_lt inertLaw(lawID, params);

	// Law for rotating-frame-fixed thrust vector
	lawID = ltlaw::GENERAL | ltlaw::CONST_F | ltlaw::CSI_VAR_M;
	params = std::vector<double> {f, Isp};
	ControlLaw_cr3bp_lt rotLaw(lawID, params);
	
	double a_desired = 0.4;

	std::vector<double> q0 {-0.6281, 0, 0, 0, -0.8723, 0, 1};
	std::vector<double> rotCtrl0 {a_desired, 0};

	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);

	Arcset_cr3bp_lt arc_rotFixed(&sys), arc_inertFixed(&sys);

	// Single propagation for arc with thrust fixed in rotating frame
	sim.runSim(q0, rotCtrl0, t0, tof, &arc_rotFixed, &rotLaw);

	Arcset_cr3bp_lt temp(&sys);
	std::vector<double> qi = q0;
	for(unsigned int s = 0; s < numSegs; s++){
		std::vector<double> inertCtrl0 {a_desired + theta0 + 
			t0 + s*seg_tof + 0.5*seg_tof, 0};

		temp.reset();
		sim.runSim(qi, inertCtrl0, t0 + s*seg_tof, seg_tof, &temp, &inertLaw);
		qi = temp.getStateByIx(-1);
		if(s == 0){
			arc_inertFixed = temp;
		}else{
			arc_inertFixed.appendSetAtNode(&temp, 
				arc_inertFixed.getNodeRefByIx(-1).getID(),
				temp.getNodeRefByIx(0).getID(), 0);
		}
	}

	arc_inertFixed.putInChronoOrder();
	arc_inertFixed.print();

	arc_rotFixed.saveToMat("arc_rotFixed.mat");
	arc_inertFixed.saveToMat("arc_inertFixed.mat");

	return EXIT_SUCCESS;
}
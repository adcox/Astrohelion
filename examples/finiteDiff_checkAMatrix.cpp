#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;
using ltlaw = astrohelion::ControlLaw_cr3bp_lt;

int main(int argc, char** argv){
	// Check the A-matrix for the CR3BP-LT with GEN_INERT
	SysData_cr3bp_lt sys("earth", "moon", 100);
	unsigned int lawID = ltlaw::GEN_INERT | ltlaw::CONST_F | ltlaw::CSI_VAR_M;
	std::vector<double> params {-3.134, 1e-2, 1500};
	ControlLaw_cr3bp_lt law(lawID, params);
	EOM_ParamStruct eomParams(&sys, &law);

	// q = {x, y, z, vx, vy, vz, m, psi, beta}
	double q0[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.123, 0.321};
	double t = 10.4;

	sys.getDynamicsModel()->finiteDiff_checkAMat(q0, t, 1e-7, &eomParams, 
		Verbosity_tp::ALL_MSG);

	return EXIT_SUCCESS;
}
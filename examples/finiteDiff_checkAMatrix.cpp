/**
 * \file finiteDiff_checkAMatrix.cpp
 * 
 * This script uses the finiteDiff_checkAMat() function to check the accuracy
 * of the A matrix for the CR3BP-LT with a specific thrust law
 */

#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include <vector>

using namespace astrohelion;
using ltlaw = astrohelion::ControlLaw_cr3bp_lt;

int main(int argc, char** argv){
	// CR3BP-LT system data; 100kg spacecraft
	SysData_cr3bp_lt sys("earth", "moon", 100);

	// Inertially-pointing, constant thrust magnitude, CSI engine
	unsigned int lawID = ltlaw::GEN_INERT | ltlaw::CONST_F | ltlaw::CSI_VAR_M;

	// params = {theta0, f, Isp}
	std::vector<double> params {-3.134, 1e-2, 1500};

	// Create the control law; store it and system data in EOM_ParamStruct
	ControlLaw_cr3bp_lt law(lawID, params);
	EOM_ParamStruct eomParams(&sys, &law);

	// q = {x, y, z, vx, vy, vz, m, psi, beta}
	double q0[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.123, 0.321};
	double t = 10.4;	// arbitrary but nonzero to test A matrix

	// Check to see if the analytically-computed A matrix matches one computed
	// via central differencing with a step of 1e-7
	bool isOk = sys.getDynamicsModel()->finiteDiff_checkAMat(q0, t, 1e-7, 
		&eomParams, Verbosity_tp::ALL_MSG);

	return EXIT_SUCCESS;
}
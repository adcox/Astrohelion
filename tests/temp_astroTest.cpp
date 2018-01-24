#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

int main(int argc, char** argv){
	(void) argc;
	(void) argv;

	double f = 1e-2;
	SysData_cr3bp_lt sys("earth", "moon", 100);

	std::vector<double> params {f, 3000};
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::Law_tp::CONST_F_C_2D_LEFT, params);

	std::vector<double> q {0.5, 0, 0, 0.5, 0.5, 0, 1};
	double H = DynamicsModel_cr3bp_lt::getHamiltonian(0, &(q[0]), &sys, &law);

	std::cout << "H = " << H << std::endl;
	
	return EXIT_SUCCESS;
}
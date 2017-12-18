#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

int main(int argc, char **argv) {
	(void) argc;
	(void) argv;

	SysData_cr3bp_lt sys("earth", "moon", 14);
	double f = 1e-3;
	std::vector<double> LPts;

	DynamicsModel_cr3bp_lt::getEquilibPt(&sys, 2, f, 1e-8, &LPts);

}
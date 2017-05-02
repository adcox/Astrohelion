#include "AsciiOutput.hpp"
#include "SimEngine.hpp"
#include "Arcset_2bp.hpp"
#include "SysData_2bp.hpp"

#include "BodyData.hpp"

#include <iostream>

using namespace astrohelion;

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

int main(){
	SysData_2bp sys("earth");

	double mu = sys.getMu();
	std::cout << "Mu = " << mu << " km^3/s^2" << std::endl;

	double r0 = 6378 + 500;
	double v0 = sqrt(mu/r0);

	double ic[] = {r0, 0, 0, 0, v0, 0};
	double period = 2*PI*sqrt(r0*r0*r0/mu);

	Arcset_2bp circOrb(&sys);
	SimEngine sim;
	sim.runSim(ic, period, &circOrb);

	circOrb.saveToMat("data/circular2BOrbit.mat");

	return EXIT_SUCCESS;
}
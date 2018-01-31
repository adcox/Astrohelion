#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

int main(int argc, char** argv){
	(void) argc;
	(void) argv;

	std::vector<double> lawParams{7e-2, 3000};

	// Create system object and control object
	const SysData_cr3bp_lt sys("earth", "moon", 1);
	ControlLaw_cr3bp_lt law(ControlLaw_cr3bp_lt::Law_tp::CONST_F_GENERAL, lawParams);

	std::vector<double> q0 {0.7, 0, 0, 0, -0.75, 0, 1};
	std::vector<double> ctrl0 {3.1, 0};

	SimEngine sim;
	Arcset_cr3bp_lt arc(&sys);

	sim.runSim(q0, ctrl0, 0, 10, &arc, &law);

	// Extract the full state history of the arc
	std::vector<double> allStates = arc.getSegRefByIx(0).getStateVector();
	unsigned int w = arc.getSegRefByIx(0).getStateWidth();
	unsigned int h = allStates.size()/w;

	// Compute initial Hamiltonian
	double H0 = DynamicsModel_cr3bp_lt::getHamiltonian(0, &(allStates[0]), &sys, &law);
	
	// Determine the biggest deviation from initial Hamiltonian
	double max_dH = 0, dH = 0;
	for(unsigned int i = 0; i < h; i++){
		printf("state %04u: ", i);
		for(unsigned int j = 0; j < 9; j++)
			printf("%f, ", allStates[w*i + j]);
		printf("\n");
		
		dH = DynamicsModel_cr3bp_lt::getHamiltonian(0, &(allStates[w*i]), &sys, &law) - H0;
		max_dH = std::abs(dH) > std::abs(max_dH) ? dH : max_dH;
	}

	printf("max dH = %e\n", max_dH);

	return EXIT_SUCCESS;
}
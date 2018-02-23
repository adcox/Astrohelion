#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

int main(int argc, char** argv){
	double errTol = 1e-14;
	std::vector<double> ic {0.82575887, 0, 0.08, 0, 0.19369725, 0, 1};
	std::vector<double> ctrl0 {1.25, 0.1};

	SysData_cr3bp_lt ltData("earth", "moon", 14);
	std::vector<double> ltParams {0.3, 1500};
	ControlLaw_cr3bp_lt control(ControlLaw_cr3bp_lt::Law_tp::CONST_F_GENERAL, ltParams);

	Arcset_cr3bp_lt ltSet(&ltData);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim_manyNodes(ic, ctrl0, 0, 2.8, 5, &ltSet, &control);
	
	// ltSet.print();
	ltSet.saveToMat("data/ltSet.mat");
	Arcset_cr3bp_lt temp(&ltData);
	std::vector<ControlLaw*> loadedLaws {};
	temp.readFromMat("data/ltSet.mat", loadedLaws);

	return EXIT_SUCCESS;
}
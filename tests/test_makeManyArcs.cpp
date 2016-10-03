/**
 * This script takes an arc, discretizes it, and then attempts to propagate
 *	lots of arcs using one simulation engine. This should test the ability of the engine
 *	to be used over and over again
 */

#include "Node.hpp"
#include "Nodeset_cr3bp.hpp"
#include "SimEngine.hpp"
#include "SysData_cr3bp.hpp"
#include "Traj_cr3bp.hpp"

using namespace astrohelion;

int main(){

	// Initial conditions at the END of the GMAT-propagated trajectory
	double IC[] = {0.959162150656959, -0.00452492253057802, -0.0125506108074828, 0.0980836599967397, 0.140670049326758, -0.274908800793288};
	double tof = 0.137;	// enough time to complete ~1 rev

	// Generate a set of nodes around the arc
	SysData_cr3bp sys("earth", "moon");
	Nodeset_cr3bp scienceOrbitNodes(&sys, IC, tof, 10);

	SimEngine sim;
	std::vector<Traj_cr3bp> allArcs;
	for(int n = 0; n < scienceOrbitNodes.getNumNodes(); n++){
		printf("Simulating Arc #%03d\n", n);
		std::vector<double> nodeIC = scienceOrbitNodes.getNodeByIx(n).getState();
		double tof = 4;
		Traj_cr3bp traj(&sys);
		sim.runSim(nodeIC, tof, &traj);
		allArcs.push_back(traj);
	}

	return EXIT_SUCCESS;
}
/**
 * This script takes an arc, discretizes it, and then attempts to propagate
 *	lots of arcs using one simulation engine. This should test the ability of the engine
 *	to be used over and over again
 */

#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"

int main(){

	// Initial conditions at the END of the GMAT-propagated trajectory
	double IC[] = {0.959162150656959, -0.00452492253057802, -0.0125506108074828, 0.0980836599967397, 0.140670049326758, -0.274908800793288};
	double tof = 0.137;	// enough time to complete ~1 rev

	// Generate a set of nodes around the arc
	tpat_sys_data_cr3bp sys("earth", "moon");
	tpat_nodeset_cr3bp scienceOrbitNodes(IC, sys, tof, 10);

	tpat_simulation_engine sim(&sys);
	std::vector<tpat_traj_cr3bp> allArcs;
	for(int n = 0; n < scienceOrbitNodes.getNumNodes(); n++){
		printf("Simulating Arc #%03d\n", n);
		std::vector<double> nodeIC = scienceOrbitNodes.getNode(n).getPosVelState();
		double tof = 4;
		sim.runSim(nodeIC, tof);
		allArcs.push_back(sim.getCR3BP_Traj());
	}

	return EXIT_SUCCESS;
}
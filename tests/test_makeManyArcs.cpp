/**
 * This script takes an arc, discretizes it, and then attempts to propagate
 *	lots of arcs using one simulation engine. This should test the ability of the engine
 *	to be used over and over again
 */

#include "tpat_node.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sim_engine.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"

int main(){

	// Initial conditions at the END of the GMAT-propagated trajectory
	double IC[] = {0.959162150656959, -0.00452492253057802, -0.0125506108074828, 0.0980836599967397, 0.140670049326758, -0.274908800793288};
	double tof = 0.137;	// enough time to complete ~1 rev

	// Generate a set of nodes around the arc
	TPAT_Sys_Data_CR3BP sys("earth", "moon");
	TPAT_Nodeset_CR3BP scienceOrbitNodes(IC, &sys, tof, 10);

	TPAT_Sim_Engine sim;
	std::vector<TPAT_Traj_CR3BP> allArcs;
	for(int n = 0; n < scienceOrbitNodes.getNumNodes(); n++){
		printf("Simulating Arc #%03d\n", n);
		std::vector<double> nodeIC = scienceOrbitNodes.getNodeByIx(n).getState();
		double tof = 4;
		TPAT_Traj_CR3BP traj(&sys);
		sim.runSim(nodeIC, tof, &traj);
		allArcs.push_back(traj);
	}

	return EXIT_SUCCESS;
}
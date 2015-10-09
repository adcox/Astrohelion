/**
 *	This script propagates an estimated science orbit, then discretizes it and
 *	runs backwards from points on the arc using the constant low-thrust model.
 */
#include "tpat_event.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_sys_data_cr3bp_ltvp.hpp"
#include "tpat_traj_cr3bp_ltvp.hpp"

int main(){

	double JC_L1 = 3.18834110539543;
	double JC_L2 = 3.17216045039482;
	double JC_L3 = 3.01214714934162;

	// Initial conditions at the END of the GMAT-propagated trajectory
	double IC[] = {0.959162150656959, -0.00452492253057802, -0.0125506108074828, 0.0980836599967397, 0.140670049326758, -0.274908800793288};
	double tof = 0.137;	// enough time to complete ~1 rev

	// Generate a set of nodes around the arc
	tpat_sys_data_cr3bp sys("earth", "moon");
	tpat_nodeset_cr3bp scienceOrbitNodes(IC, sys, tof, 10);

	tpat_sys_data_cr3bp_ltvp sys_lt("earth", "Moon", 0.0012, 2500, 8);
	tpat_simulation_engine sim(&sys_lt);
	sim.setRevTime(true);
	sim.setVerbose(true);
	
	// Add event to stop integration when a certain Jacobi is reached
	tpat_event JC_event(&sys_lt, tpat_event::JC, 0, true, &JC_L2);
	sim.addEvent(JC_event);

	std::vector<tpat_traj_cr3bp_ltvp> allArcs;
	for(int n = 0; n < scienceOrbitNodes.getNumNodes(); n++){
		printf("Simulating Arc #%03d\n", n);
		std::vector<double> nodeIC = scienceOrbitNodes.getNode(n).getPosVelState();
		double tof = 70;
		sim.runSim(nodeIC, tof);
		allArcs.push_back(sim.getCR3BP_LTVP_Traj());
		char filename[256];
		sprintf(filename, "PropFromScienceOrbit_Arcs/%04d.mat", n);
		allArcs.back().saveToMat(filename);
	}

	return EXIT_SUCCESS;
}
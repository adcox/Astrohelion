/**
 *	Generate a Lissajous trajectory "Tom's Method"
 *	
 *	Compile by running:
 *		
 *		>> make generateLisTom
 */

#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_event.hpp"
#include "tpat_family_cr3bp.hpp"
#include "tpat_family_member_cr3bp.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"

#include <cmath>
#include <cstdio>

int main(){

	tpat_sys_data_cr3bp sys("sun", "earth");

	// Load Vertical family
	tpat_family_cr3bp vertFam("../share/families_pac_checked/SE_L1_Vert_small_PAC.mat");

	// Get a vertical orbit (by inspectinv amplitdues from data)
	// tpat_family_member_cr3bp vertMember = vertFam.getMember(111);
	tpat_family_member_cr3bp vertMember = vertFam.getMember(111);

	// Move initial condition to XY plane
	tpat_simulation_engine sim(&sys);
	sim.addEvent(tpat_event::XY_PLANE, 0, true);
	sim.runSim(vertMember.getIC(), vertMember.getTOF()/2);	// tof/2 should be twice as much time as needed

	tpat_traj_cr3bp quarterArc = sim.getCR3BP_Traj();

	// Compute Nodeset from family member
	tpat_nodeset_cr3bp vertNodes(quarterArc.getState(-1), vertFam.getSysDataPtr(), vertMember.getTOF(), 7);

	// Stack the nodes on top of each other
	int numStack = 30;
	tpat_nodeset_cr3bp stackedNodes = vertNodes;
	for(int n = 0; n < numStack; n++){
		tpat_nodeset_cr3bp temp = stackedNodes + vertNodes;
		stackedNodes = temp;
	}

	// Constrain first node to have same x as original vertical, step away in y, fix z to 0
	double data[] = {vertMember.getIC().at(0), vertMember.getIC().at(1) + 0.001, 0, NAN, NAN, NAN};
	tpat_constraint initNodeCon(tpat_constraint::STATE, 0, data, 6);

	stackedNodes.addConstraint(initNodeCon);

	// Correct the entire trajectory
	tpat_correction_engine corrector;
	
	double zAmp = 0, yAmp = 0;
	tpat_nodeset_cr3bp correctedNodes(&sys);

	// Try fixing the 150th node to have a z value of 200000 km
	double maxZData[] = {NAN, NAN, 150000/sys.getCharL(), NAN, NAN, 0};
	tpat_constraint maxZCon(tpat_constraint::STATE, 149, maxZData, 6);
	stackedNodes.addConstraint(maxZCon);

	double maxY = 0;
	do{
		corrector.correct(&stackedNodes);

		correctedNodes = corrector.getCR3BP_Output();

		std::vector<double> allZ = correctedNodes.getCoord(2);
		std::vector<double> allY = correctedNodes.getCoord(1);

		zAmp = *std::max_element(allZ.begin(), allZ.end());
		yAmp = *std::max_element(allY.begin(), allY.end());

		printf("yAmp = %.4f km\n", yAmp*sys.getCharL());
		printf("zAmp = %.4f km\n", zAmp*sys.getCharL());

		std::vector<double> node = correctedNodes.getState(87);
		maxY = node[1]*sys.getCharL() + 2000;
		printf("> maxY = %.2f km, targeting %.2f km\n", node[1]*sys.getCharL(), maxY);

		// Update constraint data
		double maxYData[] = {NAN, maxY/sys.getCharL(), NAN, NAN, 0, NAN};
		tpat_constraint maxYCon(tpat_constraint::STATE, 87, maxYData, 6);
		correctedNodes.clearConstraints();
		correctedNodes.addConstraint(maxZCon);
		correctedNodes.addConstraint(maxYCon);

		// Update the set to correct
		stackedNodes = correctedNodes;
	}while(maxY < 800000);

	correctedNodes.saveToMat("data/Liss_TomMethod_Nodes.mat");
	
	tpat_traj_cr3bp lissTraj = tpat_traj_cr3bp::fromNodeset(correctedNodes);
	lissTraj.saveToMat("data/Liss_TomMethod_Traj.mat");

	// Iterate, moving initial node farther from the original vertical in the y-direction

	return EXIT_SUCCESS;
}
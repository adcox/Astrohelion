/*
 *	Try targeting the saddle point; for debugging
 */

#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_nodeset_bcr4bpr.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"


void testSP_ExactTarget(){
	tpat_sys_data_bcr4bpr semSys("Sun", "earth", "moon");
	double IC[] = {98.529563503083, 1.18065951206182, 0.00114414803935654, 0.270235338496738, -0.222413502966165, 1.71361901677108e-05};
	double t0 = 51.23235;
	double tof = 5;

	tpat_nodeset_bcr4bpr nodes0(IC, &semSys, t0, tof, 2);

	double spConData = 100/semSys.getCharL();
	printf("100 km in SEM system = %.6f non-dim units\n", spConData);
	tpat_constraint spCon(tpat_constraint::SP_RANGE, 1, &spConData, 1);

	nodes0.addConstraint(spCon);

	nodes0.print();

	tpat_correction_engine corrector;
	// corrector.setVerbose(true);
	corrector.setVarTime(false);
	corrector.correct(&nodes0);
	tpat_nodeset_bcr4bpr cNodes = corrector.getBCR4BPR_Output();

	printf("Corrected Nodes:\n");
	cNodes.print();

	// Get SP location at epoch corresponding to final node
	tpat_nodeset_bcr4bpr spNode(cNodes, 1, 1);
	printf("Saddle Point Node:\n");
	spNode.print();
}

int main(){

	testSP_ExactTarget();

	return EXIT_SUCCESS;
}
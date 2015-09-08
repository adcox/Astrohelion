/*
 *	Try targeting the saddle point; for debugging
 */

#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_nodeset_bcr4bpr.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"


int main(){

	tpat_sys_data_bcr4bpr semSys("Sun", "earth", "moon");
	double IC[] = {98.529563503083, 1.18065951206182, 0.00114414803935654, 0.270235338496738, -0.222413502966165, 1.71361901677108e-05};
	double t0 = 51.23235;
	double tof = 5;

	tpat_nodeset_bcr4bpr nodes0(IC, &semSys, t0, tof, 2);
	tpat_constraint spCon(tpat_constraint::SP, 1, IC,1);
	nodes0.addConstraint(spCon);

	nodes0.print();

	tpat_correction_engine corrector;
	// corrector.setVerbose(true);
	corrector.correct(&nodes0);

	return EXIT_SUCCESS;
}
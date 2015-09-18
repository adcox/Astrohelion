#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"

#include <vector>
#include <iostream>

int main(void){
	tpat_sys_data_cr3bp sys("earth", "moon");

	// // ICs for a 2:5 Resonant Orbit in the EM System
	// double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	// double tof = 31.00065761;

	// // Create a nodeset with 10 nodes
	// tpat_nodeset_cr3bp nodeset(IC, &sys, tof, 15);
	// tpat_traj_cr3bp traj = tpat_traj_cr3bp::fromNodeset(nodeset);

	// nodeset.saveToMat("resNodes.mat");
	// traj.saveToMat("traj.mat");

	// // Constraint node 8 to be perpendicular to XZ plane
	// double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	// tpat_constraint perpCross(tpat_constraint::STATE, 7, perpCrossData, 6);

	// // Also constraint final state to be perpendicular
	// tpat_constraint perpCrossEnd(tpat_constraint::STATE, 14, perpCrossData, 6);

	// double almostIC[] =  {IC[0], 0, 0, NAN, NAN, NAN};
	// tpat_constraint icCon(tpat_constraint::STATE, 0, almostIC, 6);

	// nodeset.addConstraint(icCon);
	// nodeset.addConstraint(perpCross);
	// nodeset.addConstraint(perpCrossEnd);

	tpat_correction_engine corrector;
	// corrector.correct(&nodeset);
	// tpat_nodeset_cr3bp correctedNodeset = corrector.getCR3BP_Output();
	// correctedNodeset.saveToMat("resNodes_Corrected.mat");

	//----------------------------------------//
	//    Try a case with equal arc times     //
	//----------------------------------------//

	double IC_halo[] = {0.9900, 0, -0.0203, 0, -1.0674, 0};
	double tof_halo = 1.8632;

	// Create a nodeset
	int halo_nodes = 6;
	tpat_nodeset_cr3bp halo(IC_halo, &sys, tof_halo, halo_nodes);

	double periodic_conData[] = {0,0,0,0,NAN,0};
	tpat_constraint periodicity(tpat_constraint::MATCH_CUST, halo_nodes-1, periodic_conData, 6);

	tpat_constraint fix_halo_x0(tpat_constraint::STATE, 0, IC_halo, 1);

	halo.addConstraint(fix_halo_x0);
	halo.addConstraint(periodicity);

	halo.print();

	corrector.setVarTime(true);
	corrector.setEqualArcTime(true);
	corrector.correct(&halo);
	tpat_nodeset_cr3bp correctedHalo = corrector.getCR3BP_Output();
	correctedHalo.saveToMat("CorrectedHalo.mat");

	return EXIT_SUCCESS;
}
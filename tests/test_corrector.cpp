#include "tpat_ascii_output.hpp"
#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_multShoot_data.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_utilities.hpp"

#include <vector>
#include <iostream>

using namespace std;
static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

void correctEMRes(){
	printColor(BOLDBLACK, "Correct EM 2:5 Resonant Orbit:\n");
	tpat_sys_data_cr3bp sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = 31.00065761;

	// Create a nodeset with
	tpat_nodeset_cr3bp nodeset(IC, &sys, tof, 15);
	tpat_traj_cr3bp traj = tpat_traj_cr3bp::fromNodeset(nodeset);

	nodeset.saveToMat("resNodes.mat");
	traj.saveToMat("traj.mat");

	// Constraint node 07 to be perpendicular to XZ plane
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	tpat_constraint perpCross(tpat_constraint::STATE, 7, perpCrossData, 6);

	// Also constraint final state to be perpendicular
	tpat_constraint perpCrossEnd(tpat_constraint::STATE, 14, perpCrossData, 6);

	double almostIC[] =  {IC[0], 0, 0, NAN, NAN, NAN};
	tpat_constraint icCon(tpat_constraint::STATE, 0, almostIC, 6);

	nodeset.addConstraint(icCon);
	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(perpCrossEnd);

	tpat_correction_engine corrector;
	tpat_nodeset_cr3bp correctedNodeset(&sys);
	
	try{
		// corrector.setVerbose(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(tpat_diverge &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(tpat_exception &e){
		cout << "Successful correction: " << FAIL << endl;
		printErr("%s\n", e.what());
	}

	correctedNodeset.saveToMat("resNodes_Corrected.mat");
}//====================================================

void correctEMRes_EqualArcTime(){
	printColor(BOLDBLACK, "Correct EM 2:5 Resonant Orbit (Equal Arc Time):\n");
	double IC_halo[] = {0.9900, 0, -0.0203, 0, -1.0674, 0};
	double tof_halo = 1.8632;
	tpat_sys_data_cr3bp sys("earth", "moon");

	// Create a nodeset
	int halo_nodes = 6;
	tpat_nodeset_cr3bp halo(IC_halo, &sys, tof_halo, halo_nodes);

	double periodic_conData[] = {0,0,0,0,NAN,0};
	tpat_constraint periodicity(tpat_constraint::MATCH_CUST, halo_nodes-1, periodic_conData, 6);

	tpat_constraint fix_halo_x0(tpat_constraint::STATE, 0, IC_halo, 1);

	halo.addConstraint(fix_halo_x0);
	halo.addConstraint(periodicity);

	tpat_correction_engine corrector;
	tpat_nodeset_cr3bp correctedHalo(&sys);

	try{
		// corrector.setVerbose(true);
		corrector.setVarTime(true);
		corrector.setEqualArcTime(true);
		corrector.multShoot(&halo, &correctedHalo);
		cout << "Successful correction: " << PASS << endl;
	}catch(tpat_diverge &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(tpat_exception &e){
		cout << "Successful correction: " << FAIL << endl;
		printErr("%s\n", e.what());
	}

	correctedHalo.saveToMat("CorrectedHalo.mat");
}//====================================================

void correctEMRes_revTime(){
	printColor(BOLDBLACK, "Correct EM 2:5 Resonant Orbit (Reverse Time):\n");
	tpat_sys_data_cr3bp sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = -31.00065761;

	// Create a nodeset
	tpat_nodeset_cr3bp nodeset(IC, &sys, tof, 15);

	// Constraint node 07 to be perpendicular to XZ plane
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	tpat_constraint perpCross(tpat_constraint::STATE, 7, perpCrossData, 6);

	// Also constraint final state to be perpendicular
	tpat_constraint perpCrossEnd(tpat_constraint::STATE, 14, perpCrossData, 6);

	double almostIC[] =  {IC[0], 0, 0, NAN, NAN, NAN};
	tpat_constraint icCon(tpat_constraint::STATE, 0, almostIC, 6);

	nodeset.addConstraint(icCon);
	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(perpCrossEnd);

	nodeset.putInChronoOrder();
	// nodeset.print();
	// nodeset.printInChrono();

	tpat_correction_engine corrector;
	tpat_nodeset_cr3bp correctedNodeset(&sys);
	
	try{
		// corrector.setVerbose(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(tpat_diverge &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(tpat_exception &e){
		cout << "Successful correction: " << FAIL << endl;
		printErr("%s\n", e.what());
	}
}//====================================================

void correctEMRes_doubleSource(){
	printColor(BOLDBLACK, "Correct EM 2:5 Resonant Orbit (Reverse Time):\n");
	tpat_sys_data_cr3bp sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = 31.00065761;

	tpat_nodeset_cr3bp posTimeArc(IC, &sys, tof/2, 8);
	tpat_nodeset_cr3bp revTimeArc(IC, &sys, -tof/2, 8);

	tpat_nodeset_cr3bp nodeset = posTimeArc;
	nodeset.appendSetAtNode(&revTimeArc, 0, 0, 0);

	// nodeset.print();
	// nodeset.printInChrono();

	// Constraint node 07 to be perpendicular to XZ plane
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	tpat_constraint perpCross(tpat_constraint::STATE, 7, perpCrossData, 6);

	// Also constraint final state to be perpendicular
	tpat_constraint perpCrossEnd(tpat_constraint::STATE, 14, perpCrossData, 6);

	double almostIC[] =  {IC[0], 0, 0, NAN, NAN, NAN};
	tpat_constraint icCon(tpat_constraint::STATE, 0, almostIC, 6);

	nodeset.addConstraint(icCon);
	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(perpCrossEnd);

	nodeset.putInChronoOrder();
	// nodeset.print();
	// nodeset.printInChrono();

	tpat_correction_engine corrector;
	tpat_nodeset_cr3bp correctedNodeset(&sys);
	
	try{
		// corrector.setVerbose(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(tpat_diverge &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(tpat_exception &e){
		cout << "Successful correction: " << FAIL << endl;
		printErr("%s\n", e.what());
	}
}//====================================================

void correctEMRes_doubleSource_irregular(){
	printColor(BOLDBLACK, "Correct EM 2:5 Resonant Orbit (Reverse Time, irrgularly spaced nodes):\n");
	tpat_sys_data_cr3bp sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = 31.00065761;

	tpat_nodeset_cr3bp posTimeArc(IC, &sys, tof/2, 8);
	tpat_nodeset_cr3bp revTimeArc(IC, &sys, -tof/2, 8);

	tpat_nodeset_cr3bp nodeset = posTimeArc;
	nodeset.appendSetAtNode(&revTimeArc, 0, 0, 0);

	// nodeset.print();
	// nodeset.printInChrono();
	// nodeset.printNodeIDMap();
	// nodeset.printSegIDMap();

	// Constraint node 07 to be perpendicular to XZ plane
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	tpat_constraint perpCross(tpat_constraint::STATE, 7, perpCrossData, 6);

	// Also constraint final state to be perpendicular
	tpat_constraint perpCrossEnd(tpat_constraint::STATE, 14, perpCrossData, 6);

	double almostIC[] =  {IC[0], 0, 0, NAN, NAN, NAN};
	tpat_constraint icCon(tpat_constraint::STATE, 0, almostIC, 6);

	nodeset.addConstraint(icCon);
	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(perpCrossEnd);

	nodeset.deleteNode(2);
	nodeset.deleteNode(5);
	nodeset.deleteNode(11);
	nodeset.deleteNode(12);

	nodeset.putInChronoOrder();
	// nodeset.print();
	// nodeset.printInChrono();
	// nodeset.printNodeIDMap();
	// nodeset.printSegIDMap();

	tpat_correction_engine corrector;
	tpat_nodeset_cr3bp correctedNodeset(&sys);
	
	try{
		// corrector.setVerbose(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(tpat_diverge &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(tpat_exception &e){
		cout << "Successful correction: " << FAIL << endl;
		printErr("%s\n", e.what());
	}
}//====================================================

int main(void){
	correctEMRes();
	correctEMRes_EqualArcTime();	
	correctEMRes_revTime();
	correctEMRes_doubleSource();
	correctEMRes_doubleSource_irregular();
	
	return EXIT_SUCCESS;
}
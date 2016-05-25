#include "tpat_ascii_output.hpp"
#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_multShoot_data.hpp"
#include "tpat_nodeset_bc4bp.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_sys_data_bc4bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_utilities.hpp"

#include <vector>
#include <iostream>

using namespace std;
static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

void correctEMRes(){
	printColor(BOLDBLACK, "CR3BP EM 2:5 Resonant Orbit:\n");
	TPAT_Sys_Data_CR3BP sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = 31.00065761;

	// Create a nodeset with
	TPAT_Nodeset_CR3BP nodeset(IC, &sys, tof, 15);
	TPAT_Traj_CR3BP traj = TPAT_Traj_CR3BP::fromNodeset(nodeset);

	nodeset.saveToMat("resNodes.mat");
	traj.saveToMat("traj.mat");

	// Constraint node 07 to be perpendicular to XZ plane
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	TPAT_Constraint perpCross(TPAT_Constraint_Tp::STATE, 7, perpCrossData, 6);

	// Also constraint final state to be perpendicular
	TPAT_Constraint perpCrossEnd(TPAT_Constraint_Tp::STATE, 14, perpCrossData, 6);

	double almostIC[] =  {IC[0], 0, 0, NAN, NAN, NAN};
	TPAT_Constraint icCon(TPAT_Constraint_Tp::STATE, 0, almostIC, 6);

	nodeset.addConstraint(icCon);
	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(perpCrossEnd);

	TPAT_Correction_Engine corrector;
	TPAT_Nodeset_CR3BP correctedNodeset(&sys);
	
	try{
		// corrector.setVerbose(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(TPAT_Diverge &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(TPAT_Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		printErr("%s\n", e.what());
	}

	correctedNodeset.saveToMat("resNodes_Corrected.mat");
}//====================================================

void correctEMRes_EqualArcTime(){
	printColor(BOLDBLACK, "CR3BP EM 2:5 Resonant Orbit (Equal Arc Time):\n");
	double IC_halo[] = {0.9900, 0, -0.0203, 0, -1.0674, 0};
	double tof_halo = 1.8632;
	TPAT_Sys_Data_CR3BP sys("earth", "moon");

	// Create a nodeset
	int halo_nodes = 6;
	TPAT_Nodeset_CR3BP halo(IC_halo, &sys, tof_halo, halo_nodes);

	double periodic_conData[] = {0,0,0,0,NAN,0};
	TPAT_Constraint periodicity(TPAT_Constraint_Tp::MATCH_CUST, halo_nodes-1, periodic_conData, 6);

	TPAT_Constraint fix_halo_x0(TPAT_Constraint_Tp::STATE, 0, IC_halo, 1);

	halo.addConstraint(fix_halo_x0);
	halo.addConstraint(periodicity);

	TPAT_Correction_Engine corrector;
	TPAT_Nodeset_CR3BP correctedHalo(&sys);

	try{
		// corrector.setVerbose(true);
		corrector.setVarTime(true);
		corrector.setEqualArcTime(true);
		corrector.multShoot(&halo, &correctedHalo);
		cout << "Successful correction: " << PASS << endl;
	}catch(TPAT_Diverge &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(TPAT_Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		printErr("%s\n", e.what());
	}

	correctedHalo.saveToMat("CorrectedHalo.mat");
}//====================================================

void correctEMRes_revTime(){
	printColor(BOLDBLACK, "CR3BP EM 2:5 Resonant Orbit (Reverse Time):\n");
	TPAT_Sys_Data_CR3BP sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = -31.00065761;

	// Create a nodeset
	TPAT_Nodeset_CR3BP nodeset(IC, &sys, tof, 15);

	// Constraint node 07 to be perpendicular to XZ plane
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	TPAT_Constraint perpCross(TPAT_Constraint_Tp::STATE, 7, perpCrossData, 6);

	// Also constraint final state to be perpendicular
	TPAT_Constraint perpCrossEnd(TPAT_Constraint_Tp::STATE, 14, perpCrossData, 6);

	double almostIC[] =  {IC[0], 0, 0, NAN, NAN, NAN};
	TPAT_Constraint icCon(TPAT_Constraint_Tp::STATE, 0, almostIC, 6);

	nodeset.addConstraint(icCon);
	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(perpCrossEnd);

	nodeset.putInChronoOrder();
	// nodeset.print();
	// nodeset.printInChrono();

	TPAT_Correction_Engine corrector;
	TPAT_Nodeset_CR3BP correctedNodeset(&sys);
	
	try{
		// corrector.setVerbose(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(TPAT_Diverge &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(TPAT_Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		printErr("%s\n", e.what());
	}
}//====================================================

void correctEMRes_doubleSource(){
	printColor(BOLDBLACK, "CR3BP EM 2:5 Resonant Orbit (Reverse Time):\n");
	TPAT_Sys_Data_CR3BP sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = 31.00065761;

	TPAT_Nodeset_CR3BP posTimeArc(IC, &sys, tof/2, 8);
	TPAT_Nodeset_CR3BP revTimeArc(IC, &sys, -tof/2, 8);

	TPAT_Nodeset_CR3BP nodeset = posTimeArc;
	nodeset.appendSetAtNode(&revTimeArc, 0, 0, 0);

	// nodeset.print();
	// nodeset.printInChrono();

	// Constraint node 07 to be perpendicular to XZ plane
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	TPAT_Constraint perpCross(TPAT_Constraint_Tp::STATE, 7, perpCrossData, 6);

	// Also constraint final state to be perpendicular
	TPAT_Constraint perpCrossEnd(TPAT_Constraint_Tp::STATE, 14, perpCrossData, 6);

	double almostIC[] =  {IC[0], 0, 0, NAN, NAN, NAN};
	TPAT_Constraint icCon(TPAT_Constraint_Tp::STATE, 0, almostIC, 6);

	nodeset.addConstraint(icCon);
	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(perpCrossEnd);

	nodeset.putInChronoOrder();
	// nodeset.print();
	// nodeset.printInChrono();

	TPAT_Correction_Engine corrector;
	TPAT_Nodeset_CR3BP correctedNodeset(&sys);
	
	try{
		// corrector.setVerbose(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(TPAT_Diverge &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(TPAT_Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		printErr("%s\n", e.what());
	}
}//====================================================

void correctEMRes_doubleSource_irregular(){
	printColor(BOLDBLACK, "CR3BP EM 2:5 Resonant Orbit (Reverse Time, irrgularly spaced nodes):\n");
	TPAT_Sys_Data_CR3BP sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = 31.00065761;

	TPAT_Nodeset_CR3BP posTimeArc(IC, &sys, tof/2, 8);
	TPAT_Nodeset_CR3BP revTimeArc(IC, &sys, -tof/2, 8);

	TPAT_Nodeset_CR3BP nodeset = posTimeArc;
	nodeset.appendSetAtNode(&revTimeArc, 0, 0, 0);

	// nodeset.print();
	// nodeset.printInChrono();
	// nodeset.printNodeIDMap();
	// nodeset.printSegIDMap();

	// Constraint node 07 to be perpendicular to XZ plane
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	TPAT_Constraint perpCross(TPAT_Constraint_Tp::STATE, 7, perpCrossData, 6);

	// Also constraint final state to be perpendicular
	TPAT_Constraint perpCrossEnd(TPAT_Constraint_Tp::STATE, 14, perpCrossData, 6);

	double almostIC[] =  {IC[0], 0, 0, NAN, NAN, NAN};
	TPAT_Constraint icCon(TPAT_Constraint_Tp::STATE, 0, almostIC, 6);

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

	TPAT_Correction_Engine corrector;
	TPAT_Nodeset_CR3BP correctedNodeset(&sys);
	
	try{
		// corrector.setVerbose(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(TPAT_Diverge &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(TPAT_Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		printErr("%s\n", e.what());
	}
}//====================================================

void correctSEMHalo(){
	printColor(BOLDBLACK, "BC4BP SEM L1 Quasi-Halo:\n");
	TPAT_Sys_Data_BC4BP sys("Sun", "earth", "moon");
	std::vector<double> haloIC {-1.144739, 0, 0.089011, 0, 0.011608, 0};
	double tof = 310;

	TPAT_Nodeset_BC4BP nodeset(haloIC, &sys, 0, tof, 7);
	nodeset.saveToMat("bc4bp_halo_raw.mat");
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	TPAT_Constraint perpCross(TPAT_Constraint_Tp::STATE, 0, perpCrossData, 6);

	double xzPlaneData[] = {NAN,0,NAN,NAN,NAN,NAN};
	TPAT_Constraint xzPlaneCon0(TPAT_Constraint_Tp::STATE, 3, xzPlaneData, 6);
	TPAT_Constraint xzPlaneConF(TPAT_Constraint_Tp::STATE, 6, xzPlaneData, 6);

	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(xzPlaneCon0);
	nodeset.addConstraint(xzPlaneConF);

	TPAT_Correction_Engine corrector;
	TPAT_Nodeset_BC4BP correctedNodeset(&sys);
	
	try{
		// corrector.setVerbose(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(TPAT_Diverge &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(TPAT_Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		printErr("%s\n", e.what());
	}

	correctedNodeset.saveToMat("bc4bp_halo.mat");
}//====================================================

void correctSEMHalo_revTime(){
	printColor(BOLDBLACK, "BC4BP SEM L1 Quasi-Halo (Reverse Time):\n");
	TPAT_Sys_Data_BC4BP sys("Sun", "earth", "moon");
	std::vector<double> haloIC {-1.144739, 0, 0.089011, 0, 0.011608, 0};
	double tof = -310;

	TPAT_Nodeset_BC4BP nodeset(haloIC, &sys, 0, tof, 7);
	nodeset.saveToMat("bc4bp_halo_raw.mat");

	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	TPAT_Constraint perpCross(TPAT_Constraint_Tp::STATE, 0, perpCrossData, 6);

	double xzPlaneData[] = {NAN,0,NAN,NAN,NAN,NAN};
	TPAT_Constraint xzPlaneCon0(TPAT_Constraint_Tp::STATE, 3, xzPlaneData, 6);
	TPAT_Constraint xzPlaneConF(TPAT_Constraint_Tp::STATE, 6, xzPlaneData, 6);

	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(xzPlaneCon0);
	nodeset.addConstraint(xzPlaneConF);

	TPAT_Correction_Engine corrector;
	TPAT_Nodeset_BC4BP correctedNodeset(&sys);
	
	try{
		// corrector.setVerbose(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(TPAT_Diverge &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(TPAT_Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		printErr("%s\n", e.what());
	}
	correctedNodeset.saveToMat("bc4bp_halo.mat");
}//====================================================

void correctSEMHalo_doubleSource(){
	printColor(BOLDBLACK, "BC4BP SEM L1 Quasi-Halo (Double Source):\n");
	TPAT_Sys_Data_BC4BP sys("Sun", "earth", "moon");
	std::vector<double> haloIC {-1.144739, 0, 0.089011, 0, 0.011608, 0};
	double tof = 310;

	TPAT_Nodeset_BC4BP posTimeArc(haloIC, &sys, 0, tof/2, 4);
	TPAT_Nodeset_BC4BP revTimeArc(haloIC, &sys, 0, -tof/2, 4);

	TPAT_Nodeset_BC4BP nodeset = posTimeArc;
	nodeset.appendSetAtNode(&revTimeArc, 0, 0, 0);

	// nodeset.print();
	// nodeset.printInChrono();

	// nodeset.saveToMat("bc4bp_halo_raw.mat");
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	TPAT_Constraint perpCross(TPAT_Constraint_Tp::STATE, 0, perpCrossData, 6);

	double xzPlaneData[] = {NAN,0,NAN,NAN,NAN,NAN};
	TPAT_Constraint xzPlaneCon0(TPAT_Constraint_Tp::STATE, 3, xzPlaneData, 6);
	TPAT_Constraint xzPlaneConF(TPAT_Constraint_Tp::STATE, 6, xzPlaneData, 6);

	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(xzPlaneCon0);
	nodeset.addConstraint(xzPlaneConF);

	TPAT_Correction_Engine corrector;
	TPAT_Nodeset_BC4BP correctedNodeset(&sys);
	
	try{
		// corrector.setVerbose(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(TPAT_Diverge &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(TPAT_Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		printErr("%s\n", e.what());
	}

	// correctedNodeset.saveToMat("bc4bp_halo.mat");
}//====================================================

void correctSEMHalo_doubleSource_irregular(){
	printColor(BOLDBLACK, "BC4BP SEM L1 Quasi-Halo (Double Source, Irregular Nodes):\n");
	TPAT_Sys_Data_BC4BP sys("Sun", "earth", "moon");
	std::vector<double> haloIC {-1.144739, 0, 0.089011, 0, 0.011608, 0};
	double tof = 310;

	TPAT_Nodeset_BC4BP posTimeArc(haloIC, &sys, 0, tof/2, 4);
	TPAT_Nodeset_BC4BP revTimeArc(haloIC, &sys, 0, -tof/2, 4);

	TPAT_Nodeset_BC4BP nodeset = posTimeArc;
	nodeset.appendSetAtNode(&revTimeArc, 0, 0, 0);

	nodeset.deleteNode(1);
	nodeset.deleteNode(5);
	// nodeset.print();
	// nodeset.printInChrono();

	// nodeset.saveToMat("bc4bp_halo_raw.mat");
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	TPAT_Constraint perpCross(TPAT_Constraint_Tp::STATE, 0, perpCrossData, 6);

	double xzPlaneData[] = {NAN,0,NAN,NAN,NAN,NAN};
	TPAT_Constraint xzPlaneCon0(TPAT_Constraint_Tp::STATE, 3, xzPlaneData, 6);
	TPAT_Constraint xzPlaneConF(TPAT_Constraint_Tp::STATE, 6, xzPlaneData, 6);

	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(xzPlaneCon0);
	nodeset.addConstraint(xzPlaneConF);

	TPAT_Correction_Engine corrector;
	TPAT_Nodeset_BC4BP correctedNodeset(&sys);
	
	try{
		// corrector.setVerbose(true);
		corrector.setTol(5e-12);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(TPAT_Diverge &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(TPAT_Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		printErr("%s\n", e.what());
	}

	// correctedNodeset.saveToMat("bc4bp_halo.mat");
}//====================================================

int main(void){
	correctEMRes();
	correctEMRes_EqualArcTime();	
	correctEMRes_revTime();
	correctEMRes_doubleSource();
	correctEMRes_doubleSource_irregular();
	
	correctSEMHalo();
	correctSEMHalo_revTime();
	correctSEMHalo_doubleSource();
	correctSEMHalo_doubleSource_irregular();
	return EXIT_SUCCESS;
}
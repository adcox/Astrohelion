#include "AsciiOutput.hpp"
#include "Constraint.hpp"
#include "CorrectionEngine.hpp"
#include "MultShootData.hpp"
#include "Nodeset_bc4bp.hpp"
#include "Nodeset_cr3bp.hpp"
#include "SysData_bc4bp.hpp"
#include "SysData_cr3bp.hpp"
#include "Traj_cr3bp.hpp"
#include "Utilities.hpp"

#include <vector>
#include <iostream>

using namespace std;
using namespace astrohelion;

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

void correctEMRes(){
	astrohelion::printColor(BOLDBLACK, "CR3BP EM 2:5 Resonant Orbit:\n");
	SysData_cr3bp sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = 31.00065761;

	// Create a nodeset with
	Nodeset_cr3bp nodeset(&sys, IC, tof, 15);
	Traj_cr3bp traj = Traj_cr3bp::fromNodeset(nodeset);

	nodeset.saveToMat("resNodes.mat");
	traj.saveToMat("traj.mat");

	// Constraint node 07 to be perpendicular to XZ plane
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	Constraint perpCross(Constraint_tp::STATE, 7, perpCrossData, 6);

	// Also constraint final state to be perpendicular
	Constraint perpCrossEnd(Constraint_tp::STATE, 14, perpCrossData, 6);

	double almostIC[] =  {IC[0], 0, 0, NAN, NAN, NAN};
	Constraint icCon(Constraint_tp::STATE, 0, almostIC, 6);

	nodeset.addConstraint(icCon);
	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(perpCrossEnd);

	CorrectionEngine corrector;
	Nodeset_cr3bp correctedNodeset(&sys);
	
	try{
		// corrector.setVerbosity(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(DivergeException &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		astrohelion::printErr("%s\n", e.what());
	}

	correctedNodeset.saveToMat("resNodes_Corrected.mat");
}//====================================================

void correctEMRes_EqualArcTime(){
	astrohelion::printColor(BOLDBLACK, "CR3BP EM 2:5 Resonant Orbit (Equal Arc Time):\n");
	double IC_halo[] = {0.9900, 0, -0.0203, 0, -1.0674, 0};
	double tof_halo = 1.8632;
	SysData_cr3bp sys("earth", "moon");

	// Create a nodeset
	int halo_nodes = 6;
	Nodeset_cr3bp halo(&sys, IC_halo, tof_halo, halo_nodes);

	double periodic_conData[] = {0,0,0,0,NAN,0};
	Constraint periodicity(Constraint_tp::MATCH_CUST, halo_nodes-1, periodic_conData, 6);

	Constraint fix_halo_x0(Constraint_tp::STATE, 0, IC_halo, 1);

	halo.addConstraint(fix_halo_x0);
	halo.addConstraint(periodicity);

	CorrectionEngine corrector;
	Nodeset_cr3bp correctedHalo(&sys);

	try{
		// corrector.setVerbosity(true);
		corrector.setVarTime(true);
		corrector.setEqualArcTime(true);
		corrector.multShoot(&halo, &correctedHalo);
		cout << "Successful correction: " << PASS << endl;
	}catch(DivergeException &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		astrohelion::printErr("%s\n", e.what());
	}

	correctedHalo.saveToMat("CorrectedHalo.mat");
}//====================================================

void correctEMRes_revTime(){
	astrohelion::printColor(BOLDBLACK, "CR3BP EM 2:5 Resonant Orbit (Reverse Time):\n");
	SysData_cr3bp sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = -31.00065761;

	// Create a nodeset
	Nodeset_cr3bp nodeset(&sys, IC, tof, 15);

	// Constraint node 07 to be perpendicular to XZ plane
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	Constraint perpCross(Constraint_tp::STATE, 7, perpCrossData, 6);

	// Also constraint final state to be perpendicular
	Constraint perpCrossEnd(Constraint_tp::STATE, 14, perpCrossData, 6);

	double almostIC[] =  {IC[0], 0, 0, NAN, NAN, NAN};
	Constraint icCon(Constraint_tp::STATE, 0, almostIC, 6);

	nodeset.addConstraint(icCon);
	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(perpCrossEnd);

	nodeset.putInChronoOrder();
	// nodeset.print();
	// nodeset.printInChrono();

	CorrectionEngine corrector;
	Nodeset_cr3bp correctedNodeset(&sys);
	
	try{
		// corrector.setVerbosity(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(DivergeException &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		astrohelion::printErr("%s\n", e.what());
	}
}//====================================================

void correctEMRes_doubleSource(){
	astrohelion::printColor(BOLDBLACK, "CR3BP EM 2:5 Resonant Orbit (Reverse Time):\n");
	SysData_cr3bp sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = 31.00065761;

	Nodeset_cr3bp posTimeArc(&sys, IC, tof/2, 8);
	Nodeset_cr3bp revTimeArc(&sys, IC, -tof/2, 8);

	Nodeset_cr3bp nodeset = posTimeArc;
	nodeset.appendSetAtNode(&revTimeArc, 0, 0, 0);

	// nodeset.print();
	// nodeset.printInChrono();

	// Constraint node 07 to be perpendicular to XZ plane
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	Constraint perpCross(Constraint_tp::STATE, 7, perpCrossData, 6);

	// Also constraint final state to be perpendicular
	Constraint perpCrossEnd(Constraint_tp::STATE, 14, perpCrossData, 6);

	double almostIC[] =  {IC[0], 0, 0, NAN, NAN, NAN};
	Constraint icCon(Constraint_tp::STATE, 0, almostIC, 6);

	nodeset.addConstraint(icCon);
	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(perpCrossEnd);

	nodeset.putInChronoOrder();
	// nodeset.print();
	// nodeset.printInChrono();

	CorrectionEngine corrector;
	Nodeset_cr3bp correctedNodeset(&sys);
	
	try{
		// corrector.setVerbosity(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(DivergeException &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		astrohelion::printErr("%s\n", e.what());
	}
}//====================================================

void correctEMRes_doubleSource_irregular(){
	astrohelion::printColor(BOLDBLACK, "CR3BP EM 2:5 Resonant Orbit (Reverse Time, irrgularly spaced nodes):\n");
	SysData_cr3bp sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = 31.00065761;

	Nodeset_cr3bp posTimeArc(&sys, IC, tof/2, 8);
	Nodeset_cr3bp revTimeArc(&sys, IC, -tof/2, 8);

	Nodeset_cr3bp nodeset = posTimeArc;
	nodeset.appendSetAtNode(&revTimeArc, 0, 0, 0);

	// nodeset.print();
	// nodeset.printInChrono();
	// nodeset.printNodeIDMap();
	// nodeset.printSegIDMap();

	// Constraint node 07 to be perpendicular to XZ plane
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	Constraint perpCross(Constraint_tp::STATE, 7, perpCrossData, 6);

	// Also constraint final state to be perpendicular
	Constraint perpCrossEnd(Constraint_tp::STATE, 14, perpCrossData, 6);

	double almostIC[] =  {IC[0], 0, 0, NAN, NAN, NAN};
	Constraint icCon(Constraint_tp::STATE, 0, almostIC, 6);

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

	CorrectionEngine corrector;
	Nodeset_cr3bp correctedNodeset(&sys);
	
	try{
		// corrector.setVerbosity(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(DivergeException &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		astrohelion::printErr("%s\n", e.what());
	}
}//====================================================

void correctSEMHalo(){
	astrohelion::printColor(BOLDBLACK, "BC4BP SEM L1 Quasi-Halo:\n");
	SysData_bc4bp sys("Sun", "earth", "moon");
	std::vector<double> haloIC {-1.144739, 0, 0.089011, 0, 0.011608, 0};
	double tof = 310;

	Nodeset_bc4bp nodeset(&sys, haloIC, 0, tof, 7);
	nodeset.saveToMat("bc4bp_halo_raw.mat");
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	Constraint perpCross(Constraint_tp::STATE, 0, perpCrossData, 6);

	double xzPlaneData[] = {NAN,0,NAN,NAN,NAN,NAN};
	Constraint xzPlaneCon0(Constraint_tp::STATE, 3, xzPlaneData, 6);
	Constraint xzPlaneConF(Constraint_tp::STATE, 6, xzPlaneData, 6);

	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(xzPlaneCon0);
	nodeset.addConstraint(xzPlaneConF);

	CorrectionEngine corrector;
	Nodeset_bc4bp correctedNodeset(&sys);
	
	try{
		// corrector.setVerbosity(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(DivergeException &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		astrohelion::printErr("%s\n", e.what());
	}

	correctedNodeset.saveToMat("bc4bp_halo.mat");
}//====================================================

void correctSEMHalo_revTime(){
	astrohelion::printColor(BOLDBLACK, "BC4BP SEM L1 Quasi-Halo (Reverse Time):\n");
	SysData_bc4bp sys("Sun", "earth", "moon");
	std::vector<double> haloIC {-1.144739, 0, 0.089011, 0, 0.011608, 0};
	double tof = -310;

	Nodeset_bc4bp nodeset(&sys, haloIC, 0, tof, 7);
	nodeset.saveToMat("bc4bp_halo_raw.mat");

	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	Constraint perpCross(Constraint_tp::STATE, 0, perpCrossData, 6);

	double xzPlaneData[] = {NAN,0,NAN,NAN,NAN,NAN};
	Constraint xzPlaneCon0(Constraint_tp::STATE, 3, xzPlaneData, 6);
	Constraint xzPlaneConF(Constraint_tp::STATE, 6, xzPlaneData, 6);

	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(xzPlaneCon0);
	nodeset.addConstraint(xzPlaneConF);

	CorrectionEngine corrector;
	Nodeset_bc4bp correctedNodeset(&sys);
	
	try{
		// corrector.setVerbosity(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(DivergeException &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		astrohelion::printErr("%s\n", e.what());
	}
	correctedNodeset.saveToMat("bc4bp_halo.mat");
}//====================================================

void correctSEMHalo_doubleSource(){
	astrohelion::printColor(BOLDBLACK, "BC4BP SEM L1 Quasi-Halo (Double Source):\n");
	SysData_bc4bp sys("Sun", "earth", "moon");
	std::vector<double> haloIC {-1.144739, 0, 0.089011, 0, 0.011608, 0};
	double tof = 310;

	Nodeset_bc4bp posTimeArc(&sys, haloIC, 0, tof/2, 4);
	Nodeset_bc4bp revTimeArc(&sys, haloIC, 0, -tof/2, 4);

	Nodeset_bc4bp nodeset = posTimeArc;
	nodeset.appendSetAtNode(&revTimeArc, 0, 0, 0);

	// nodeset.print();
	// nodeset.printInChrono();

	// nodeset.saveToMat("bc4bp_halo_raw.mat");
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	Constraint perpCross(Constraint_tp::STATE, 0, perpCrossData, 6);

	double xzPlaneData[] = {NAN,0,NAN,NAN,NAN,NAN};
	Constraint xzPlaneCon0(Constraint_tp::STATE, 3, xzPlaneData, 6);
	Constraint xzPlaneConF(Constraint_tp::STATE, 6, xzPlaneData, 6);

	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(xzPlaneCon0);
	nodeset.addConstraint(xzPlaneConF);

	CorrectionEngine corrector;
	Nodeset_bc4bp correctedNodeset(&sys);
	
	try{
		// corrector.setVerbosity(true);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(DivergeException &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		astrohelion::printErr("%s\n", e.what());
	}

	// correctedNodeset.saveToMat("bc4bp_halo.mat");
}//====================================================

void correctSEMHalo_doubleSource_irregular(){
	astrohelion::printColor(BOLDBLACK, "BC4BP SEM L1 Quasi-Halo (Double Source, Irregular Nodes):\n");
	SysData_bc4bp sys("Sun", "earth", "moon");
	std::vector<double> haloIC {-1.144739, 0, 0.089011, 0, 0.011608, 0};
	double tof = 310;

	Nodeset_bc4bp posTimeArc(&sys, haloIC, 0, tof/2, 4);
	Nodeset_bc4bp revTimeArc(&sys, haloIC, 0, -tof/2, 4);

	Nodeset_bc4bp nodeset = posTimeArc;
	nodeset.appendSetAtNode(&revTimeArc, 0, 0, 0);

	nodeset.deleteNode(1);
	nodeset.deleteNode(5);
	// nodeset.print();
	// nodeset.printInChrono();

	// nodeset.saveToMat("bc4bp_halo_raw.mat");
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	Constraint perpCross(Constraint_tp::STATE, 0, perpCrossData, 6);

	double xzPlaneData[] = {NAN,0,NAN,NAN,NAN,NAN};
	Constraint xzPlaneCon0(Constraint_tp::STATE, 3, xzPlaneData, 6);
	Constraint xzPlaneConF(Constraint_tp::STATE, 6, xzPlaneData, 6);

	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(xzPlaneCon0);
	nodeset.addConstraint(xzPlaneConF);

	CorrectionEngine corrector;
	Nodeset_bc4bp correctedNodeset(&sys);
	
	try{
		// corrector.setVerbosity(true);
		corrector.setTol(5e-12);
		corrector.multShoot(&nodeset, &correctedNodeset);
		cout << "Successful correction: " << PASS << endl;
	}catch(DivergeException &e){
		cout << "Successful correction: " << FAIL << endl;
	}catch(Exception &e){
		cout << "Successful correction: " << FAIL << endl;
		astrohelion::printErr("%s\n", e.what());
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
#define BOOST_TEST_MODULE MultShootEngine

#include <boost/test/unit_test.hpp>

#include <vector>
#include <iostream>

#include "Arcset_bc4bp.hpp"
#include "Arcset_cr3bp.hpp"
#include "Constraint.hpp"
#include "MultShootEngine.hpp"
#include "MultShootData.hpp"
#include "SimEngine.hpp"
#include "SysData_bc4bp.hpp"
#include "SysData_cr3bp.hpp"
#include "Utilities.hpp"

using namespace std;
using namespace astrohelion;

BOOST_AUTO_TEST_CASE(CR3BP_EM_IntermediateSolution){
	SysData_cr3bp sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = 31.00065761;

	SimEngine sim;
	Arcset_cr3bp nodeset(&sys), correctedNodeset(&sys);
	sim.runSim_manyNodes(IC, tof, 15, &nodeset);

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

	// Set up multiple shooting engine to only do a few iterations so that error is still large
	// Ignore the large error and return the "corrected" arcset
	MultShootEngine corrector;
	corrector.setMaxIts(2);
	corrector.setIgnoreDiverge(true);
	BOOST_CHECK_NO_THROW(corrector.multShoot(&nodeset, &correctedNodeset));

	// Check to make sure that segments and nodes are consistent - previous bugs had discontinuities
	// between segments and nodes because the free variable vector was updated AFTER the error was calculated
	for(unsigned int s = 0; s < correctedNodeset.getNumSegs(); s++){
		Segment &refSeg = correctedNodeset.getSegRefByIx(s);
		std::vector<double> originState = correctedNodeset.getNodeRef(refSeg.getOrigin()).getState();
		std::vector<double> segState = refSeg.getStateByRow(0);

		for(unsigned int i = 0; i < 6; i++){
			BOOST_CHECK_CLOSE(originState[i], segState[i], 1e-12);
		}

		BOOST_CHECK_CLOSE(correctedNodeset.getNodeRef(refSeg.getOrigin()).getEpoch(), refSeg.getTimeByIx(0), 1e-12);
	}
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_Resonant){
	SysData_cr3bp sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = 31.00065761;

	SimEngine sim;
	Arcset_cr3bp nodeset(&sys), correctedNodeset(&sys);
	sim.runSim_manyNodes(IC, tof, 15, &nodeset);

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

	MultShootEngine corrector;
	
	BOOST_CHECK_NO_THROW(corrector.multShoot(&nodeset, &correctedNodeset));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_Halo_EqualArcTime){
	double IC_halo[] = {0.9900, 0, -0.0203, 0, -1.0674, 0};
	double tof_halo = 1.8632;
	SysData_cr3bp sys("earth", "moon");

	// Create a nodeset
	int halo_nodes = 6;
	Arcset_cr3bp halo(&sys), correctedHalo(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(IC_halo, tof_halo, halo_nodes, &halo);

	double periodic_conData[] = {0,0,0,0,NAN,0};
	Constraint periodicity(Constraint_tp::MATCH_CUST, halo_nodes-1, periodic_conData, 6);

	Constraint fix_halo_x0(Constraint_tp::STATE, 0, IC_halo, 1);

	halo.addConstraint(fix_halo_x0);
	halo.addConstraint(periodicity);

	MultShootEngine corrector;
	corrector.setVarTime(true);
	corrector.setEqualArcTime(true);
	BOOST_CHECK_NO_THROW(corrector.multShoot(&halo, &correctedHalo));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_Resonant_RevTime){
	SysData_cr3bp sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = -31.00065761;

	SimEngine sim;
	Arcset_cr3bp nodeset(&sys), correctedNodeset(&sys);
	sim.runSim_manyNodes(IC, tof, 15, &nodeset);

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

	MultShootEngine corrector;
	BOOST_CHECK_NO_THROW(corrector.multShoot(&nodeset, &correctedNodeset));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_Resonant_doubleSource){
	SysData_cr3bp sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = 31.00065761;

	Arcset_cr3bp posTimeArc(&sys), revTimeArc(&sys), correctedNodeset(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(IC, tof/2, 8, &posTimeArc);
	sim.setRevTime(true);
	sim.runSim_manyNodes(IC, -tof/2, 8, &revTimeArc);

	Arcset_cr3bp nodeset = posTimeArc;
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

	MultShootEngine corrector;
	BOOST_CHECK_NO_THROW(corrector.multShoot(&nodeset, &correctedNodeset));
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_Resonant_Irregular){
	SysData_cr3bp sys("earth", "moon");

	// ICs for a 2:5 Resonant Orbit in the EM System
	double IC[] = {0.6502418226, 0, 0, 0, 0.9609312003, 0};	
	double tof = 31.00065761;

	Arcset_cr3bp posTimeArc(&sys), revTimeArc(&sys), correctedNodeset(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(IC, tof/2, 8, &posTimeArc);
	sim.setRevTime(true);
	sim.runSim_manyNodes(IC, -tof/2, 8, &revTimeArc);

	Arcset_cr3bp nodeset = posTimeArc;
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

	MultShootEngine corrector;
	BOOST_CHECK_NO_THROW(corrector.multShoot(&nodeset, &correctedNodeset));
}//====================================================

BOOST_AUTO_TEST_CASE(BC4BP_SEM_Halo){
	SysData_bc4bp sys("Sun", "earth", "moon");
	std::vector<double> haloIC {-1.144739, 0, 0.089011, 0, 0.011608, 0};
	double tof = 310;

	Arcset_bc4bp nodeset(&sys), correctedNodeset(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(haloIC, 0, tof, 7, &nodeset);

	
	// nodeset.saveToMat("bc4bp_halo_raw.mat");
	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	Constraint perpCross(Constraint_tp::STATE, 0, perpCrossData, 6);

	double xzPlaneData[] = {NAN,0,NAN,NAN,NAN,NAN};
	Constraint xzPlaneCon0(Constraint_tp::STATE, 3, xzPlaneData, 6);
	Constraint xzPlaneConF(Constraint_tp::STATE, 6, xzPlaneData, 6);

	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(xzPlaneCon0);
	nodeset.addConstraint(xzPlaneConF);

	MultShootEngine corrector;
	BOOST_CHECK_NO_THROW(corrector.multShoot(&nodeset, &correctedNodeset));
}//====================================================

BOOST_AUTO_TEST_CASE(BC4BP_SEM_Halo_RevTime){
	SysData_bc4bp sys("Sun", "earth", "moon");
	std::vector<double> haloIC {-1.144739, 0, 0.089011, 0, 0.011608, 0};
	double tof = -310;

	Arcset_bc4bp nodeset(&sys), correctedNodeset(&sys);
	SimEngine sim;
	sim.setRevTime(true);
	sim.runSim_manyNodes(haloIC, 0, tof, 7, &nodeset);
	// nodeset.saveToMat("bc4bp_halo_raw.mat");

	double perpCrossData[] = {NAN,0,NAN,0,NAN,0};
	Constraint perpCross(Constraint_tp::STATE, 0, perpCrossData, 6);

	double xzPlaneData[] = {NAN,0,NAN,NAN,NAN,NAN};
	Constraint xzPlaneCon0(Constraint_tp::STATE, 3, xzPlaneData, 6);
	Constraint xzPlaneConF(Constraint_tp::STATE, 6, xzPlaneData, 6);

	nodeset.addConstraint(perpCross);
	nodeset.addConstraint(xzPlaneCon0);
	nodeset.addConstraint(xzPlaneConF);

	MultShootEngine corrector;
	BOOST_CHECK_NO_THROW(corrector.multShoot(&nodeset, &correctedNodeset));
}//====================================================

BOOST_AUTO_TEST_CASE(BC4BP_SEM_doubleSource){
	SysData_bc4bp sys("Sun", "earth", "moon");
	std::vector<double> haloIC {-1.144739, 0, 0.089011, 0, 0.011608, 0};
	double tof = 310;

	Arcset_bc4bp posTimeArc(&sys), revTimeArc(&sys), correctedNodeset(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(haloIC, 0, tof/2, 4, &posTimeArc);
	sim.setRevTime(true);
	sim.runSim_manyNodes(haloIC, 0, -tof/2, 4, &revTimeArc);

	Arcset_bc4bp nodeset = posTimeArc;
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

	MultShootEngine corrector;
	BOOST_CHECK_NO_THROW(corrector.multShoot(&nodeset, &correctedNodeset));
}//====================================================

BOOST_AUTO_TEST_CASE(BC4BP_SEM_Halo_DoubleSource_Irregular){
	SysData_bc4bp sys("Sun", "earth", "moon");
	std::vector<double> haloIC {-1.144739, 0, 0.089011, 0, 0.011608, 0};
	double tof = 310;

	Arcset_bc4bp posTimeArc(&sys), revTimeArc(&sys), correctedNodeset(&sys);
	SimEngine sim;
	sim.runSim_manyNodes(haloIC, 0, tof/2, 4, &posTimeArc);
	sim.setRevTime(true);
	sim.runSim_manyNodes(haloIC, 0, -tof/2, 4, &revTimeArc);

	Arcset_bc4bp nodeset = posTimeArc;
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

	MultShootEngine corrector;
	corrector.setTol(5e-12);
	BOOST_CHECK_NO_THROW(corrector.multShoot(&nodeset, &correctedNodeset));
}//====================================================
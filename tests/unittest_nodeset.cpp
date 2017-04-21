#define BOOST_TEST_MODULE Nodeset

#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iostream>

#include "AsciiOutput.hpp"
#include "Constraint.hpp"
#include "ControlLaw_cr3bp_lt.hpp"
#include "MultShootEngine.hpp"
#include "Event.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "Node.hpp"
#include "Nodeset_bc4bp.hpp"
#include "Nodeset_cr3bp.hpp"
#include "Nodeset_cr3bp_lt.hpp"
#include "SysData_bc4bp.hpp"
#include "SysData_cr3bp.hpp"
#include "SysData_cr3bp_lt.hpp"
#include "Traj_bc4bp.hpp"
#include "Utilities.hpp"

using namespace astrohelion;

Nodeset_bc4bp *bcSet;

Nodeset_cr3bp test_createCR3BPNodeset(SysData_cr3bp*);
void test_createBCR4BPRNodeset(SysData_bc4bp*);
bool dummy_predicate(Exception&);

bool dummy_predicate( Exception const& ex ) { return true; }

BOOST_AUTO_TEST_CASE(Concat_CR3BP){
	SysData_cr3bp sys("Saturn", "Titan");
	SysData_cr3bp emSys("earth", "moon");
	Nodeset_cr3bp set1(&sys);
	Nodeset_cr3bp set2(&sys);
	Nodeset_cr3bp set3(&sys);
	Nodeset_cr3bp set4(&emSys);

	double state1[] = {1,0,0,0,0,0};
	double state2[] = {2,1,0,0,0,0};
	double state3[] = {3,0,1,0,0,0};
	double state4[] = {4,0,0,1,0,0};

	set1.addNode(Node(state1, 6, 0));
	set1.addNode(Node(state2, 6, 1.1));
	set1.addSeg(Segment(0, 1, 1.1));

	set2.addNode(Node(state3, 6, 0));
	set2.addNode(Node(state4, 6, 2.2));
	set2.addSeg(Segment(0, 1, 2.2));

	set3.addNode(Node(state2, 6, 3.3));
	set3.addNode(Node(state3, 6, 4.4));
	set3.addNode(Node(state4, 6, 5.5));
	set3.addSeg(Segment(0, 1, 1.1));
	set3.addSeg(Segment(1, 2, 1.1));

	set4.addNode(Node(state4, 6, 0));

	Nodeset_cr3bp sum1 = set1 + set2;
	BOOST_CHECK(sum1.getStateByIx(0)[0] == 1);
	BOOST_CHECK(sum1.getStateByIx(1)[0] == 2);
	BOOST_CHECK(sum1.getStateByIx(2)[0] == 4);

	Nodeset_cr3bp sum2 = set1;
	sum2 += set3;
	BOOST_CHECK(sum2.getStateByIx(0)[0] == 1);
	BOOST_CHECK(sum2.getStateByIx(1)[0] == 2);
	BOOST_CHECK(sum2.getStateByIx(2)[0] == 3);
	BOOST_CHECK(sum2.getStateByIx(3)[0] == 4);

	// Sum of different systems
	Nodeset_cr3bp sum3 = set1;
	BOOST_CHECK_EXCEPTION(sum3 += set4, Exception, dummy_predicate);	
}//=======================================

BOOST_AUTO_TEST_CASE(CR3BP_NodesAtEvents){
	SysData_cr3bp sys("earth", "moon");
	double emDRO_ic[] = {0.66703088566639, 0, 0, 0, 0.763253816058075, 0};
	double emDRO_T = 5.18136624737627;

	// Second test case: Generate orbit, use createNodesAtEvent and check the functionality, TOF computation, etc.
	Nodeset_cr3bp set2(&sys, emDRO_ic, emDRO_T, 2);
	// set2.saveToMat("emDRO_2Nodes.mat");
	// cout << "CR3BP Nodeset generated from ICs (saved to emDRO_2Nodes.mat):" << endl;
	BOOST_CHECK(set2.getNumNodes() == 2);
	BOOST_CHECK(set2.getTOFByIx(0) == emDRO_T);

	std::vector<double> xMoonData = {1 - sys.getMu()};
	Event xMoonEvt(Event_tp::YZ_PLANE, 0, true, xMoonData);
	Event xzPlaneEvt(Event_tp::XZ_PLANE, 0, true);
	std::vector<Event> events {xMoonEvt, xzPlaneEvt};

	set2.createNodesAtEvents(0, events);
	set2.putInChronoOrder();
	// set2.printInChrono();
	// set2.saveToMat("emDRO_newNodes.mat");
	// cout << "CR3BP createNodesAtEvents (saved to emDRO_newNodes.mat):" << endl;
	BOOST_CHECK(set2.getNumNodes() == 5);
	BOOST_CHECK(set2.getStateByIx(1)[0] == xMoonData[0]);
	BOOST_CHECK(set2.getStateByIx(2)[1] == 0);
	BOOST_CHECK(set2.getStateByIx(3)[0] == xMoonData[0]);
	BOOST_CHECK(set2.getTotalTOF() == emDRO_T);
	// set2.print();
}//====================================================

BOOST_AUTO_TEST_CASE(BC4BP_NodesetFromICs){
	SysData_bc4bp bcSys("sun", "earth", "moon");
	double qho_ic[] = {-0.86464955943628, -0.523239865136876, -0.0309591111054232, -0.00352683110021282, -0.00217207557203108, 0.00179392516522105};
	double qho_T0 = 100;
	double qho_Period = 360;
	Nodeset_bc4bp set3(&bcSys, qho_ic, qho_T0, qho_Period, 2);
	// cout << "BC4BP Nodeset generated from ICs:" << endl;
	BOOST_CHECK(set3.getNumNodes() == 2);
	BOOST_CHECK(set3.getTOFByIx(0) == qho_Period);
	BOOST_CHECK(set3.getEpochByIx(0) == qho_T0);
	BOOST_CHECK(set3.getEpochByIx(1) == qho_T0 + qho_Period);
	// set3.print();

	Event sem_xzPlaneEvt(Event_tp::XZ_PLANE, 0, false);
	set3.createNodesAtEvent(0, sem_xzPlaneEvt);
	set3.putInChronoOrder();
	// set3.printInChrono();
	// set3.saveToMat("semQHO_newNodes.mat");
	// cout << "BC4BP createNodesAtEvent (saved to semQHO_newNodes.mat):" << endl;
	BOOST_CHECK(set3.getNumNodes() == 4);
	BOOST_CHECK(std::abs(set3.getStateByIx(1)[1]) < 1e-12);
	BOOST_CHECK(std::abs(set3.getEpochByIx(1) - qho_T0 - set3.getTOFByIx(0)) < 1e-12);
	BOOST_CHECK(std::abs(set3.getStateByIx(2)[1]) < 1e-12);
	BOOST_CHECK(std::abs(set3.getEpochByIx(2) - qho_T0 - set3.getTOFByIx(0)) < 1e-12 + set3.getTOFByIx(1));
	BOOST_CHECK(set3.getTotalTOF() == qho_Period);
	// set3.print();
}//====================================================

Nodeset_cr3bp test_createCR3BPNodeset(SysData_cr3bp *emData){
	// Define system and IC
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	// Create a node set from the IC and sysDdata
	Nodeset_cr3bp crSet(emData, ic, 2.77, 5, Nodeset::TIME);
	
	int nodes[] = {3,4};
	std::vector<int> velCon(nodes, nodes+2);
	// crSet->allowDV_at(velCon);

	// Add a constraint
	// double data[] = {1,1,1,NAN,NAN,NAN};
	double data[] = {1.5};
	Constraint crCon1(Constraint_tp::MAX_DELTA_V, 0, data, 1);
	// crSet->addConstraint(crCon1);

	// crSet.print();

	return crSet;
}//====================================================

void test_createBCR4BPRNodeset(SysData_bc4bp *semData){
	double ic2[] = {82.575887, 0, 8.0, 0, 0.19369725, 0};

	bcSet = new Nodeset_bc4bp(semData, ic2, 0, 40, 5, Nodeset::TIME);

	// Add a constraint
	// double data[] = {82.576, 0, 8.001, NAN, NAN, NAN, NAN};
	// double data[] = {5,5,5,NAN,NAN,NAN,NAN};
	double data[] = {1.5};
	Constraint bcCon1(Constraint_tp::MAX_DELTA_V, 0, data, 1);
	bcSet->addConstraint(bcCon1);

	int nodes[] = {2,3};
	std::vector<int> velCon(nodes, nodes+2);
	bcSet->allowDV_at(velCon);

	// bcSet->print();
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_Nodeset_Save_Load){
	MultShootEngine corrector;

	SysData_cr3bp emData("earth", "moon");
	Nodeset_cr3bp crSet = test_createCR3BPNodeset(&emData);
	// corrector.setVerbosity(Verbosity_tp::ALL_MSG);
	corrector.multShoot(&crSet, NULL);

	crSet.saveToMat("data/crSet.mat");
	Nodeset_cr3bp crTemp(&emData);
	crTemp.readFromMat("data/crSet.mat");

	BOOST_CHECK(crSet.getStateByIx(-1) == crTemp.getStateByIx(-1));
	BOOST_CHECK(crSet.getTOFByIx(-1) == crTemp.getTOFByIx(-1));
}//====================================================

BOOST_AUTO_TEST_CASE(BC4BP_Nodeset_Save_Load){
	MultShootEngine corrector;
	
	SysData_bc4bp semData("sun", "earth", "moon");
	test_createBCR4BPRNodeset(&semData);
	// corrector.setVerbosity(Verbosity_tp::ALL_MSG);
	corrector.multShoot(bcSet, NULL);

	bcSet->saveToMat("data/bcSet.mat");
	Nodeset_bc4bp bcTemp(&semData);
	bcTemp.readFromMat("data/bcSet.mat");

	BOOST_CHECK(bcSet->getStateByIx(-1) == bcTemp.getStateByIx(-1));
	BOOST_CHECK(bcSet->getTOFByIx(-1) == bcTemp.getTOFByIx(-1));
	BOOST_CHECK(bcSet->getEpochByIx(-1) == bcTemp.getEpochByIx(-1));

	delete bcSet;
}//====================================================

BOOST_AUTO_TEST_CASE(CR3BP_LT_Nodeset_Save_Load){
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0, 1};

	MultShootEngine corrector;

	SysData_cr3bp_lt ltData("earth", "moon", 12e-3, 1500, 14);
	Nodeset_cr3bp_lt ltSet(&ltData, ic, 2.77, 5, Nodeset::TIME, ControlLaw_cr3bp_lt::Law_tp::CONST_C_2D_LEFT);

	corrector.multShoot(&ltSet, NULL);

	ltSet.saveToMat("data/ltSet.mat");
	Nodeset_cr3bp_lt temp(&ltData);
	temp.readFromMat("data/ltSet.mat");

	BOOST_CHECK(ltSet.getStateByIx(-1) == temp.getStateByIx(-1));
	BOOST_CHECK(temp.getCtrlLawByIx(0) == ltSet.getCtrlLawByIx(0));
	BOOST_CHECK(ltSet.getCtrlLawByIx(0) == ControlLaw_cr3bp_lt::Law_tp::CONST_C_2D_LEFT);
	BOOST_CHECK(ltSet.getTOFByIx(-1) == temp.getTOFByIx(-1));
}//====================================================






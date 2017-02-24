/**
 *	Test the nodeset object
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of Astrohelion
 *
 *  Astrohelion is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Astrohelion is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Astrohelion.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "AsciiOutput.hpp"
#include "Constraint.hpp"
#include "CorrectionEngine.hpp"
#include "Event.hpp"
#include "Exceptions.hpp"
#include "MultShootData.hpp"
#include "Node.hpp"
#include "Nodeset_cr3bp.hpp"
#include "Nodeset_bc4bp.hpp"
#include "SysData_bc4bp.hpp"
#include "SysData_cr3bp.hpp"
#include "Traj_bc4bp.hpp"
#include "Utilities.hpp"

#include <cmath>
#include <iostream>

using namespace std;
using namespace astrohelion;

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

Nodeset_cr3bp *crSet;
Nodeset_bc4bp *bcSet;

void test_concat_CR3BP(){
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
	bool checkSum1 = sum1.getStateByIx(0)[0] == 1 && sum1.getStateByIx(1)[0] == 2 && sum1.getStateByIx(2)[0] == 4;
	cout << "CR3BP Nodeset operator+ : " << (checkSum1 ? PASS : FAIL) << endl;
	// sum1.print();

	Nodeset_cr3bp sum2 = set1;
	sum2 += set3;
	bool checkSum2 = sum2.getStateByIx(0)[0] == 1 && sum2.getStateByIx(1)[0] == 2 && 
		sum2.getStateByIx(2)[0] == 3 && sum2.getStateByIx(3)[0] == 4;
	cout << "CR3BP Nodeset operator+= : " << (checkSum2 ? PASS : FAIL) << endl;
	// sum2.print();
	
	try{
		cout << "Testing sum of different systems: ";
		Nodeset_cr3bp sum3 = set1;
		sum3 += set4;
		cout << FAIL << endl;
	}catch(Exception &e){
		cout << PASS << endl;
	}catch(exception &e){
		cout << FAIL << endl;
	}
}//=======================================

void test_nodeManip(){
	SysData_cr3bp sys("earth", "moon");
	Nodeset_cr3bp set(&sys);
	double emDRO_ic[] = {0.66703088566639, 0, 0, 0, 0.763253816058075, 0};
	double emDRO_T = 5.18136624737627;

	double node1[] = {1,0,0,0,0,0};
	double node2[] = {2,0,0,0,0,0};
	double node3[] = {3,0,0,0,0,0};
	double node4[] = {4,0,0,0,0,0};

	set.addNode(Node(node1, 6, 0.1));
	set.addNode(Node(node2, 6, 0.2));
	set.addNode(Node(node3, 6, 0.3));
	set.addNode(Node(node4, 6, 0.4));

	// Second test case: Generate orbit, use createNodesAtEvent and check the functionality, TOF computation, etc.
	Nodeset_cr3bp set2(&sys, emDRO_ic, emDRO_T, 2);
	set2.saveToMat("emDRO_2Nodes.mat");
	cout << "CR3BP Nodeset generated from ICs (saved to emDRO_2Nodes.mat):" << endl;
	cout << "  Correct number of nodes: " << (set2.getNumNodes() == 2 ? PASS : FAIL) << endl;
	cout << "  Correct TOFs: " << (set2.getTOFByIx(0) == emDRO_T ? PASS : FAIL) << endl;	

	std::vector<double> xMoonData = {1 - sys.getMu()};
	Event xMoonEvt(Event_tp::YZ_PLANE, 0, true, xMoonData);
	Event xzPlaneEvt(Event_tp::XZ_PLANE, 0, true);
	std::vector<Event> events {xMoonEvt, xzPlaneEvt};

	set2.createNodesAtEvents(0, events);
	set2.putInChronoOrder();
	// set2.printInChrono();
	set2.saveToMat("emDRO_newNodes.mat");
	cout << "CR3BP createNodesAtEvents (saved to emDRO_newNodes.mat):" << endl;
	cout << "  Correct number of nodes: " << (set2.getNumNodes() == 5 ? PASS : FAIL) << endl;
	cout << "  Correct node(1) state: " << (set2.getStateByIx(1)[0] == xMoonData[0] ? PASS : FAIL) << endl;
	cout << "  Correct node(2) state: " << (set2.getStateByIx(2)[1] == 0 ? PASS : FAIL) << endl;
	cout << "  Correct node(3) state: " << (set2.getStateByIx(3)[0] == xMoonData[0] ? PASS : FAIL) << endl;
	cout << "  Correct total TOF: " << (set2.getTotalTOF() == emDRO_T ? PASS : FAIL) << endl;
	set2.print();

	SysData_bc4bp bcSys("sun", "earth", "moon");
	double qho_ic[] = {-0.86464955943628, -0.523239865136876, -0.0309591111054232, -0.00352683110021282, -0.00217207557203108, 0.00179392516522105};
	double qho_T0 = 100;
	double qho_Period = 360;
	Nodeset_bc4bp set3(&bcSys, qho_ic, qho_T0, qho_Period, 2);
	cout << "BC4BP Nodeset generated from ICs:" << endl;
	cout << "  Correct number of nodes: " << (set3.getNumNodes() == 2 ? PASS : FAIL) << endl;
	cout << "  Correct TOFs: " << (set3.getTOFByIx(0) == qho_Period ? PASS : FAIL) << endl;
	cout << "  Correct Epochs: " << (set3.getEpochByIx(0) == qho_T0 && set3.getEpochByIx(1) == qho_T0 + qho_Period ? PASS : FAIL) << endl;
	set3.print();

	Event sem_xzPlaneEvt(Event_tp::XZ_PLANE, 0, false);
	set3.createNodesAtEvent(0, sem_xzPlaneEvt);
	set3.putInChronoOrder();
	// set3.printInChrono();
	set3.saveToMat("semQHO_newNodes.mat");
	cout << "BC4BP createNodesAtEvent (saved to semQHO_newNodes.mat):" << endl;
	cout << "  Correct number of nodes: " << (set3.getNumNodes() == 4 ? PASS : FAIL) << endl;
	cout << "  Correct node(1) state: " << (std::abs(set3.getStateByIx(1)[1]) < 1e-12 ? PASS : FAIL) << endl;
	cout << "  Correct node(1) epoch: " << (std::abs(set3.getEpochByIx(1) - qho_T0 - set3.getTOFByIx(0)) < 1e-12 ? PASS : FAIL) << endl;
	cout << "  Correct node(2) state: " << (std::abs(set3.getStateByIx(2)[1]) < 1e-12 ? PASS : FAIL) << endl;
	cout << "  Correct node(2) epoch: " << (std::abs(set3.getEpochByIx(2) - qho_T0 - set3.getTOFByIx(0)) < 1e-12 + set3.getTOFByIx(1) ? PASS : FAIL) << endl;
	cout << "  Correct total TOF: " << (set3.getTotalTOF() == qho_Period ? PASS : FAIL) << endl;
	cout << "Total TOF = " << set3.getTotalTOF() << endl;
	set3.print();
}//==============================================

Nodeset_cr3bp test_createCR3BPNodeset(SysData_cr3bp *emData){
	// Define system and IC
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	// Create a node set from the IC and sysDdata
	Nodeset_cr3bp crSet(emData, ic, 2.77, 5, Nodeset::TIME);
	
	int nodes[] = {3,4};
	vector<int> velCon(nodes, nodes+2);
	// crSet->allowDV_at(velCon);

	// Add a constraint
	// double data[] = {1,1,1,NAN,NAN,NAN};
	double data[] = {1.5};
	Constraint crCon1(Constraint_tp::MAX_DELTA_V, 0, data, 1);
	// crSet->addConstraint(crCon1);

	crSet.print();

	return crSet;
}//==============================================

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
	vector<int> velCon(nodes, nodes+2);
	bcSet->allowDV_at(velCon);

	bcSet->print();
}//==============================================

int main(void){
	CorrectionEngine corrector;

	SysData_cr3bp emData("earth", "moon");
	Nodeset_cr3bp crSet = test_createCR3BPNodeset(&emData);
	// corrector.setVerbosity(Verbosity_tp::ALL_MSG);
	corrector.multShoot(&crSet, NULL);

	crSet.saveToMat("data/crSet.mat");
	Nodeset_cr3bp crTemp(&emData);
	crTemp.readFromMat("data/crSet.mat");
	
	SysData_bc4bp semData("sun", "earth", "moon");
	test_createBCR4BPRNodeset(&semData);
	// corrector.setVerbosity(Verbosity_tp::ALL_MSG);
	corrector.multShoot(bcSet, NULL);

	bcSet->saveToMat("data/bcSet.mat");
	Nodeset_bc4bp bcTemp(&semData);
	bcTemp.readFromMat("data/bcSet.mat");

	test_concat_CR3BP();

	printf("Testing Save/Read functions on CR3BP Nodeset\n");
	cout << "Same Final State: " << (crSet.getStateByIx(-1) == crTemp.getStateByIx(-1) ? PASS : FAIL) << endl;
	cout << "Same Final TOF: " << (crSet.getTOFByIx(-1) == crTemp.getTOFByIx(-1) ? PASS : FAIL) << endl;

	printf("Testing Save/Read functions on BC4BP Nodeset\n");
	cout << "Same Final State: " << (bcSet->getStateByIx(-1) == bcTemp.getStateByIx(-1) ? PASS : FAIL) << endl;
	cout << "Same Final TOF: " << (bcSet->getTOFByIx(-1) == bcTemp.getTOFByIx(-1) ? PASS : FAIL) << endl;
	cout << "Same Final Epoch: " << (bcSet->getEpochByIx(-1) == bcTemp.getEpochByIx(-1) ? PASS : FAIL) << endl;
	
	printf("Testing Node Insert/Delete/InsertAtEvent\n");
	test_nodeManip();

	// Memory clean-up
	delete bcSet;
	return 0;
}
/**
 *	Test the nodeset object
 */
/*
 *	Trajectory Propagation and Analysis Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
 *
 *  TPAT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  TPAT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with TPAT.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "tpat_ascii_output.hpp"
#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_event.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_multShoot_data.hpp"
#include "tpat_node.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_nodeset_bcr4bp.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_bcr4bp.hpp"
#include "tpat_utilities.hpp"

#include <cmath>
#include <iostream>

using namespace std;

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

tpat_nodeset_cr3bp *crSet;
tpat_nodeset_bcr4bp *bcSet;

void test_concat_CR3BP(){
	tpat_sys_data_cr3bp sys("Saturn", "Titan");
	tpat_sys_data_cr3bp emSys("earth", "moon");
	tpat_nodeset_cr3bp set1(&sys);
	tpat_nodeset_cr3bp set2(&sys);
	tpat_nodeset_cr3bp set3(&sys);
	tpat_nodeset_cr3bp set4(&emSys);

	double state1[] = {1,0,0,0,0,0};
	double state2[] = {2,1,0,0,0,0};
	double state3[] = {3,0,1,0,0,0};
	double state4[] = {4,0,0,1,0,0};

	set1.addNode(tpat_node(state1, 0));
	set1.addNode(tpat_node(state2, 1.1));
	set1.addSeg(tpat_segment(0, 1, 1.1));

	set2.addNode(tpat_node(state3, 0));
	set2.addNode(tpat_node(state4, 2.2));
	set2.addSeg(tpat_segment(0, 1, 2.2));

	set3.addNode(tpat_node(state2, 3.3));
	set3.addNode(tpat_node(state3, 4.4));
	set3.addNode(tpat_node(state4, 5.5));
	set3.addSeg(tpat_segment(0, 1, 1.1));
	set3.addSeg(tpat_segment(1, 2, 1.1));

	set4.addNode(tpat_node(state4, 0));

	tpat_nodeset_cr3bp sum1 = set1 + set2;
	bool checkSum1 = sum1.getStateByIx(0)[0] == 1 && sum1.getStateByIx(1)[0] == 2 && sum1.getStateByIx(2)[0] == 4;
	cout << "CR3BP Nodeset operator+ : " << (checkSum1 ? PASS : FAIL) << endl;
	// sum1.print();

	tpat_nodeset_cr3bp sum2 = set1;
	sum2 += set3;
	bool checkSum2 = sum2.getStateByIx(0)[0] == 1 && sum2.getStateByIx(1)[0] == 2 && 
		sum2.getStateByIx(2)[0] == 3 && sum2.getStateByIx(3)[0] == 4;
	cout << "CR3BP Nodeset operator+= : " << (checkSum2 ? PASS : FAIL) << endl;
	// sum2.print();
	
	try{
		cout << "Testing sum of different systems: ";
		tpat_nodeset_cr3bp sum3 = set1;
		sum3 += set4;
		cout << FAIL << endl;
	}catch(tpat_exception &e){
		cout << PASS << endl;
	}catch(exception &e){
		cout << FAIL << endl;
	}
}//=======================================

void test_nodeManip(){
	tpat_sys_data_cr3bp sys("earth", "moon");
	tpat_nodeset_cr3bp set(&sys);
	double emDRO_ic[] = {0.66703088566639, 0, 0, 0, 0.763253816058075, 0};
	double emDRO_T = 5.18136624737627;

	double node1[] = {1,0,0,0,0,0};
	double node2[] = {2,0,0,0,0,0};
	double node3[] = {3,0,0,0,0,0};
	double node4[] = {4,0,0,0,0,0};

	set.addNode(tpat_node(node1, 0.1));
	set.addNode(tpat_node(node2, 0.2));
	set.addNode(tpat_node(node3, 0.3));
	set.addNode(tpat_node(node4, 0.4));

	// Second test case: Generate orbit, use createNodesAtEvent and check the functionality, TOF computation, etc.
	tpat_nodeset_cr3bp set2(emDRO_ic, &sys, emDRO_T, 2);
	set2.saveToMat("emDRO_2Nodes.mat");
	cout << "CR3BP Nodeset generated from ICs (saved to emDRO_2Nodes.mat):" << endl;
	cout << "  Correct number of nodes: " << (set2.getNumNodes() == 2 ? PASS : FAIL) << endl;
	cout << "  Correct TOFs: " << (set2.getTOFByIx(0) == emDRO_T ? PASS : FAIL) << endl;	

	double xMoonData = 1 - sys.getMu();
	tpat_event xMoonEvt(&sys, tpat_event::YZ_PLANE, 0, true, &xMoonData);
	tpat_event xzPlaneEvt(&sys, tpat_event::XZ_PLANE, 0, true);
	std::vector<tpat_event> events {xMoonEvt, xzPlaneEvt};

	set2.createNodesAtEvents(0, events);
	set2.putInChronoOrder();
	// set2.printInChrono();
	set2.saveToMat("emDRO_newNodes.mat");
	cout << "CR3BP createNodesAtEvents (saved to emDRO_newNodes.mat):" << endl;
	cout << "  Correct number of nodes: " << (set2.getNumNodes() == 5 ? PASS : FAIL) << endl;
	cout << "  Correct node(1) state: " << (set2.getStateByIx(1)[0] == xMoonData ? PASS : FAIL) << endl;
	cout << "  Correct node(2) state: " << (set2.getStateByIx(2)[1] == 0 ? PASS : FAIL) << endl;
	cout << "  Correct node(3) state: " << (set2.getStateByIx(3)[0] == xMoonData ? PASS : FAIL) << endl;
	cout << "  Correct total TOF: " << (set2.getTotalTOF() == emDRO_T ? PASS : FAIL) << endl;
	set2.print();

	tpat_sys_data_bcr4bpr bcSys("sun", "earth", "moon");
	double qho_ic[] = {-0.86464955943628, -0.523239865136876, -0.0309591111054232, -0.00352683110021282, -0.00217207557203108, 0.00179392516522105};
	double qho_T0 = 100;
	double qho_Period = 360;
	tpat_nodeset_bcr4bp set3(qho_ic, &bcSys, qho_T0, qho_Period, 2);
	cout << "BC4BP Nodeset generated from ICs:" << endl;
	cout << "  Correct number of nodes: " << (set3.getNumNodes() == 2 ? PASS : FAIL) << endl;
	cout << "  Correct TOFs: " << (set3.getTOFByIx(0) == qho_Period ? PASS : FAIL) << endl;
	cout << "  Correct Epochs: " << (set3.getEpochByIx(0) == qho_T0 && set3.getEpochByIx(1) == qho_T0 + qho_Period ? PASS : FAIL) << endl;
	set3.print();

	tpat_event sem_xzPlaneEvt(&bcSys, tpat_event::XZ_PLANE, 0, false);
	set3.createNodesAtEvent(0, sem_xzPlaneEvt);
	set3.putInChronoOrder();
	// set3.printInChrono();
	set3.saveToMat("semQHO_newNodes.mat");
	cout << "BC4BP createNodesAtEvent (saved to semQHO_newNodes.mat):" << endl;
	cout << "  Correct number of nodes: " << (set3.getNumNodes() == 4 ? PASS : FAIL) << endl;
	cout << "  Correct node(1) state: " << (set3.getStateByIx(1)[1] == 0 ? PASS : FAIL) << endl;
	cout << "  Correct node(1) epoch: " << (set3.getEpochByIx(1) == qho_T0 + set3.getTOFByIx(0) ? PASS : FAIL) << endl;
	cout << "  Correct node(2) state: " << (set3.getStateByIx(2)[1] == 0 ? PASS : FAIL) << endl;
	cout << "  Correct node(2) epoch: " << (set3.getEpochByIx(2) == qho_T0 + set3.getTOFByIx(0) + set3.getTOFByIx(1) ? PASS : FAIL) << endl;
	cout << "  Correct total TOF: " << (set3.getTotalTOF() == qho_Period ? PASS : FAIL) << endl;
	cout << "Total TOF = " << set3.getTotalTOF() << endl;
	set3.print();
}//==============================================

tpat_nodeset_cr3bp test_createCR3BPNodeset(tpat_sys_data_cr3bp *emData){
	// Define system and IC
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	// Create a node set from the IC and sysDdata
	tpat_nodeset_cr3bp crSet(ic, emData, 2.77, 5, tpat_nodeset::DISTRO_TIME);
	
	int nodes[] = {3,4};
	vector<int> velCon(nodes, nodes+2);
	// crSet->allowDV_at(velCon);

	// Add a constraint
	// double data[] = {1,1,1,NAN,NAN,NAN};
	double data[] = {1.5};
	tpat_constraint crCon1(tpat_constraint::MAX_DELTA_V, 0, data, 1);
	// crSet->addConstraint(crCon1);

	crSet.print();

	return crSet;
}//==============================================

void test_createBCR4BPRNodeset(tpat_sys_data_bcr4bpr *semData){
	double ic2[] = {82.575887, 0, 8.0, 0, 0.19369725, 0};

	bcSet = new tpat_nodeset_bcr4bp(ic2, semData, 0, 40, 5, tpat_nodeset::DISTRO_TIME);

	// Add a constraint
	// double data[] = {82.576, 0, 8.001, NAN, NAN, NAN, NAN};
	// double data[] = {5,5,5,NAN,NAN,NAN,NAN};
	double data[] = {1.5};
	tpat_constraint bcCon1(tpat_constraint::MAX_DELTA_V, 0, data, 1);
	bcSet->addConstraint(bcCon1);

	int nodes[] = {2,3};
	vector<int> velCon(nodes, nodes+2);
	bcSet->allowDV_at(velCon);

	bcSet->print();
}//==============================================

int main(void){
	tpat_correction_engine corrector;

	tpat_sys_data_cr3bp emData("earth", "moon");
	tpat_nodeset_cr3bp crSet = test_createCR3BPNodeset(&emData);
	// corrector.setVerbose(ALL_MSG);
	corrector.multShoot(&crSet, NULL);

	crSet.saveToMat("data/crSet.mat");
	tpat_nodeset_cr3bp crTemp(&emData);
	crTemp.readFromMat("data/crSet.mat");
	
	tpat_sys_data_bcr4bpr semData("sun", "earth", "moon");
	test_createBCR4BPRNodeset(&semData);
	// corrector.setVerbose(ALL_MSG);
	corrector.multShoot(bcSet, NULL);

	bcSet->saveToMat("data/bcSet.mat");
	tpat_nodeset_bcr4bp bcTemp(&semData);
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
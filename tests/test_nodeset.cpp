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
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_nodeset_bcr4bpr.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_exceptions.hpp"

#include <cmath>
#include <iostream>

using namespace std;

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

tpat_nodeset_cr3bp *crSet;
tpat_nodeset_bcr4bpr *bcSet;

void test_concat_CR3BP(){
	tpat_sys_data_cr3bp sys("Saturn", "Titan");
	tpat_sys_data_cr3bp emSys("earth", "moon");
	tpat_nodeset_cr3bp set1(sys);
	tpat_nodeset_cr3bp set2(sys);
	tpat_nodeset_cr3bp set3(sys);
	tpat_nodeset_cr3bp set4(emSys);

	double node1[] = {1,0,0,0,0,0};
	double node2[] = {2,1,0,0,0,0};
	double node3[] = {3,0,1,0,0,0};
	double node4[] = {4,0,0,1,0,0};

	set1.appendNode(node1);
	set1.appendNode(node2);
	
	set2.appendNode(node3);
	set2.appendNode(node4);
	
	set3.appendNode(node2);
	set3.appendNode(node3);
	set3.appendNode(node4);

	set4.appendNode(node4);

	tpat_nodeset_cr3bp sum1 = set1 + set2;
	printf("Concat nodesets: Should have nodes with x from 1 to 4\n");
	sum1.print();

	tpat_nodeset_cr3bp sum2 = set1 + set3;
	printf("Concat nodesets: Should have nodes with x from 1 to 4\n");
	sum2.print();
	
	try{
		cout << "Testing sum of different systems: ";
		tpat_nodeset_cr3bp sum3 = set1 + set4;
		cout << FAIL << endl;
	}catch(tpat_exception &e){
		cout << PASS << endl;
	}
	catch(...){
		cout << FAIL << endl;
	}
}

void test_createCR3BPNodeset(){
	// Define system and IC
	tpat_sys_data_cr3bp sysData("earth", "moon");
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	// Create a node set from the IC and sysDdata
	crSet = new tpat_nodeset_cr3bp(ic, sysData, 2.77, 5, tpat_nodeset::DISTRO_TIME);
	
	int nodes[] = {3,4};
	vector<int> velCon(nodes, nodes+2);
	// crSet->setVelConNodes_allBut(velCon);

	// Add a constraint
	// double data[] = {1,1,1,NAN,NAN,NAN};
	double data[] = {1.5};
	tpat_constraint crCon1(tpat_constraint::MAX_DELTA_V, 0, data, 1);
	// crSet->addConstraint(crCon1);

	printf("CR3BP Nodeset:\n Nodes: %d\n", crSet->getNumNodes());
	for (int n = 0; n < crSet->getNumNodes(); n++){
		vector<double> node = crSet->getNode(n);
		printf("  %02d: %13.8f %13.8f %13.8f %13.8f %13.8f %13.8f", n+1,
			node.at(0), node.at(1), node.at(2), node.at(3), node.at(4), node.at(5));
		if(n < crSet->getNumNodes()-1){
			printf("   TOF = %.8f\n", crSet->getTOF(n));
		}else{
			printf("\n");
		}
	}
	cout << " Constraints:" << endl;
	for(int n = 0; n < crSet->getNumCons(); n++){
		crSet->getConstraint(n).print();
	}
}//==============================================

void test_createBCR4BPRNodeset(){
	tpat_sys_data_bcr4bpr semData("sun", "earth", "moon");
	double ic2[] = {82.575887, 0, 8.0, 0, 0.19369725, 0};

	bcSet = new tpat_nodeset_bcr4bpr(ic2, semData, 0, 40, 5, tpat_nodeset::DISTRO_TIME);

	// Add a constraint
	// double data[] = {82.576, 0, 8.001, NAN, NAN, NAN, NAN};
	// double data[] = {5,5,5,NAN,NAN,NAN,NAN};
	double data[] = {1.5};
	tpat_constraint bcCon1(tpat_constraint::MAX_DELTA_V, 0, data, 1);
	bcSet->addConstraint(bcCon1);

	int nodes[] = {2,3};
	vector<int> velCon(nodes, nodes+2);
	bcSet->setVelConNodes_allBut(velCon);

	printf("BCR4BPR Nodeset:\n Nodes: %d\n", bcSet->getNumNodes());
	for (int n = 0; n < bcSet->getNumNodes(); n++){
		vector<double> node = bcSet->getNode(n);
		printf("  %02d: %9.4f -- %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n", n+1,
			bcSet->getEpoch(n), node.at(0), node.at(1), node.at(2), node.at(3),
			node.at(4), node.at(5));
	}
	cout << " Constraints:" << endl;
	for(int n = 0; n < bcSet->getNumCons(); n++){
		bcSet->getConstraint(n).print();
	}
}//==============================================


int main(void){
	tpat_correction_engine corrector;

	test_createCR3BPNodeset();
	corrector.setVerbose(true);
	corrector.correct_cr3bp(crSet);

	test_createBCR4BPRNodeset();
	corrector.setVerbose(true);
	corrector.correct_bcr4bpr(bcSet);

	// Memory clean-up
	delete crSet;
	delete bcSet;

	test_concat_CR3BP();

	return 0;
}
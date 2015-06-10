/**
 *	Test the nodeset object
 */
/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (TPAT).
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

#include "tpat_constraint.hpp"
#include "tpat_correction_engine.hpp"
#include "tpat_cr3bp_nodeset.hpp"
#include "tpat_bcr4bpr_nodeset.hpp"
#include "tpat_cr3bp_sys_data.hpp"

#include <cmath>
#include <iostream>

using namespace std;

tpat_cr3bp_nodeset *crSet;
tpat_bcr4bpr_nodeset *bcSet;

void test_createCR3BPNodeset(){
	// Define system and IC
	tpat_cr3bp_sys_data sysData("earth", "moon");
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	// Create a node set from the IC and sysDdata
	crSet = new tpat_cr3bp_nodeset(ic, sysData, 2.77, 5, tpat_nodeset::TIME);
	
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
	tpat_bcr4bpr_sys_data semData("sun", "earth", "moon");
	double ic2[] = {82.575887, 0, 8.0, 0, 0.19369725, 0};

	bcSet = new tpat_bcr4bpr_nodeset(ic2, semData, 0, 40, 5, tpat_nodeset::TIME);

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
	corrector.setVerbose(false);
	corrector.correct_cr3bp(crSet);

	test_createBCR4BPRNodeset();
	corrector.setVerbose(false);
	corrector.correct_bcr4bpr(bcSet);

	// Memory clean-up
	delete crSet;
	delete bcSet;

	return 0;
}
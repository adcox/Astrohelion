/**
 *	Test the nodeset object
 */
/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (ADTK).
 *
 *  ADTK is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ADTK is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ATDK.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "adtk_correction_engine.hpp"
#include "adtk_cr3bp_constraint.hpp"
#include "adtk_bcr4bpr_constraint.hpp"
#include "adtk_cr3bp_nodeset.hpp"
#include "adtk_bcr4bpr_nodeset.hpp"
#include "adtk_cr3bp_sys_data.hpp"

#include <cmath>
#include <iostream>

using namespace std;

adtk_cr3bp_nodeset *crSet;
adtk_bcr4bpr_nodeset *bcSet;

void test_createCR3BPNodeset(){
	// Define system and IC
	adtk_cr3bp_sys_data sysData("earth", "moon");
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	// Create a node set from the IC and sysDdata
	crSet = new adtk_cr3bp_nodeset(ic, sysData, 2.77, 5, adtk_nodeset::TIME);
	
	int nodes[] = {3,4};
	vector<int> velCon(nodes, nodes+2);
	// crSet->setVelConNodes_allBut(velCon);

	// Add a constraint
	// double data[] = {1,1,1,NAN,NAN,NAN};
	double data[] = {1.5};
	adtk_cr3bp_constraint crCon1(adtk_constraint::MAX_DELTA_V, 0, data);
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
	adtk_bcr4bpr_sys_data semData("sun", "earth", "moon");
	double ic2[] = {82.575887, 0, 8.0, 0, 0.19369725, 0};

	bcSet = new adtk_bcr4bpr_nodeset(ic2, semData, 0, 40, 5, adtk_nodeset::TIME);

	// Add a constraint
	// double data[] = {82.576, 0, 8.001, NAN, NAN, NAN, NAN};
	// double data[] = {5,5,5,NAN,NAN,NAN,NAN};
	double data[] = {1.5, NAN, NAN, NAN, NAN, NAN, NAN};
	adtk_bcr4bpr_constraint bcCon1(adtk_constraint::MAX_DELTA_V, 0, data);
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
	
	// test_createCR3BPNodeset();
	test_createBCR4BPRNodeset();

	adtk_correction_engine corrector;
	// corrector.correct_cr3bp(crSet);

	corrector.correct_bcr4bpr(bcSet);

	// Memory clean-up
	delete crSet;
	delete bcSet;

	return 0;
}
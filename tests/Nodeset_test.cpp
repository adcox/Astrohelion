/**
 *	Test the nodeset object
 */

#include "adtk_cr3bp_constraint.hpp"
#include "adtk_bcr4bpr_constraint.hpp"
#include "adtk_cr3bp_nodeset.hpp"
#include "adtk_bcr4bpr_nodeset.hpp"
#include "adtk_cr3bp_sys_data.hpp"

#include <cmath>
#include <iostream>

using namespace std;

int main(void){
	// Define system and IC
	adtk_cr3bp_sys_data sysData("earth", "moon");
	double ic[] = {0.82575887, 0, 0.08, 0, 0.19369725, 0};

	// Create a node set from the IC and sysDdata
	adtk_cr3bp_nodeset crSet(ic, sysData, 2.77, 5, adtk_nodeset::TIME);
	
	// Add a constraint
	double data[] = {1,1,1,NAN,NAN,NAN};
	adtk_cr3bp_constraint crCon1(adtk_constraint::MATCH_ALL, 3, data);
	crSet.addConstraint(crCon1);

	adtk_bcr4bpr_nodeset bcSet;

	cout << "CR3BP Node set has nodes with " << crSet.getNodeSize() << " states"<<endl;
	crSet.getConstraint(0).print();
	cout << "BCR4BPR Node set has nodes with " << bcSet.getNodeSize() << " states"<<endl;
	return 0;
}
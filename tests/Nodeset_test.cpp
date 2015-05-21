/**
 *	Test the nodeset object
 */

#include "adtk_cr3bp_constraint.hpp"
#include "adtk_bcr4bpr_constraint.hpp"
#include "adtk_cr3bp_nodeset.hpp"
#include "adtk_bcr4bpr_nodeset.hpp"

#include <cmath>
#include <iostream>

using namespace std;

int main(void){

	adtk_cr3bp_nodeset crSet;
	double data[] = {1,1,1,NAN,NAN,NAN};
	// vector<double> crCon1_data (data, data + sizeof(data)/sizeof(data[0]));
	adtk_cr3bp_constraint crCon1(adtk_constraint::MATCH_ALL, 3, data);
	crSet.addConstraint(crCon1);

	adtk_bcr4bpr_nodeset bcSet;

	cout << "CR3BP Node set has nodes with " << crSet.getNodeSize() << " states"<<endl;
	crSet.getConstraint(0).print();
	cout << "BCR4BPR Node set has nodes with " << bcSet.getNodeSize() << " states"<<endl;
	return 0;
}
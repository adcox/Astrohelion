/**
 *	Test the nodeset object
 */

#include "adtk_cr3bp_nodeset.hpp"
#include "adtk_bcr4bpr_nodeset.hpp"

#include <iostream>

using namespace std;

int main(void){

	adtk_cr3bp_nodeset crSet;
	adtk_bcr4bpr_nodeset bcSet;

	cout << "CR3BP Node set has nodes with " << crSet.getNodeSize() << " states"<<endl;
	cout << "BCR4BPR Node set has nodes with " << bcSet.getNodeSize() << " states"<<endl;
	return 0;
}
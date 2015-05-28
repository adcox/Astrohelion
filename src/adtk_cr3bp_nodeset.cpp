/**
 *	@file adtk_cr3bp_nodeset.cpp
 */

#include "adtk_cr3bp_nodeset.hpp"

#include "adtk_cr3bp_sys_data.hpp"
 
#include <cmath>
#include <iostream>

/**
 *	Create a general CR3BP nodeset
 */
adtk_cr3bp_nodeset::adtk_cr3bp_nodeset() : adtk_nodeset(6){}

/**
 *	Create a nodeset with specified system data
 *	@param data system data object
 */
adtk_cr3bp_nodeset::adtk_cr3bp_nodeset(adtk_cr3bp_sys_data data) : adtk_nodeset(6){
	sysData = data;
}

/**
 *	Compute a set of nodes by integrating from initial conditions for some time, then split the
 *	integrated trajectory into pieces (nodes).
 *
 *	@param IC a set of initial conditions, non-dimensional units
 *	@param data a pointer to a system data object that describes the model to integrate in
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC
 *	@param type node distribution type
 */
adtk_cr3bp_nodeset::adtk_cr3bp_nodeset(double IC[6], adtk_cr3bp_sys_data data, double tof, int numNodes, 
		node_distro_t type) : adtk_nodeset(6){
	sysData = data;

	initSetFromICs(IC, &sysData, 0, tof, numNodes, type);
}//======================================================================

/**
 *	Copy input nodeset
 *	@param n a nodeset
 */
adtk_cr3bp_nodeset::adtk_cr3bp_nodeset(const adtk_cr3bp_nodeset& n) : adtk_nodeset(6){
	sysData = n.sysData;
}

/**
 *	Make this nodeset equal to the input nodeset
 *	@param n a nodeset
 */
adtk_cr3bp_nodeset& adtk_cr3bp_nodeset::operator =(const adtk_cr3bp_nodeset& n){
	adtk_nodeset::operator =(n);
	sysData = n.sysData;
	return *this;
}

/**
 *	Retrieve a specific constraint
 *	@param i consraint index (begins with zero)
 *	@return a constraint
 */
adtk_cr3bp_constraint adtk_cr3bp_nodeset::getConstraint(int i) const{
	adtk_cr3bp_constraint temp(constraints.at(i)); 
	return temp;
}

/**
 *	Add a constraint to the nodeset
 *	@param c a constraint to add
 */
void adtk_cr3bp_nodeset::addConstraint(adtk_cr3bp_constraint c){ constraints.push_back(c); }

/**
 *	@return the number of constraints stored in this nodeset
 */
int adtk_cr3bp_nodeset::getNumCons() const { return constraints.size(); }

/**
 *	@return a pointer to the system data object stored in this nodeset
 */
adtk_sys_data* adtk_cr3bp_nodeset::getSysData() { return &sysData; }

void adtk_cr3bp_nodeset::print() {
	printf("CR3BP Nodeset:\n  Nodes:\n");
	for(int n = 0; n < getNumNodes(); n++){
		printf("  > %02d -> [%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f]\n", n,
			nodes[n*nodeSize+0], nodes[n*nodeSize+1], nodes[n*nodeSize+2], 
			nodes[n*nodeSize+3], nodes[n*nodeSize+4], nodes[n*nodeSize+5]);
	}
	for(int c = 0; c < getNumCons(); c++){
		constraints[c].print();
	}
}
/**
 *	@file adtk_cr3bp_nodeset.cpp
 */

#include "adtk_cr3bp_nodeset.hpp"

#include "adtk_cr3bp_sys_data.hpp"
 
#include <cmath>
#include <iostream>

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

adtk_cr3bp_nodeset::adtk_cr3bp_nodeset(const adtk_cr3bp_nodeset& n) : adtk_nodeset(6){
	sysData = n.sysData;
}

adtk_cr3bp_nodeset& adtk_cr3bp_nodeset::operator =(const adtk_cr3bp_nodeset& n){
	sysData = n.sysData;
	return *this;
}

/**
 *	Retrieve a specific constraint
 *	@param i consraint index (begins with zero)
 *	@return a constraint
 */
adtk_cr3bp_constraint adtk_cr3bp_nodeset::getConstraint(int i){ return constraints.at(i); }

/**
 *	Add a constraint to the nodeset
 *	@param c a constraint to add
 */
void adtk_cr3bp_nodeset::addConstraint(adtk_cr3bp_constraint c){ constraints.push_back(c); }
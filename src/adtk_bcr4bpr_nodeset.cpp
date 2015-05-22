/**
 *	@file adtk_bcr4bpr_nodeset.cpp
 */

#include "adtk_bcr4bpr_nodeset.hpp"
#include "adtk_bcr4bpr_sys_data.hpp"

/**
 *	Compute a set of nodes by integrating from initial conditions for some time, then split the
 *	integrated trajectory into pieces (nodes).
 *
 *	@param IC a set of initial conditions, non-dimensional units
 *	@param data a pointer to a system data object that describes the model to integrate in
 *	@param t0 time that corresponds to IC
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC
 *	@param type node distribution type
 */
adtk_bcr4bpr_nodeset::adtk_bcr4bpr_nodeset(double IC[6], adtk_bcr4bpr_sys_data data, double t0,
		double tof, int numNodes, node_distro_t type) : adtk_nodeset(6){

	// initSetFromICs(IC, data, t0, tof, numNodes, type);

	// Compute epoch times for each node
	epochs.reserve(numNodes);

	double ellapsed = 0;
	epochs.push_back(t0);
	for(int n = 0; n < numNodes-1; n++){
		ellapsed += tofs.at(n);
		epochs.push_back(ellapsed);
	}
}//======================================================================

/**
 *	Retrieve a specific constraint
 *	@param i consraint index (begins with zero)
 *	@return a constraint
 */
adtk_bcr4bpr_constraint adtk_bcr4bpr_nodeset::getConstraint(int i){ return constraints.at(i); }

/**
 *	Add a constraint to the nodeset
 *	@param c a constraint to add
 */
void adtk_bcr4bpr_nodeset::addConstraint(adtk_bcr4bpr_constraint c){ constraints.push_back(c); }
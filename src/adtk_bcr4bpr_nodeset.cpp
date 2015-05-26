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
adtk_bcr4bpr_nodeset::adtk_bcr4bpr_nodeset(double IC[6], adtk_bcr4bpr_sys_data data, 
		double t0, double tof, int numNodes, node_distro_t type) : adtk_nodeset(6){

	sysData = data;

	initSetFromICs(IC, &sysData, t0, tof, numNodes, type);

	// Compute epoch times for each node
	epochs.reserve(numNodes);

	double ellapsed = 0;
	epochs.push_back(t0);
	for(int n = 0; n < numNodes-1; n++){
		ellapsed += tofs.at(n);
		epochs.push_back(ellapsed);
	}
}//======================================================================

std::vector<double>* adtk_bcr4bpr_nodeset::getEpochs(){ return &epochs; }

/**
 *	Retrieve a specific constraint
 *	@param i consraint index (begins with zero)
 *	@return a constraint
 */
adtk_bcr4bpr_constraint adtk_bcr4bpr_nodeset::getConstraint(int i) const {
	return constraints.at(i);
}//=====================================

/**
 *	@return the number of constraints stored in this nodeset
 */
int adtk_bcr4bpr_nodeset::getNumCons() const{
	return constraints.size();
}//=====================================

/**
 *	Retrieve a specifi epoch
 *	@param i epoch index (begins with zero)
 *	@return the epoch, non-dimensional units
 */
double adtk_bcr4bpr_nodeset::getEpoch(int i) const {
	return epochs.at(i);
}//=====================================

adtk_sys_data* adtk_bcr4bpr_nodeset::getSysData() { return &sysData; }

/**
 *	Add a constraint to the nodeset
 *	@param c a constraint to add
 */
void adtk_bcr4bpr_nodeset::addConstraint(adtk_bcr4bpr_constraint c){
	constraints.push_back(c);
}//=====================================

/**
 *	Add an epoch to the nodeset
 *	@param d an epoch (non-dimensional time) to add
 */
void adtk_bcr4bpr_nodeset::appendEpoch(double d){
	epochs.push_back(d);
}//=====================================
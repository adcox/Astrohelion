/**
 *	@file adtk_bcr4bpr_nodeset.cpp
 */

#include "adtk_bcr4bpr_nodeset.hpp"
#include "adtk_bcr4bpr_sys_data.hpp"

/**
 *	@brief Compute a set of nodes by integrating from initial conditions for some time, then split the
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

/**
 *	@brief Copy constructor - calls base class copy constructor to handle basic copy
 *	@param n a BCR4BPR nodeset
 */
adtk_bcr4bpr_nodeset::adtk_bcr4bpr_nodeset(const adtk_bcr4bpr_nodeset &n) : adtk_nodeset(n){
	sysData = n.sysData;
	epochs = n.epochs;
}//=========================================

/**
 *	@brief Assignment operator - calls base class assignment operator to handle basic assignment
 *	@param n an input BCR4BPR nodeset
 *	@return this nodeset, made equal to the input nodeset
 */
adtk_bcr4bpr_nodeset& adtk_bcr4bpr_nodeset::operator =(const adtk_bcr4bpr_nodeset &n){
	adtk_nodeset::operator =(n);
	sysData = n.sysData;
	epochs = n.epochs;
	return *this;
}//=========================================

/**
 *	@brief Retrieve a pointer to the vector of epochs
 *	@return a pointer to the beginning of the epochs vector
 */
std::vector<double>* adtk_bcr4bpr_nodeset::getEpochs(){ return &epochs; }

/**
 *	@brief Retrieve a specifi epoch
 *	@param i epoch index (begins with zero)
 *	@return the epoch, non-dimensional units
 */
double adtk_bcr4bpr_nodeset::getEpoch(int i) const {
	return epochs.at(i);
}//=====================================

/**
 *	@brief  Retrieve a pointer to the system data object
 *	@return a pointer to the system data object for this nodeset
 */
adtk_sys_data* adtk_bcr4bpr_nodeset::getSysData() { return &sysData; }

/**
 *	@brief Add an epoch to the nodeset
 *	@param d an epoch (non-dimensional time) to add
 */
void adtk_bcr4bpr_nodeset::appendEpoch(double d){
	epochs.push_back(d);
}//=====================================

/**
 *	@brief Print a textual representation of this object to the standard output
 */
void adtk_bcr4bpr_nodeset::print() const {
	printf("BCR4BPR Nodeset:\n  Nodes:\n");
	for(int n = 0; n < getNumNodes(); n++){
		printf("  > %02d -> [%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f]\n", n,
			nodes[n*nodeSize+0], nodes[n*nodeSize+1], nodes[n*nodeSize+2], 
			nodes[n*nodeSize+3], nodes[n*nodeSize+4], nodes[n*nodeSize+5]);
	}
	for(int c = 0; c < getNumCons(); c++){
		constraints[c].print();
	}
}
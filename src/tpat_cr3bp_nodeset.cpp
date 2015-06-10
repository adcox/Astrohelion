/**
 *	@file tpat_cr3bp_nodeset.cpp
 */

#include "tpat_cr3bp_nodeset.hpp"

#include "tpat_cr3bp_sys_data.hpp"
#include "tpat_utilities.hpp"
 
#include <cmath>
#include <iostream>

/**
 *	@brief Create a general CR3BP nodeset
 */
tpat_cr3bp_nodeset::tpat_cr3bp_nodeset() : tpat_nodeset(6){}

/**
 *	@brief Create a nodeset with specified system data
 *	@param data system data object
 */
tpat_cr3bp_nodeset::tpat_cr3bp_nodeset(tpat_cr3bp_sys_data data) : tpat_nodeset(6){
	sysData = data;
}

/**
 *	@brief Compute a set of nodes by integrating from initial conditions for some time, then split the
 *	integrated trajectory into pieces (nodes).
 *
 *	@param IC a set of initial conditions, non-dimensional units
 *	@param data a pointer to a system data object that describes the model to integrate in
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC
 *	@param type node distribution type
 */
tpat_cr3bp_nodeset::tpat_cr3bp_nodeset(double IC[6], tpat_cr3bp_sys_data data, double tof, int numNodes, 
		node_distro_t type) : tpat_nodeset(6){
	sysData = data;

	initSetFromICs(IC, &sysData, 0, tof, numNodes, type);
}//======================================================================

/**
 *	@brief Copy input nodeset. 
 *
 *	This function calls the base-class copy constructor to
 *	handle copying the generic fields like state and tofs
 *	@param n a nodeset
 */
tpat_cr3bp_nodeset::tpat_cr3bp_nodeset(const tpat_cr3bp_nodeset& n) : tpat_nodeset(n){
	sysData = n.sysData;
}

/**
 *	@brief Make this nodeset equal to the input nodeset
 *	@param n a nodeset
 */
tpat_cr3bp_nodeset& tpat_cr3bp_nodeset::operator =(const tpat_cr3bp_nodeset& n){
	tpat_nodeset::operator =(n);
	sysData = n.sysData;
	return *this;
}//==============================

/**
 *	@return a pointer to the system data object stored in this nodeset
 */
tpat_sys_data* tpat_cr3bp_nodeset::getSysData() { return &sysData; }

/**
 *	@brief Print a textual representation of this object to the standard output
 */
void tpat_cr3bp_nodeset::print() const{
	printf("CR3BP Nodeset:\n  Nodes:\n");
	for(int n = 0; n < getNumNodes(); n++){
		printf("  > %02d -> [%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f]\n", n,
			nodes[n*nodeSize+0], nodes[n*nodeSize+1], nodes[n*nodeSize+2], 
			nodes[n*nodeSize+3], nodes[n*nodeSize+4], nodes[n*nodeSize+5]);
	}
	for(int c = 0; c < getNumCons(); c++){
		constraints[c].print();
	}
}//==================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_cr3bp_nodeset::saveToMat(const char* filename){
	// TODO: Check for propper file extension, add if necessary

	/*	Create a new Matlab MAT file with the given name and optional
	 *	header string. If no header string is given, the default string 
	 *	used containing the software, version, and date in it. If a header
	 *	string is specified, at most the first 116 characters are written to
	 *	the file. Arguments are:
	 *	const char *matname 	- 	the name of the file
	 *	const char *hdr_str 	- 	the 116 byte header string
	 *	enum mat_ft 			- 	matlab file @version MAT_FT_MAT5 or MAT_FT_MAT4
	 */
	mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		printErr("Error creating MAT file\n");
	}else{
		saveNodes(matfp);
		saveTOFs(matfp);
		sysData.saveToMat(matfp);
		// TODO: Add these functions:
		// saveCons(matfp);
		// saveVelCon(matfp);
	}

	Mat_Close(matfp);
}//========================================



/**
 *  @file tpat_nodeset_bcr4bpr.cpp
 *	@brief Derivative of tpat_nodeset, specific to BCR4BPR
 *
 *	@author Andrew Cox
 *	@version September 2, 2015
 *	@copyright GNU GPL v3.0
 */
 
/*
 *  Trajectory Propagation and Analysis Toolkit 
 *  Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
 *
 *  TPAT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  TPAT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with TPAT.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "tpat.hpp"

#include "tpat_nodeset_bcr4bpr.hpp"

#include "tpat_node.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_utilities.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------
/**
 *	@brief Construct a nodeset with no data other than the system
 *	@param data system data object describing the system the nodes exist in
 */
tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(tpat_sys_data_bcr4bpr *data) : tpat_nodeset(data){
	initExtraParam();
}

/**
 *	@brief Compute a set of nodes by integrating from initial conditions for some time, then split the
 *	integrated trajectory into pieces (nodes).
 *
 *	The type is automatically set to splitting the trajectory equally in TIME
 *
 *	@param IC a set of initial conditions, non-dimensional units
 *	@param data a pointer to a system data object that describes the model to integrate in
 *	@param t0 time that corresponds to IC
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC
 */
tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(double IC[6], tpat_sys_data_bcr4bpr *data, 
	double t0, double tof, int numNodes) : tpat_nodeset(data){

	initExtraParam();
	initSetFromICs(IC, data, t0, tof, numNodes, tpat_nodeset::DISTRO_TIME);
	initEpochs(t0);
}//======================================================================

/**
 *	@brief Compute a set of nodes by integrating from initial conditions for some time, then split the
 *	integrated trajectory into pieces (nodes).
 *
 *	The type is automatically set to splitting the trajectory equally in TIME
 *
 *	@param IC a set of initial conditions, non-dimensional units
 *	@param data a pointer to a system data object that describes the model to integrate in
 *	@param t0 time that corresponds to IC
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC
 */
tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(std::vector<double> IC, tpat_sys_data_bcr4bpr *data, 
	double t0, double tof, int numNodes) : tpat_nodeset(data){

	initExtraParam();
	initSetFromICs(&(IC[0]), data, t0, tof, numNodes, tpat_nodeset::DISTRO_TIME);
	initEpochs(t0);
}//======================================================================

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
tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(double IC[6], tpat_sys_data_bcr4bpr *data, 
	double t0, double tof, int numNodes, node_distro_t type) : tpat_nodeset(data){

	initExtraParam();
	initSetFromICs(IC, data, t0, tof, numNodes, type);
	initEpochs(t0);
}//======================================================================

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
tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(std::vector<double> IC, tpat_sys_data_bcr4bpr *data, 
	double t0, double tof, int numNodes, node_distro_t type) : tpat_nodeset(data){

	initExtraParam();
	initSetFromICs(&(IC[0]), data, t0, tof, numNodes, type);
	initEpochs(t0);
}//======================================================================

/**
 *	@brief Create a nodeset as a subset of another
 *	@param orig Original nodeset
 *	@param first index of the first node to be included in the new nodeset
 *	@param last index of the last node to be included in the new nodeset
 */
tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(const tpat_nodeset_bcr4bpr &orig, int first,
	int last) : tpat_nodeset(orig, first, last){}

/**
 *	@brief Copy input nodeset. 
 *
 *	This function calls the base-class copy constructor to
 *	handle copying the generic fields like state and tofs
 *	@param n a nodeset
 */
tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(const tpat_nodeset_bcr4bpr& n) : tpat_nodeset(n) {}

tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(const tpat_arc_data &a) : tpat_nodeset(a) {}

/**
 *	@brief Auto-generate epochs for all nodes
 *
 *	Using the times-of-flight for each node an an initial 
 *	time, compute the epoch for each node assuming time
 *	flows continuously through all nodes
 */
void tpat_nodeset_bcr4bpr::initEpochs(double t0){
	
	// Compute epoch times for each node
	double ellapsed = t0;
	for(size_t n = 0; n < steps.size(); n++){
		steps[n].setExtraParam(1, ellapsed);
		tpat_node *node = static_cast<tpat_node*>(&steps[n]);
		ellapsed += node->getTOF();
	}
}//=========================================


//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

double tpat_nodeset_bcr4bpr::getEpoch(int ix) const {
	if(ix < 0)
		ix += steps.size();

	return steps[ix].getExtraParam(1);
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

void tpat_nodeset_bcr4bpr::initExtraParam(){
	// This function in tpat_nodeset was already called, so 
	// numExtraParam has been set to 1 and a row size has
	// been appended for the TOF variable

	// Add another variable for Epoch Time
	numExtraParam = 2;
	extraParamRowSize.push_back(1);
}//====================================================

void tpat_nodeset_bcr4bpr::saveToMat(const char *filename){
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
		saveState(matfp, "Nodes");
		saveTOFs(matfp);
		saveEpochs(matfp);
		sysData->saveToMat(matfp);
		// TODO: Add these functions:
		// saveCons(matfp);
		// saveVelCon(matfp);
	}

	Mat_Close(matfp);
}//====================================================

/**
 *	@brief Save the epoch values to a file
 *	@param matFile a pointer to the destination matlab file
 */
void tpat_nodeset_bcr4bpr::saveEpochs(mat_t *matFile){
	std::vector<double> epochs;
	for(size_t n = 0; n < steps.size(); n++){
		epochs.push_back(steps[n].getExtraParam(1));
	}
	size_t dims[2] = {epochs.size(), 1};
	matvar_t *matvar = Mat_VarCreate("Epochs", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(epochs[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "Epochs", MAT_COMPRESSION_NONE);
}//====================================================


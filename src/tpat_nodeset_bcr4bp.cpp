/**
 *  @file tpat_nodeset_bcr4bp.cpp
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

#include "tpat_nodeset_bcr4bp.hpp"

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
tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(tpat_sys_data_bcr4bpr *data) : tpat_nodeset(data){
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
tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(double IC[6], tpat_sys_data_bcr4bpr *data, 
	double t0, double tof, int numNodes) : tpat_nodeset(data){

	initExtraParam();
	initSetFromICs(IC, data, t0, tof, numNodes, tpat_nodeset::DISTRO_TIME);
	initEpochs(t0);
}//====================================================

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
tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(std::vector<double> IC, tpat_sys_data_bcr4bpr *data, 
	double t0, double tof, int numNodes) : tpat_nodeset(data){

	initExtraParam();
	initSetFromICs(&(IC[0]), data, t0, tof, numNodes, tpat_nodeset::DISTRO_TIME);
	initEpochs(t0);
}//====================================================

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
tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(double IC[6], tpat_sys_data_bcr4bpr *data, 
	double t0, double tof, int numNodes, node_distro_t type) : tpat_nodeset(data){

	initExtraParam();
	initSetFromICs(IC, data, t0, tof, numNodes, type);
	initEpochs(t0);
}//====================================================

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
tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(std::vector<double> IC, tpat_sys_data_bcr4bpr *data, 
	double t0, double tof, int numNodes, node_distro_t type) : tpat_nodeset(data){

	initExtraParam();
	initSetFromICs(&(IC[0]), data, t0, tof, numNodes, type);
	initEpochs(t0);
}//====================================================

/**
 *	@brief Create a nodeset as a subset of another
 *	@param orig Original nodeset
 *	@param first index of the first node to be included in the new nodeset
 *	@param last index of the last node to be included in the new nodeset
 */
tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(const tpat_nodeset_bcr4bp &orig, int first,
	int last) : tpat_nodeset(orig, first, last){}

/**
 *	@brief Copy input nodeset. 
 *
 *	This function calls the base-class copy constructor to
 *	handle copying the generic fields like state and tofs
 *	@param n a nodeset reference
 */
tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(const tpat_nodeset_bcr4bp& n) : tpat_nodeset(n) {}

/**
 *	@brief Create a BCR4BPR nodeset from its base class
 *	@param a an arc data reference
 */
tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(const tpat_arc_data &a) : tpat_nodeset(a) {}

/**
 *	@brief Auto-generate epochs for all nodes
 *
 *	Using the times-of-flight for each node an an initial 
 *	time, compute the epoch for each node assuming time
 *	flows continuously through all nodes
 *
 *	@param t0 the epoch for the first node
 */
void tpat_nodeset_bcr4bp::initEpochs(double t0){
	
	// Compute epoch times for each node
	double ellapsed = t0;
	for(size_t n = 0; n < steps.size(); n++){
		steps[n].setExtraParam(1, ellapsed);
		tpat_node *node = static_cast<tpat_node*>(&steps[n]);
		ellapsed += node->getTOF();
	}
}//====================================================

/**
 *	@brief Auto-generate epochs for all nodes
 *
 *	Assumes the first node has a properly set initial epoch
 *
 *	Using the times-of-flight for each node an an initial 
 *	time, compute the epoch for each node assuming time
 *	flows continuously through all nodes
 */
void tpat_nodeset_bcr4bp::initEpochs(){
	double ellapsed = 0;
	for(size_t n = 0; n < steps.size(); n++){
		tpat_node *node = static_cast<tpat_node*>(&(steps[n]));

		if(n == 0)
			ellapsed = node->getExtraParam(1);

		node->setExtraParam(1, ellapsed);
		ellapsed += node->getTOF();
	}
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Retrieve the epoch time at a specific node
 *	@param ix node index; if < 0, will count backwards from the end of the nodeset
 *	return the non-dimensional epoch time associated with the specified node
 */
double tpat_nodeset_bcr4bp::getEpoch(int ix) const {
	if(ix < 0)
		ix += steps.size();

	return steps[ix].getExtraParam(1);
}//====================================================

/**
 *	@brief Append a node to the nodeset
 *
 *	Epoch times are automatically generated and updated
 *	for the entire nodeset
 *
 *	@param node a new node to append to the end of the set
 */
void tpat_nodeset_bcr4bp::appendNode(tpat_node node){
	tpat_nodeset::appendNode(node);
	initEpochs();
}//====================================================

/**
 *	@brief Insert a node to the nodeset
 *
 *	Epoch times are automatically generated and updated
 *	for the entire nodeset
 *
 *	@param ix node index; if < 0, will count backwards from end of nodeset
 *	@param node a new node to insert
 */
void tpat_nodeset_bcr4bp::insertNode(int ix, tpat_node node){
	tpat_nodeset::insertNode(ix, node);
	initEpochs();
}//====================================================

/**
 *	@brief Delete a node from the nodeset
 *
 *	Epoch times are automatically generated and updated
 *	for the entire nodeset
 *
 *	@param ix node index; if < 0, will count backwards from end of nodeset
 */
void tpat_nodeset_bcr4bp::deleteNode(int ix){
	tpat_nodeset::deleteNode(ix);
	initEpochs();
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Initialize the extra param vector to hold info specific to this nodeset
 */
void tpat_nodeset_bcr4bp::initExtraParam(){
	// This function in tpat_nodeset was already called, so 
	// numExtraParam has been set to 1 and a row size has
	// been appended for the TOF variable

	// Add another variable for Epoch Time
	numExtraParam = 2;
	extraParamRowSize.push_back(1);
}//====================================================

/**
 *	@brief Save this nodeset to a mat file
 *	@param filename a relative or absolute filepath to a mat file. Note
 *	the extension MUST be ".mat"
 */
void tpat_nodeset_bcr4bp::saveToMat(const char *filename){
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
		saveExtraParam(matfp, 1, "Epochs");
		sysData->saveToMat(matfp);
		// TODO: Add these functions:
		// saveCons(matfp);
		// saveVelCon(matfp);
	}

	Mat_Close(matfp);
}//====================================================

/**
 *	@brief Display a textual representation of this object in the standard output
 */
void tpat_nodeset_bcr4bp::print() const{
	printf("%s Nodeset:\n Nodes: %zu\n", sysData->getTypeStr().c_str(), steps.size());
	for (size_t n = 0; n < steps.size(); n++){
		std::vector<double> node = steps[n].getPosVelState();
		printf("  %02lu: @ %.2f, %13.8f %13.8f %13.8f %13.8f %13.8f %13.8f", n+1, getEpoch(n),
			node.at(0), node.at(1), node.at(2), node.at(3), node.at(4), node.at(5));
		if(n < steps.size()-1){
			printf("   TOF = %.8f\n", getTOF(n));
		}else{
			printf("\n");
		}
	}
	printf(" Constraints:\n");
	for(size_t n = 0; n < steps.size(); n++){
		std::vector<tpat_constraint> nodeCons = getNodeCons(n);
		for(size_t c = 0; c < nodeCons.size(); c++){
			nodeCons[c].print();
		}
	}
}//====================================================


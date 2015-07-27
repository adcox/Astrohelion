/**
 *	@file tpat_nodeset_bcr4bpr.cpp
 *	@brief Derivative of tpat_nodeset, specific to BCR4BPR
 */
/*
 *	Trajectory Propagation and Analysis Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_utilities.hpp"

/**
 *	@brief Construct a nodeset with no data other than the system
 *	@param data system data object describing the system the nodes exist in
 */
tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(tpat_sys_data_bcr4bpr data){
	sysData = data;
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
tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(double IC[6], tpat_sys_data_bcr4bpr data, 
	double t0, double tof, int numNodes){

	sysData = data;

	initSetFromICs(IC, &sysData, t0, tof, numNodes, tpat_nodeset::DISTRO_TIME);
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
tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(std::vector<double> IC, tpat_sys_data_bcr4bpr data, 
	double t0, double tof, int numNodes){

	sysData = data;

	initSetFromICs(&(IC[0]), &sysData, t0, tof, numNodes, tpat_nodeset::DISTRO_TIME);
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
tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(double IC[6], tpat_sys_data_bcr4bpr data, 
	double t0, double tof, int numNodes, node_distro_t type){

	sysData = data;

	initSetFromICs(IC, &sysData, t0, tof, numNodes, type);
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
tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(std::vector<double> IC, tpat_sys_data_bcr4bpr data, 
	double t0, double tof, int numNodes, node_distro_t type){

	sysData = data;

	initSetFromICs(&(IC[0]), &sysData, t0, tof, numNodes, type);
	initEpochs(t0);
}//======================================================================

/**
 *	@brief Create a nodeset as a subset of another
 *	@param orig Original nodeset
 *	@param first index of the first node to be included in the new nodeset
 *	@param last index of the last node to be included in the new nodeset
 */
tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(const tpat_nodeset_bcr4bpr &orig, int first,
	int last) : tpat_nodeset(orig, first, last){

	sysData = orig.sysData;
}//=======================================================

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
	for(size_t n = 0; n < nodes.size(); n++){
		nodes[n].setExtraParam(0, ellapsed);
		ellapsed += nodes[n].getTOF();
	}
}//=========================================

/**
 *	@brief Destructor
 */
tpat_nodeset_bcr4bpr::~tpat_nodeset_bcr4bpr(){}

/**
 *	@brief Copy constructor - calls base class copy constructor to handle basic copy
 *	@param n a BCR4BPR nodeset
 */
tpat_nodeset_bcr4bpr::tpat_nodeset_bcr4bpr(const tpat_nodeset_bcr4bpr &n) : tpat_nodeset(n){
	sysData = n.sysData;
}//=========================================

/**
 *	@brief Assignment operator - calls base class assignment operator to handle basic assignment
 *	@param n an input BCR4BPR nodeset
 *	@return this nodeset, made equal to the input nodeset
 */
tpat_nodeset_bcr4bpr& tpat_nodeset_bcr4bpr::operator =(const tpat_nodeset_bcr4bpr &n){
	tpat_nodeset::operator =(n);
	sysData = n.sysData;
	return *this;
}//=========================================

/**
 *	@brief Concatenate two nodesets
 *
 *	If the final node in <tt>lhs</tt> is the same as the first node in <tt>rhs</tt>, the
 *	concatenation will delete one occurence of the node to achieve continuity. Otherwise
 *	the nodes from <tt>rhs</tt> are concatenated to the end of <tt>lhs</tt>. The velocity
 *	continuity specifications and constraints for <tt>rhs</tt> are updated to reflect the 
 *	new indices of the nodes they describe.
 *
 *	@param lhs
 *	@param rhs
 *	@return a nodeset containing the concatenated input nodesets
 */
tpat_nodeset_bcr4bpr operator +(const tpat_nodeset_bcr4bpr &lhs, const tpat_nodeset_bcr4bpr &rhs){
	if(lhs.sysData != rhs.sysData){
		throw tpat_exception("Cannot add nodesets from different systems; please transform them to be in the same system");
	}

	tpat_nodeset_bcr4bpr temp(lhs.sysData);
	tpat_nodeset_bcr4bpr::basicConcat(lhs, rhs, &temp);

	return temp;
}//=====================================================

/**
 *	@brief Retrieve a specifi epoch
 *	@param i epoch index (begins with zero). If i is negative, the epoch will be
 *	selected from the end of the vector, i.e. -1 will give the last value, -2 will
 *	give the second to last, etc.
 *	@return the epoch, non-dimensional units
 */
double tpat_nodeset_bcr4bpr::getEpoch(int i) const {
	if(i < 0)
		i += nodes.size();
	return nodes[i].getExtraParam(0);
}//=====================================

/**
 *	@brief  Retrieve a pointer to the system data object
 *	@return a pointer to the system data object for this nodeset
 */
tpat_sys_data* tpat_nodeset_bcr4bpr::getSysData() { return &sysData; }

/**
 *	@brief Print a textual representation of this object to the standard output
 */
void tpat_nodeset_bcr4bpr::print() const {
	printf("BCR4BPR Nodeset:\n  Nodes:\n");
	for(size_t n = 0; n < nodes.size(); n++){
		std::vector<double> state = nodes[n].getPosVelState();
		printf("  > %02zu -> @ %9.4f [%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f]", n,
			nodes[n].getExtraParam(0), state[0], state[1], state[2], state[3],
			state[4], state[5]);

		if(n + 1 < nodes.size())
			printf("  TOF = %.4f\n", nodes[n].getTOF());
		else
			printf("\n");
	}

	printf("  Constraints:%s", getNumCons() > 0 ? "\n" : " None\n");
	for(int c = 0; c < getNumCons(); c++){
		constraints[c].print();
	}
}//========================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_nodeset_bcr4bpr::saveToMat(const char* filename){
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
		saveEpochs(matfp);
		sysData.saveToMat(matfp);
		// TODO: Add these functions:
		// saveCons(matfp);
		// saveVelCon(matfp);
	}

	Mat_Close(matfp);
}//========================================

/**
 *	@brief Save the epoch values to a file
 *	@param matFile a pointer to the destination matlab file
 */
void tpat_nodeset_bcr4bpr::saveEpochs(mat_t *matFile){
	std::vector<double> epochs;
	for(size_t n = 0; n < nodes.size(); n++){
		epochs.push_back(nodes[n].getExtraParam(0));
	}
	size_t dims[2] = {epochs.size(), 1};
	matvar_t *matvar = Mat_VarCreate("Epochs", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(epochs[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "Epochs", MAT_COMPRESSION_NONE);
}//=========================================




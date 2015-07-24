/**
 *	@file tpat_nodeset_cr3bp.cpp
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

#include "tpat_nodeset_cr3bp.hpp"

#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_utilities.hpp"
 
#include <cmath>

/**
 *	@brief Create a general CR3BP nodeset
 */
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(){}

/**
 *	@brief Create a nodeset with specified system data
 *	@param data system data object
 */
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(tpat_sys_data_cr3bp data){
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
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(double IC[6], tpat_sys_data_cr3bp data, double tof,
	int numNodes, node_distro_t type){

	sysData = data;

	initSetFromICs(IC, &sysData, 0, tof, numNodes, type);
}//======================================================================

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
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(std::vector<double> IC, tpat_sys_data_cr3bp data, double tof,
	int numNodes, node_distro_t type){

	sysData = data;
	initSetFromICs(&(IC[0]), &sysData, 0, tof, numNodes, type);
}//=====================================================================

/**
 *	@brief Compute a set of nodes by integrating from initial conditions for some time, then split the
 *	integrated trajectory into pieces (nodes).
 *
 *	The type is automatically specified as splitting the trajectory equally in TIME
 *
 *	@param IC a set of initial conditions, non-dimensional units
 *	@param data a pointer to a system data object that describes the model to integrate in
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC
 */
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(double IC[6], tpat_sys_data_cr3bp data, double tof, 
	int numNodes){
	sysData = data;

	initSetFromICs(IC, &sysData, 0, tof, numNodes, tpat_nodeset::DISTRO_TIME);
}//======================================================================

/**
 *	@brief Compute a set of nodes by integrating from initial conditions for some time, then split the
 *	integrated trajectory into pieces (nodes).
 *
 *	The type is automatically specified as splitting the trajectory equally in TIME
 *
 *	@param IC a set of initial conditions, non-dimensional units
 *	@param data a pointer to a system data object that describes the model to integrate in
 *	@param tof duration of the simulation, non-dimensional
 *	@param numNodes number of nodes to create, including IC
 */
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(std::vector<double> IC, tpat_sys_data_cr3bp data, double tof, 
	int numNodes){
	sysData = data;

	initSetFromICs(&(IC[0]), &sysData, 0, tof, numNodes, tpat_nodeset::DISTRO_TIME);
}//======================================================================

/**
 *	@brief Create a noteset by splitting a trajectory into pieces (nodes)
 *
 *	The node distribution type is automatically specified as splitting the trajectory equally in TIME
 *
 *	@param traj the trajectory to split
 *	@param numNodes the number of nodes.
 */
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(tpat_traj_cr3bp traj, int numNodes){
	sysData = traj.getSysData();
	initSetFromTraj(traj, &sysData, numNodes, tpat_nodeset::DISTRO_TIME);
}//===========================================

/**
 *	@brief Create a noteset by splitting a trajectory into pieces (nodes)
 *
 *	@param traj the trajectory to split
 *	@param numNodes the number of nodes.
 *	@param type the node distribution type
 */
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(tpat_traj_cr3bp traj, int numNodes,
	node_distro_t type){
	
	sysData = traj.getSysData();
	initSetFromTraj(traj, &sysData, numNodes, type);
}//===========================================

/**
 *	@brief Create a nodeset as a subset of another
 *	@param orig Original nodeset
 *	@param index of the first node to be included in the new nodeset
 *	@param last index of the last node to be included in the new nodeset
 */
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(const tpat_nodeset_cr3bp &orig, int first,
	int last) : tpat_nodeset(orig, first, last){
	
	sysData = orig.sysData;
}//========================================

/**
 *	@brief Copy input nodeset. 
 *
 *	This function calls the base-class copy constructor to
 *	handle copying the generic fields like state and tofs
 *	@param n a nodeset
 */
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(const tpat_nodeset_cr3bp& n) : tpat_nodeset(n){
	sysData = n.sysData;
}//========================================

/**
 *	@brief Make this nodeset equal to the input nodeset
 *	@param n a nodeset
 */
tpat_nodeset_cr3bp& tpat_nodeset_cr3bp::operator =(const tpat_nodeset_cr3bp& n){
	tpat_nodeset::operator =(n);
	sysData = n.sysData;
	return *this;
}//==============================

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
tpat_nodeset_cr3bp operator +(const tpat_nodeset_cr3bp &lhs, const tpat_nodeset_cr3bp &rhs){
	if(lhs.sysData != rhs.sysData){
		throw tpat_exception("Cannot add nodesets from different systems; please transform them to be in the same system");
	}

	tpat_nodeset_cr3bp temp(lhs.sysData);
	tpat_nodeset_cr3bp::basicConcat(lhs, rhs, &temp);
	return temp;
}//=====================================================

/**
 *	@return a pointer to the system data object stored in this nodeset
 */
tpat_sys_data* tpat_nodeset_cr3bp::getSysData() { return &sysData; }

/**
 *	@brief Print a textual representation of this object to the standard output
 */
void tpat_nodeset_cr3bp::print() const{
	printf("CR3BP Nodeset:\n  Nodes:\n");
	for(size_t n = 0; n < nodes.size(); n++){
		std::vector<double> state = nodes[n].getPosVelState();
		printf("  > %02zu -> [%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f]", n,
			state[0], state[1], state[2], state[3], state[4], state[5]);

		if(n + 1 < nodes.size())
			printf("  TOF = %.4f\n", nodes[n].getTOF());
		else
			printf("\n");
	}

	printf("  Constraints:%s", getNumCons() > 0 ? "\n" : " None\n");
	for(int c = 0; c < getNumCons(); c++){
		constraints[c].print();
	}
}//==================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_nodeset_cr3bp::saveToMat(const char* filename){
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



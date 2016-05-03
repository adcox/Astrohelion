/**
 *  @file tpat_nodeset_cr3bp.cpp
 *	@brief Derivative of tpat_nodeset, specific to CR3BP
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

#include "tpat_nodeset_cr3bp.hpp"

#include "tpat_exceptions.hpp"
#include "tpat_traj_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------
/**
 *	@brief Create a nodeset with specified system data
 *	@param data system data object
 */
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(const tpat_sys_data_cr3bp *data) : tpat_nodeset(data){
	initExtraParam();
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
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(const double IC[6], const tpat_sys_data_cr3bp *data, double tof,
	int numNodes, tpat_nodeDistro_tp type) : tpat_nodeset(data){

	initExtraParam();
	initFromICs(IC, 0, tof, numNodes, type);
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
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(std::vector<double> IC, const tpat_sys_data_cr3bp *data, double tof,
	int numNodes, tpat_nodeDistro_tp type) : tpat_nodeset(data){

	initExtraParam();
	initFromICs(&(IC[0]), 0, tof, numNodes, type);
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
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(const double IC[6], const tpat_sys_data_cr3bp *data, double tof, 
	int numNodes) : tpat_nodeset(data){

	initExtraParam();
	initFromICs(IC, 0, tof, numNodes, tpat_nodeset::DISTRO_TIME);
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
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(std::vector<double> IC, const tpat_sys_data_cr3bp *data, double tof, 
	int numNodes) : tpat_nodeset(data){

	initExtraParam();
	initFromICs(&(IC[0]), 0, tof, numNodes, tpat_nodeset::DISTRO_TIME);
}//======================================================================

/**
 *	@brief Create a noteset by splitting a trajectory into pieces (nodes)
 *
 *	The node distribution type is automatically specified as splitting the trajectory equally in TIME
 *	This function relies on integration to generate nodes, so do not use this if
 *	the trajectory was created by a linearization or other method that uses something
 *	other than the non-linear dynamics to compute points along the trajectory.
 *
 *	@param traj the trajectory to split
 *	@param numNodes the number of nodes.
 */
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(tpat_traj_cr3bp traj, int numNodes) : tpat_nodeset(traj.getSysData()){

	initExtraParam();
	initFromTraj(traj, numNodes, tpat_nodeset::DISTRO_TIME);
}//===========================================

/**
 *	@brief Create a noteset by splitting a trajectory into pieces (nodes)
 *
 *	This function relies on integration to generate nodes, so do not use this if
 *	the trajectory was created by a linearization or other method that uses something
 *	other than the non-linear dynamics to compute points along the trajectory.
 *
 *	@param traj the trajectory to split
 *	@param numNodes the number of nodes.
 *	@param type the node distribution type
 */
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(tpat_traj_cr3bp traj, int numNodes,
	tpat_nodeDistro_tp type) : tpat_nodeset(traj.getSysData()){
	
	initExtraParam();
	initFromTraj(traj, numNodes, type);
}//===========================================

/**
 *	@brief Create a nodeset as a subset of another
 *	@param orig Original nodeset
 *	@param first index of the first node to be included in the new nodeset
 *	@param last index of the last node to be included in the new nodeset
 */
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(const tpat_nodeset_cr3bp &orig, int first,
	int last) : tpat_nodeset(orig, first, last){}

/**
 *	@brief Copy input nodeset. 
 *
 *	This function calls the base-class copy constructor to
 *	handle copying the generic fields like state and tofs
 *	@param n a nodeset
 */
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(const tpat_nodeset_cr3bp& n) : tpat_nodeset(n) {
	initExtraParam();
}//====================================================

/**
 *	@brief Create a CR3BP nodeset from its base class
 *	@param a an arc data reference
 */
tpat_nodeset_cr3bp::tpat_nodeset_cr3bp(const tpat_arc_data &a) : tpat_nodeset(a) {
	initExtraParam();
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *  @brief Get the Jacobi constant value associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @return the Jacobi constant value
 */
double tpat_nodeset_cr3bp::getJacobi(int id) const{
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_nodeset_cr3bp::getJacobi: Node ID out of range");

	return nodes[nodeIDMap[id]].getExtraParam(0);
}//====================================================

/**
 *	@brief Retrieve the value of Jacobi's Constant at the specified step or node
 *	@param ix step index; if < 0, counts backwards from end of nodeset
 *	@return Jacobi at the specified step or node
 */
double tpat_nodeset_cr3bp::getJacobiByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > ((int)nodes.size()))
		throw tpat_exception("tpat_nodeset_cr3bp::getJacobi: invalid index");

	return nodes[ix].getExtraParam(0);
}//====================================================

/**
 *  @brief Set the Jacobi constant value associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @param jacobi Jacobi constant value
 */
void tpat_nodeset_cr3bp::setJacobi(int id, double jacobi){
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw tpat_exception("tpat_nodeset_cr3bp::setJacobi: Node ID out of range");

	nodes[nodeIDMap[id]].setExtraParam(0, jacobi);
}//====================================================

/**
 *	@brief Set Jacobi at the specified step or node
 *	@param ix step index; if < 0, counts backwards from end of nodeset
 *	@param val value of Jacobi
 */
void tpat_nodeset_cr3bp::setJacobiByIx(int ix, double val){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > ((int)nodes.size()))
		throw tpat_exception("tpat_nodeset_cr3bp::setJacobi: invalid index");

	nodes[ix].setExtraParam(0, val);
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Initialize the extra param vector for info specific to this nodeset
 */
void tpat_nodeset_cr3bp::initExtraParam(){
	// Add another variable for Jacobi Constant
	numExtraParam = 1;
	extraParamRowSize.push_back(1);
}//====================================================
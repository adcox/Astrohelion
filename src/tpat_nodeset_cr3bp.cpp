/**
 *  @file tpat_nodeset_cr3bp.cpp
 *	@brief Derivative of TPAT_Nodeset, specific to CR3BP
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
TPAT_Nodeset_CR3BP::TPAT_Nodeset_CR3BP(const TPAT_Sys_Data_CR3BP *data) : TPAT_Nodeset(data){
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
TPAT_Nodeset_CR3BP::TPAT_Nodeset_CR3BP(const double IC[6], const TPAT_Sys_Data_CR3BP *data, double tof,
	int numNodes, tpat_nodeDistro_tp type) : TPAT_Nodeset(data){

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
TPAT_Nodeset_CR3BP::TPAT_Nodeset_CR3BP(std::vector<double> IC, const TPAT_Sys_Data_CR3BP *data, double tof,
	int numNodes, tpat_nodeDistro_tp type) : TPAT_Nodeset(data){

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
TPAT_Nodeset_CR3BP::TPAT_Nodeset_CR3BP(const double IC[6], const TPAT_Sys_Data_CR3BP *data, double tof, 
	int numNodes) : TPAT_Nodeset(data){

	initExtraParam();
	initFromICs(IC, 0, tof, numNodes, TPAT_Nodeset::DISTRO_TIME);
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
TPAT_Nodeset_CR3BP::TPAT_Nodeset_CR3BP(std::vector<double> IC, const TPAT_Sys_Data_CR3BP *data, double tof, 
	int numNodes) : TPAT_Nodeset(data){

	initExtraParam();
	initFromICs(&(IC[0]), 0, tof, numNodes, TPAT_Nodeset::DISTRO_TIME);
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
TPAT_Nodeset_CR3BP::TPAT_Nodeset_CR3BP(TPAT_Traj_CR3BP traj, int numNodes) : TPAT_Nodeset(traj.getSysData()){

	initExtraParam();
	initFromTraj(traj, numNodes, TPAT_Nodeset::DISTRO_TIME);
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
TPAT_Nodeset_CR3BP::TPAT_Nodeset_CR3BP(TPAT_Traj_CR3BP traj, int numNodes,
	tpat_nodeDistro_tp type) : TPAT_Nodeset(traj.getSysData()){
	
	initExtraParam();
	initFromTraj(traj, numNodes, type);
}//===========================================

/**
 *	@brief Create a nodeset as a subset of another
 *	@param orig Original nodeset
 *	@param first index of the first node to be included in the new nodeset
 *	@param last index of the last node to be included in the new nodeset
 */
TPAT_Nodeset_CR3BP::TPAT_Nodeset_CR3BP(const TPAT_Nodeset_CR3BP &orig, int first,
	int last) : TPAT_Nodeset(orig, first, last){}

/**
 *	@brief Copy input nodeset. 
 *
 *	This function calls the base-class copy constructor to
 *	handle copying the generic fields like state and tofs
 *	@param n a nodeset
 */
TPAT_Nodeset_CR3BP::TPAT_Nodeset_CR3BP(const TPAT_Nodeset_CR3BP& n) : TPAT_Nodeset(n) {
	initExtraParam();
}//====================================================

/**
 *	@brief Create a CR3BP nodeset from its base class
 *	@param a an arc data reference
 */
TPAT_Nodeset_CR3BP::TPAT_Nodeset_CR3BP(const TPAT_Base_Arcset &a) : TPAT_Nodeset(a) {
	initExtraParam();
}//====================================================

/**
 *  @brief Create a new nodeset object on the stack
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @param sys pointer to a system data object; should be a 
 *  CR3BP system as the pointer will be cast to that derived class
 *  @return a pointer to the newly created nodeset
 */
baseArcsetPtr TPAT_Nodeset_CR3BP::create( const TPAT_Sys_Data *sys) const{
	const TPAT_Sys_Data_CR3BP *crSys = static_cast<const TPAT_Sys_Data_CR3BP*>(sys);
	return baseArcsetPtr(new TPAT_Nodeset_CR3BP(crSys));
}//====================================================

/**
 *  @brief Create a new nodeset object on the stack that is a 
 *  duplicate of this object
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @return a pointer to the newly cloned nodeset
 */
baseArcsetPtr TPAT_Nodeset_CR3BP::clone() const{
	return baseArcsetPtr(new TPAT_Nodeset_CR3BP(*this));
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
 *  @throw TPAT_Exception if <tt>id</tt> is out of bounds
 */
double TPAT_Nodeset_CR3BP::getJacobi(int id) const{
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw TPAT_Exception("TPAT_Nodeset_CR3BP::getJacobi: Node ID out of range");

	return nodes[nodeIDMap[id]].getExtraParam(0);
}//====================================================

/**
 *	@brief Retrieve the value of Jacobi's Constant at the specified step or node
 *	@param ix step index; if < 0, counts backwards from end of nodeset
 *	@return Jacobi at the specified step or node
 *	@throw TPAT_Exception if <tt>ix</tt> is out of bounds
 */
double TPAT_Nodeset_CR3BP::getJacobiByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > ((int)nodes.size()))
		throw TPAT_Exception("TPAT_Nodeset_CR3BP::getJacobi: invalid index");

	return nodes[ix].getExtraParam(0);
}//====================================================

/**
 *  @brief Set the Jacobi constant value associated with a node
 *  with the specified ID
 * 
 *  @param id the ID of a node
 *  @param jacobi Jacobi constant value
 *  @throw TPAT_Exception if <tt>id</tt> is out of bounds
 */
void TPAT_Nodeset_CR3BP::setJacobi(int id, double jacobi){
	if(id < 0 || id >= (int)(nodeIDMap.size()))
		throw TPAT_Exception("TPAT_Nodeset_CR3BP::setJacobi: Node ID out of range");

	nodes[nodeIDMap[id]].setExtraParam(0, jacobi);
}//====================================================

/**
 *	@brief Set Jacobi at the specified step or node
 *	@param ix step index; if < 0, counts backwards from end of nodeset
 *	@param val value of Jacobi
 *	@throw TPAT_Exception if <tt>ix</tt> is out of bounds
 */
void TPAT_Nodeset_CR3BP::setJacobiByIx(int ix, double val){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > ((int)nodes.size()))
		throw TPAT_Exception("TPAT_Nodeset_CR3BP::setJacobi: invalid index");

	nodes[ix].setExtraParam(0, val);
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Initialize the extra param vector for info specific to this nodeset
 */
void TPAT_Nodeset_CR3BP::initExtraParam(){
	// Add another variable for Jacobi Constant
	numExtraParam = 1;
	extraParamRowSize.push_back(1);
}//====================================================
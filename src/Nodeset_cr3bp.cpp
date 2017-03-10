/**
 *  \file Nodeset_cr3bp.cpp
 *	\brief Derivative of Nodeset, specific to CR3BP
 *
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of Astrohelion
 *
 *  Astrohelion is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Astrohelion is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Astrohelion.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Nodeset_cr3bp.hpp"

#include "Exceptions.hpp"
#include "Traj_cr3bp.hpp"
#include "SysData_cr3bp.hpp"

namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	\brief Create a nodeset with specified system data
 *	\param pData system data object
 */
Nodeset_cr3bp::Nodeset_cr3bp(const SysData_cr3bp *pData) : Nodeset(pData){}

/**
 *	\brief Compute a set of nodes by integrating from initial conditions for some time, then split the
 *	integrated trajectory into pieces (nodes).
 *
 *	\param pData a pointer to a system data object that describes the model to integrate in
 *	\param IC a set of initial conditions, non-dimensional units
 *	\param tof duration of the simulation, non-dimensional
 *	\param numNodes number of nodes to create, including IC
 *	\param type node distribution type (default is NodeDistro_tp::TIME)
 */
Nodeset_cr3bp::Nodeset_cr3bp(const SysData_cr3bp *pData, const double IC[6], double tof,
	int numNodes, NodeDistro_tp type) : Nodeset(pData){

	std::vector<double> ic(IC, IC+6);
	initFromICs(ic, 0, tof, numNodes, type);
}//======================================================================

/**
 *	\brief Compute a set of nodes by integrating from initial conditions for some time, then split the
 *	integrated trajectory into pieces (nodes).
 *
 *	\param pData a pointer to a system data object that describes the model to integrate in
 *	\param IC a set of initial conditions, non-dimensional units
 *	\param tof duration of the simulation, non-dimensional
 *	\param numNodes number of nodes to create, including IC
 *	\param type node distribution type (default is NodeDistro_tp::TIME)
 */
Nodeset_cr3bp::Nodeset_cr3bp(const SysData_cr3bp *pData, std::vector<double> IC, double tof,
	int numNodes, NodeDistro_tp type) : Nodeset(pData){

	initFromICs(IC, 0, tof, numNodes, type);
}//=====================================================================

/**
 *	\brief Create a noteset by splitting a trajectory into pieces (nodes)
 *
 *	This function relies on integration to generate nodes, so do not use this if
 *	the trajectory was created by a linearization or other method that uses something
 *	other than the non-linear dynamics to compute points along the trajectory.
 *
 *	\param traj the trajectory to split
 *	\param numNodes the number of nodes.
 *	\param type the node distribution type
 */
Nodeset_cr3bp::Nodeset_cr3bp(Traj_cr3bp traj, int numNodes,
	NodeDistro_tp type) : Nodeset(traj.getSysData()){
	
	initFromTraj(traj, numNodes, type);
}//===========================================

/**
 *	\brief Create a nodeset as a subset of another
 *	\param orig Original nodeset
 *	\param first index of the first node to be included in the new nodeset
 *	\param last index of the last node to be included in the new nodeset
 */
Nodeset_cr3bp::Nodeset_cr3bp(const Nodeset_cr3bp &orig, int first,
	int last) : Nodeset(orig, first, last){}

/**
 *	\brief Copy input nodeset. 
 *
 *	This function calls the base-class copy constructor to
 *	handle copying the generic fields like state and tofs
 *	\param n a nodeset
 */
Nodeset_cr3bp::Nodeset_cr3bp(const Nodeset_cr3bp& n) : Nodeset(n) {}

/**
 *	\brief Create a CR3BP nodeset from its base class
 *	\param a an arc data reference
 */
Nodeset_cr3bp::Nodeset_cr3bp(const BaseArcset &a) : Nodeset(a) {}

/**
 *  \brief Create a new nodeset object on the stack
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  \param pSys pointer to a system data object; should be a 
 *  CR3BP system as the pointer will be cast to that derived class
 *  @return a pointer to the newly created nodeset
 */
baseArcsetPtr Nodeset_cr3bp::create( const SysData *pSys) const{
	const SysData_cr3bp *crSys = static_cast<const SysData_cr3bp*>(pSys);
	return baseArcsetPtr(new Nodeset_cr3bp(crSys));
}//====================================================

/**
 *  \brief Create a new nodeset object on the stack that is a 
 *  duplicate of this object
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @return a pointer to the newly cloned nodeset
 */
baseArcsetPtr Nodeset_cr3bp::clone() const{
	return baseArcsetPtr(new Nodeset_cr3bp(*this));
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *  \brief Get the Jacobi constant value associated with a node
 *  with the specified ID
 * 
 *  \param id the ID of a node
 *  @return the Jacobi constant value
 *  @throw Exception if <tt>id</tt> is out of bounds
 */
double Nodeset_cr3bp::getJacobi(int id) const{
	if(nodeIDMap.count(id) == 0)
		throw Exception("Nodeset_cr3bp::getJacobi: Node ID out of range");

	return nodes[nodeIDMap.at(id)].getExtraParam("J");
}//====================================================

/**
 *	\brief Retrieve the value of Jacobi's Constant at the specified step or node
 *	\param ix step index; if < 0, counts backwards from end of nodeset
 *	@return Jacobi at the specified step or node
 *	@throw Exception if <tt>ix</tt> is out of bounds
 */
double Nodeset_cr3bp::getJacobiByIx(int ix) const{
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > static_cast<int>(nodes.size()))
		throw Exception("Nodeset_cr3bp::getJacobi: invalid index");

	return nodes[ix].getExtraParam("J");
}//====================================================

/**
 *  \brief Set the Jacobi constant value associated with a node
 *  with the specified ID
 * 
 *  \param id the ID of a node
 *  \param jacobi Jacobi constant value
 *  @throw Exception if <tt>id</tt> is out of bounds
 */
void Nodeset_cr3bp::setJacobi(int id, double jacobi){
	if(nodeIDMap.count(id) == 0)
		throw Exception("Nodeset_cr3bp::setJacobi: Node ID out of range");

	nodes[nodeIDMap.at(id)].setExtraParam("J", jacobi);
}//====================================================

/**
 *	\brief Set Jacobi at the specified step or node
 *	\param ix step index; if < 0, counts backwards from end of nodeset
 *	\param val value of Jacobi
 *	@throw Exception if <tt>ix</tt> is out of bounds
 */
void Nodeset_cr3bp::setJacobiByIx(int ix, double val){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > static_cast<int>(nodes.size()))
		throw Exception("Nodeset_cr3bp::setJacobi: invalid index");

	nodes[ix].setExtraParam("J", val);
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

}// END of Astrohelion namespace
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

#include "tpat_nodeset_bcr4bp.hpp"

#include "tpat_event.hpp"
#include "tpat_node.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_utilities.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------
/**
 *	@brief Construct a nodeset with no data other than the system
 *	@param data system data object describing the system the nodes exist in
 */
tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(const tpat_sys_data_bcr4bpr *data) : tpat_nodeset(data){
	initExtraParam();
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
tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(const double IC[6], const tpat_sys_data_bcr4bpr *data, 
	double t0, double tof, int numNodes) : tpat_nodeset(data){

	initExtraParam();
	initFromICs(IC, t0, tof, numNodes, tpat_nodeset::DISTRO_TIME);
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
tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(std::vector<double> IC, const tpat_sys_data_bcr4bpr *data, 
	double t0, double tof, int numNodes) : tpat_nodeset(data){

	initExtraParam();
	initFromICs(&(IC[0]), t0, tof, numNodes, tpat_nodeset::DISTRO_TIME);
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
tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(const double IC[6], const tpat_sys_data_bcr4bpr *data, 
	double t0, double tof, int numNodes, tpat_nodeDistro_tp type) : tpat_nodeset(data){

	initExtraParam();
	initFromICs(IC, t0, tof, numNodes, type);
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
tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(std::vector<double> IC, const tpat_sys_data_bcr4bpr *data, 
	double t0, double tof, int numNodes, tpat_nodeDistro_tp type) : tpat_nodeset(data){

	initExtraParam();
	initFromICs(&(IC[0]), t0, tof, numNodes, type);
}//====================================================

/**
 *	@brief Create a nodeset as a subset of another
 *	@param orig Original nodeset
 *	@param first index of the first node to be included in the new nodeset
 *	@param last index of the last node to be included in the new nodeset
 */
// tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(const tpat_nodeset_bcr4bp &orig, int first,
// 	int last) : tpat_nodeset(orig, first, last){
// }

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
tpat_nodeset_bcr4bp::tpat_nodeset_bcr4bp(const tpat_base_arcset &a) : tpat_nodeset(a) {}

/**
 *  @brief Create a new nodeset object on the stack
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @param sys pointer to a system data object; should be a 
 *  BCR4BPR system as the pointer will be cast to that derived class
 *  @return a pointer to the newly created nodeset
 */
tpat_nodeset_bcr4bp* tpat_nodeset_bcr4bp::create( const tpat_sys_data *sys) const{
	const tpat_sys_data_bcr4bpr *bcSys = static_cast<const tpat_sys_data_bcr4bpr*>(sys);
	return new tpat_nodeset_bcr4bp(bcSys);
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
tpat_nodeset_bcr4bp* tpat_nodeset_bcr4bp::clone() const{
	return new tpat_nodeset_bcr4bp(*this);
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Initialize the extra param vector to hold info specific to this nodeset
 */
void tpat_nodeset_bcr4bp::initExtraParam(){
	// Nothing to do here!
}//====================================================


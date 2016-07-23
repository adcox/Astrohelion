/**
 *  @file Nodeset_bc4bp.cpp
 *	@brief Derivative of Nodeset, specific to BCR4BPR
 *
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "Nodeset_bc4bp.hpp"

#include "Event.hpp"
#include "Node.hpp"
#include "SimEngine.hpp"
#include "SysData_bc4bp.hpp"
#include "Utilities.hpp"

namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------
/**
 *	@brief Construct a nodeset with no data other than the system
 *	@param data system data object describing the system the nodes exist in
 */
Nodeset_bc4bp::Nodeset_bc4bp(const SysData_bc4bp *data) : Nodeset(data){
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
Nodeset_bc4bp::Nodeset_bc4bp(const double IC[6], const SysData_bc4bp *data, 
	double t0, double tof, int numNodes) : Nodeset(data){

	initExtraParam();
	initFromICs(IC, t0, tof, numNodes, Nodeset::DISTRO_TIME);
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
Nodeset_bc4bp::Nodeset_bc4bp(std::vector<double> IC, const SysData_bc4bp *data, 
	double t0, double tof, int numNodes) : Nodeset(data){

	initExtraParam();
	initFromICs(&(IC[0]), t0, tof, numNodes, Nodeset::DISTRO_TIME);
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
Nodeset_bc4bp::Nodeset_bc4bp(const double IC[6], const SysData_bc4bp *data, 
	double t0, double tof, int numNodes, tpat_nodeDistro_tp type) : Nodeset(data){

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
Nodeset_bc4bp::Nodeset_bc4bp(std::vector<double> IC, const SysData_bc4bp *data, 
	double t0, double tof, int numNodes, tpat_nodeDistro_tp type) : Nodeset(data){

	initExtraParam();
	initFromICs(&(IC[0]), t0, tof, numNodes, type);
}//====================================================

/**
 *	@brief Copy input nodeset. 
 *
 *	This function calls the base-class copy constructor to
 *	handle copying the generic fields like state and tofs
 *	@param n a nodeset reference
 */
Nodeset_bc4bp::Nodeset_bc4bp(const Nodeset_bc4bp& n) : Nodeset(n) {}

/**
 *	@brief Create a BCR4BPR nodeset from its base class
 *	@param a an arc data reference
 */
Nodeset_bc4bp::Nodeset_bc4bp(const BaseArcset &a) : Nodeset(a) {}

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
baseArcsetPtr Nodeset_bc4bp::create( const SysData *sys) const{
	const SysData_bc4bp *bcSys = static_cast<const SysData_bc4bp*>(sys);
	return baseArcsetPtr(new Nodeset_bc4bp(bcSys));
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
baseArcsetPtr Nodeset_bc4bp::clone() const{
	return baseArcsetPtr(new Nodeset_bc4bp(*this));
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
void Nodeset_bc4bp::initExtraParam(){
	// Nothing to do here!
}//====================================================

}// END of Astrohelion namespace
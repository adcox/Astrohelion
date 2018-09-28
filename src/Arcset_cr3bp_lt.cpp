/**
 *  @file Arcset_cr3bp_lt.cpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 1, 2017
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of Astrohelion
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

#include "Arcset_cr3bp_lt.hpp"

#include "Exceptions.hpp"
#include "SysData_cr3bp_lt.hpp"

namespace astrohelion{

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Create an arcset with specified system data
 *	@param pData system data object
 */
Arcset_cr3bp_lt::Arcset_cr3bp_lt(const SysData_cr3bp_lt *pData) : Arcset_cr3bp(pData){}

/**
 *	@brief Copy input arcset. 
 *
 *	This function calls the base-class copy constructor to
 *	handle copying the generic fields like state and tofs
 *	@param n a arcset
 */
Arcset_cr3bp_lt::Arcset_cr3bp_lt(const Arcset_cr3bp_lt& n) : Arcset_cr3bp(n) {}

/**
 *	@brief Create a CR3BP arcset from its base class
 *	@param a an arc data reference
 */
Arcset_cr3bp_lt::Arcset_cr3bp_lt(const BaseArcset &a) : Arcset_cr3bp(a) {}

/**
 *  @brief Create a new arcset object on the stack
 *  @details the `delete` function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @param pSys pointer to a system data object; should be a 
 *  CR3BP system as the pointer will be cast to that derived class
 *  @return a pointer to the newly created arcset
 */
baseArcsetPtr Arcset_cr3bp_lt::create( const SysData *pSys) const{
	const SysData_cr3bp_lt *bcSys = static_cast<const SysData_cr3bp_lt*>(pSys);
	return baseArcsetPtr(new Arcset_cr3bp_lt(bcSys));
}//====================================================

/**
 *  @brief Create a new arcset object on the stack that is a 
 *  duplicate of this object
 *  @details the `delete` function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @return a pointer to the newly cloned arcset
 */
baseArcsetPtr Arcset_cr3bp_lt::clone() const{
	return baseArcsetPtr(new Arcset_cr3bp_lt(*this));
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Retrieve the value of Jacobi's Constant at the specified step
 *	@details If the specified node does not include a precomputed Jacobi
 *	constant value, the value is computed and stored in the node for
 *	future reference.
 *	
 *	@param ix step index; if < 0, counts backwards from end of trajectory
 *	@return Jacobi at the specified step
 *	@throws Exception if `ix` is out of bounds
 */
double Arcset_cr3bp_lt::getHltByIx(int ix){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > static_cast<int>(nodes.size())){
		char msg[64];
		snprintf(msg, 64, "Arcset_cr3bp_lt::getHltByIx: Index %d out of range", ix);
		throw Exception(msg);
	}

	// Find a segment that originates at this node
	int l1 = nodes[ix].getLink(0), l2 = nodes[ix].getLink(1);
	int segID = Linkable::INVALID_ID;
	if(l1 != Linkable::INVALID_ID && 
		segs[segIDMap[l1]].getOrigin() == nodes[ix].getID()){

		segID = l1;
	}else if(l2 != Linkable::INVALID_ID && 
		segs[segIDMap[l2]].getOrigin() == nodes[ix].getID()){

		segID = l2;
	}else{
		// This node is node an origin? Supposed to be impossible
		char msg[64];
		snprintf(msg, 64, "Arcset_cr3bp_lt::getHltByIx: Node %d (ix) is not the "
			"origin of a segment; cannot get H_lt w/o a segment and associated"
			" control law", ix);
		throw Exception(msg);
	}
	
	std::vector<double> state = segs[segIDMap[segID]].getStateByRow(0);

	return DynamicsModel_cr3bp_lt::getHamiltonian(\
		segs[segIDMap[segID]].getTimeByIx(0), &(state[0]), 
		static_cast<const SysData_cr3bp_lt *>(pSysData),
		static_cast<const ControlLaw_cr3bp_lt *>\
		(segs[segIDMap[segID]].getCtrlLaw()));
}//====================================================


//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

}// End of astrohelion namespace
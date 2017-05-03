/**
 *  \file Arcset_cr3bp_lt.cpp
 *	\brief 
 *	
 *	\author Andrew Cox
 *	\version May 1, 2017
 *	\copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "SysData_cr3bp_lt.hpp"

namespace astrohelion{

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	\brief Create an arcset with specified system data
 *	\param pData system data object
 */
Arcset_cr3bp_lt::Arcset_cr3bp_lt(const SysData_cr3bp_lt *pData) : Arcset_cr3bp(pData){}

/**
 *	\brief Copy input arcset. 
 *
 *	This function calls the base-class copy constructor to
 *	handle copying the generic fields like state and tofs
 *	\param n a arcset
 */
Arcset_cr3bp_lt::Arcset_cr3bp_lt(const Arcset_cr3bp_lt& n) : Arcset_cr3bp(n) {}

/**
 *	\brief Create a CR3BP arcset from its base class
 *	\param a an arc data reference
 */
Arcset_cr3bp_lt::Arcset_cr3bp_lt(const BaseArcset &a) : Arcset_cr3bp(a) {}

/**
 *  \brief Create a new arcset object on the stack
 *  \details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  \param pSys pointer to a system data object; should be a 
 *  CR3BP system as the pointer will be cast to that derived class
 *  \return a pointer to the newly created arcset
 */
baseArcsetPtr Arcset_cr3bp_lt::create( const SysData *pSys) const{
	const SysData_cr3bp_lt *bcSys = static_cast<const SysData_cr3bp_lt*>(pSys);
	return baseArcsetPtr(new Arcset_cr3bp_lt(bcSys));
}//====================================================

/**
 *  \brief Create a new arcset object on the stack that is a 
 *  duplicate of this object
 *  \details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  \return a pointer to the newly cloned arcset
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

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

}// End of astrohelion namespace
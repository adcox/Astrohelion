/**
 *  \file Arcset_2bp.cpp
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

#include "Arcset_2bp.hpp"

#include "SysData_2bp.hpp"

namespace astrohelion{

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	\brief Create an arcset with specified system data
 *	\param pSys system data object
 */
Arcset_2bp::Arcset_2bp(const SysData_2bp *pSys) : Arcset(pSys){}

/**
 *	\brief Copy input arcset. 
 *
 *	This function calls the base-class copy constructor to
 *	handle copying the generic fields like state and tofs
 *	\param n a arcset
 */
Arcset_2bp::Arcset_2bp(const Arcset_2bp& n) : Arcset(n) {}

/**
 *	\brief Create a CR3BP arcset from its base class
 *	\param a an arc data reference
 */
Arcset_2bp::Arcset_2bp(const BaseArcset &a) : Arcset(a) {}

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
baseArcsetPtr Arcset_2bp::create( const SysData *pSys) const{
	const SysData_2bp *bcSys = static_cast<const SysData_2bp*>(pSys);
	return baseArcsetPtr(new Arcset_2bp(bcSys));
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
baseArcsetPtr Arcset_2bp::clone() const{
	return baseArcsetPtr(new Arcset_2bp(*this));
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
 *  \brief Execute commands to save data to a Matlab file
 *  \param pMatFile pointer to an open Matlab file
 */
void Arcset_2bp::saveCmds(mat_t* pMatFile, Save_tp saveTp) const{
	Arcset::saveCmds(pMatFile, saveTp);

	saveNodeExtraParam(pMatFile, PARAMKEY_SMA, saveTp, VARNAME_SMA);
	saveNodeExtraParam(pMatFile, PARAMKEY_ECC, saveTp, VARNAME_ECC);
	saveNodeExtraParam(pMatFile, PARAMKEY_ANGMOM, saveTp, VARNAME_ANGMOM);
	saveNodeExtraParam(pMatFile, PARAMKEY_2BP_ENERGY, saveTp, VARNAME_2BP_ENERGY);
}//====================================================

/**
 *  \brief Execute commands to read data from a Matlab file
 *  \param pMatFile pointer to an open Matlab file
 *  \param refLaws Reference to a vector of ControlLaw pointers. As control laws are read
 *  from the Matlab file, unique control laws are constructed and allocated on the stack.
 *  The user must manually delete the ControlLaw objects to avoid memory leaks.
 */
void Arcset_2bp::readCmds(mat_t *pMatFile, std::vector<ControlLaw*> &refLaws){
	Arcset::readCmds(pMatFile, refLaws);

	readNodeExtraParamFromMat(pMatFile, PARAMKEY_SMA, VARNAME_SMA);
	readNodeExtraParamFromMat(pMatFile, PARAMKEY_ECC, VARNAME_ECC);
	readNodeExtraParamFromMat(pMatFile, PARAMKEY_ANGMOM, VARNAME_ANGMOM);
	readNodeExtraParamFromMat(pMatFile, PARAMKEY_2BP_ENERGY, VARNAME_2BP_ENERGY);
}//====================================================
}// End of astrohelion namespace
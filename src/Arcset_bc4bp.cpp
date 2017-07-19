/**
 *  \file Arcset_bc4bp.cpp
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

#include "Arcset_bc4bp.hpp"

#include "Exceptions.hpp"
#include "SysData_bc4bp.hpp"

namespace astrohelion{

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	\brief Create an arcset with specified system data
 *	\param pData system data object
 */
Arcset_bc4bp::Arcset_bc4bp(const SysData_bc4bp *pData) : Arcset(pData){}

/**
 *	\brief Copy input arcset. 
 *
 *	This function calls the base-class copy constructor to
 *	handle copying the generic fields like state and tofs
 *	\param n a arcset
 */
Arcset_bc4bp::Arcset_bc4bp(const Arcset_bc4bp& n) : Arcset(n) {}

/**
 *	\brief Create a CR3BP arcset from its base class
 *	\param a an arc data reference
 */
Arcset_bc4bp::Arcset_bc4bp(const BaseArcset &a) : Arcset(a) {}

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
baseArcsetPtr Arcset_bc4bp::create( const SysData *pSys) const{
	const SysData_bc4bp *bcSys = static_cast<const SysData_bc4bp*>(pSys);
	return baseArcsetPtr(new Arcset_bc4bp(bcSys));
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
baseArcsetPtr Arcset_bc4bp::clone() const{
	return baseArcsetPtr(new Arcset_bc4bp(*this));
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------


//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	\return the angle between the P1/P2 line and the inertial x-axis, radians
 */
double Arcset_bc4bp::getTheta0(){
	const SysData_bc4bp *bcSys = static_cast<const SysData_bc4bp *>(pSysData);
	return bcSys->getTheta0();
}//====================================================

/**
 *	\return the angle between the P2/P3 line (projected into the inertial XY plane)
 *	and the inertial x-axis, radians
 */
double Arcset_bc4bp::getPhi0(){
	const SysData_bc4bp *bcSys = static_cast<const SysData_bc4bp *>(pSysData);
	return bcSys->getPhi0();
}//====================================================

/**
 *	\return the inclination of the P2/P3 orbital plane relative to the P1/P2 orbital
 *	plane, radians
 */
double Arcset_bc4bp::getGamma(){
	const SysData_bc4bp *bcSys = static_cast<const SysData_bc4bp *>(pSysData);
	return bcSys->getGamma();
}//====================================================

/**
 *	\param ix the index of the dqdT vector to retrieve
 *	\return the i'th 6-element dqdT vector. If ix is negative, the count
 *	will proceed from the end of the vector, i.e. -1 will return the final time, 
 *	-2 will give the second to last value, etc.
 *	\throws Exception if <tt>ix</tt> is out of bounds
 */
std::vector<double> Arcset_bc4bp::get_dqdTByIx(int ix){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > static_cast<int>(nodes.size()))
		throw Exception("Arcset_bc4bp::getdqdT: invalid index");

	return getExtraParamVecByIx(ix, PARAMKEY_STATE_EPOCH_DERIV);
}//====================================================

/**
 *	\brief Set the value of the dqdT vector for the specified step
 *	\param ix the index of the step; if < 0, it will count backwards from the end
 *	\param dqdT a pointer to the dqdT vector; this MUST have at least 6 elements,
 *	or the function will read unallocated memory.
 *	\throws Exception if <tt>ix</tt> is out of bounds
 */
void Arcset_bc4bp::set_dqdTByIx(int ix, const double *dqdT){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > static_cast<int>(nodes.size()))
		throw Exception("Arcset_bc4bp::setdqdT: invalid index");

	std::vector<double> dqdT_vec(dqdT, dqdT+6);
	nodes[ix].setExtraParamVec(PARAMKEY_STATE_EPOCH_DERIV, dqdT_vec);
}//====================================================

/**
 *	\brief Set the value of the dqdT vector for the specified step
 *	\param ix the index of the step; if < 0, it will count backwards from the end
 *	\param dqdT a vector (6 elements) representing the dqdT vector
 *	\throws Exception if <tt>ix</tt> is out of bounds
 */
void Arcset_bc4bp::set_dqdTByIx(int ix, std::vector<double> dqdT){
	if(dqdT.size() != 6)
		throw Exception("Arcset_bc4bp::set_dqdT: Cannot accept a dqdT with anything other than 6 elements");

	set_dqdTByIx(ix, &(dqdT[0]));
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *  \brief Execute commands to save data to a Matlab file
 *  \param pMatFile pointer to an open Matlab file
 */
void Arcset_bc4bp::saveCmds(mat_t* pMatFile, Save_tp saveTp) const{
	Arcset::saveCmds(pMatFile, saveTp);

	saveNodeExtraParamVec(pMatFile, PARAMKEY_STATE_EPOCH_DERIV, 6, VARNAME_STATE_EPOCH_DERIV);
}//====================================================

/**
 *  \brief Execute commands to read data from a Matlab file
 *  \param pMatFile pointer to an open Matlab file
 *  \param refLaws Reference to a vector of ControlLaw pointers. As control laws are read
 *  from the Matlab file, unique control laws are constructed and allocated on the stack.
 *  The user must manually delete the ControlLaw objects to avoid memory leaks.
 */
void Arcset_bc4bp::readCmds(mat_t *pMatFile, std::vector<ControlLaw*> &refLaws){
	Arcset::readCmds(pMatFile, refLaws);
	readNodeExtraParamVecFromMat(pMatFile, PARAMKEY_STATE_EPOCH_DERIV, 6, VARNAME_STATE_EPOCH_DERIV);
}//====================================================
}// End of astrohelion namespace
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
#include "Utilities.hpp"

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
 *  \details the `delete` function must be called to 
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
 *  \details the `delete` function must be called to 
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
 *  \param saveTp describes how much data to save
 */
void Arcset_2bp::saveCmds_toFile(mat_t* pMatFile, Save_tp saveTp) const{
	Arcset::saveCmds_toFile(pMatFile, saveTp);

	matvar_t *pSMA = createVar_NodeExtraParam(PARAMKEY_SMA, saveTp, VARNAME_SMA);
	saveVar(pMatFile, pSMA, VARNAME_SMA, MAT_COMPRESSION_NONE);

	matvar_t *pEcc = createVar_NodeExtraParam(PARAMKEY_ECC, saveTp, VARNAME_ECC);
	saveVar(pMatFile, pEcc, VARNAME_ECC, MAT_COMPRESSION_NONE);

	matvar_t *pAngMom = createVar_NodeExtraParam(PARAMKEY_ANGMOM, saveTp, VARNAME_ANGMOM);
	saveVar(pMatFile, pAngMom, VARNAME_ANGMOM, MAT_COMPRESSION_NONE);

	matvar_t *pEnergy = createVar_NodeExtraParam(PARAMKEY_2BP_ENERGY, saveTp, VARNAME_2BP_ENERGY);
	saveVar(pMatFile, pEnergy, VARNAME_2BP_ENERGY, MAT_COMPRESSION_NONE);
}//====================================================

/**
 *  \brief Execute commands to save data to a structure array
 * 
 *  \param pStruct pointer to a structure array
 *  \param ix index of this arcset within the structure array
 *  \param saveTp Describes how much data to save
 */
void Arcset_2bp::saveCmds_toStruct(matvar_t *pStruct, unsigned int ix, Save_tp saveTp) const{
	Arcset::saveCmds_toStruct(pStruct, ix, saveTp);

	if(matvar_t *pSMA = createVar_NodeExtraParam(PARAMKEY_SMA, saveTp, nullptr))
		Mat_VarSetStructFieldByName(pStruct, VARNAME_SMA, ix, pSMA);

	if(matvar_t *pEcc = createVar_NodeExtraParam(PARAMKEY_ECC, saveTp, nullptr))
		Mat_VarSetStructFieldByName(pStruct, VARNAME_ECC, ix, pEcc);

	if(matvar_t *pAngMom = createVar_NodeExtraParam(PARAMKEY_ANGMOM, saveTp, nullptr))
		Mat_VarSetStructFieldByName(pStruct, VARNAME_ANGMOM, ix, pAngMom);

	if(matvar_t *pEnergy = createVar_NodeExtraParam(PARAMKEY_2BP_ENERGY, saveTp, nullptr))
		Mat_VarSetStructFieldByName(pStruct, VARNAME_2BP_ENERGY, ix, pEnergy);
}//====================================================

/**
 *  \brief Execute commands to read data from a Matlab file
 *  \param pMatFile pointer to an open Matlab file
 *  \param refLaws Reference to a vector of ControlLaw pointers. As control laws are read
 *  from the Matlab file, unique control laws are constructed and allocated on the stack.
 *  The user must manually delete the ControlLaw objects to avoid memory leaks.
 */
void Arcset_2bp::readCmds_fromFile(mat_t *pMatFile, std::vector<ControlLaw*> &refLaws){
	Arcset::readCmds_fromFile(pMatFile, refLaws);

	Save_tp saveTp = Save_tp::SAVE_ALL;	// not used in these functions, but provided for forward compatibility
	matvar_t *pSMA = Mat_VarRead(pMatFile, VARNAME_SMA);
	matvar_t *pEcc = Mat_VarRead(pMatFile, VARNAME_ECC);
	matvar_t *pAngMom = Mat_VarRead(pMatFile, VARNAME_ANGMOM);
	matvar_t *pEnergy = Mat_VarRead(pMatFile, VARNAME_2BP_ENERGY);

	if(readVar_NodeExtraParam(pSMA, PARAMKEY_SMA, saveTp)){ Mat_VarFree(pSMA); }
	if(readVar_NodeExtraParam(pEcc, PARAMKEY_ECC, saveTp)){ Mat_VarFree(pEcc); }
	if(readVar_NodeExtraParam(pAngMom, PARAMKEY_ANGMOM, saveTp)){ Mat_VarFree(pAngMom); }
	if(readVar_NodeExtraParam(pEnergy, PARAMKEY_2BP_ENERGY, saveTp)){ Mat_VarFree(pEnergy); }
}//====================================================

/**
 *  \brief Execute commands to read from data from a structure array
 * 
 *  \param pStruct Pointer to the structure array variable
 *  \param ix index of this arcset within the structure array
 * 	\param refLaws Reference to a vector of ControlLaw pointers. As control laws are read
 *  from the Matlab file, unique control laws are constructed and allocated on the stack.
 *  The user must manually delete the ControlLaw objects to avoid memory leaks.
 */
void Arcset_2bp::readCmds_fromStruct(matvar_t *pStruct, unsigned int ix, std::vector<ControlLaw*> &refLaws){
	Arcset::readCmds_fromStruct(pStruct, ix, refLaws);

	Save_tp saveTp = Save_tp::SAVE_ALL;	// not used in these functions, but provided for forward compatibility
	matvar_t *pSMA = Mat_VarGetStructFieldByName(pStruct, VARNAME_SMA, ix);
	matvar_t *pEcc = Mat_VarGetStructFieldByName(pStruct, VARNAME_ECC, ix);
	matvar_t *pAngMom = Mat_VarGetStructFieldByName(pStruct, VARNAME_ANGMOM, ix);
	matvar_t *pEnergy = Mat_VarGetStructFieldByName(pStruct, VARNAME_2BP_ENERGY, ix);

	if(readVar_NodeExtraParam(pSMA, PARAMKEY_SMA, saveTp)){ Mat_VarFree(pSMA); }
	if(readVar_NodeExtraParam(pEcc, PARAMKEY_ECC, saveTp)){ Mat_VarFree(pEcc); }
	if(readVar_NodeExtraParam(pAngMom, PARAMKEY_ANGMOM, saveTp)){ Mat_VarFree(pAngMom); }
	if(readVar_NodeExtraParam(pEnergy, PARAMKEY_2BP_ENERGY, saveTp)){ Mat_VarFree(pEnergy); }
}//====================================================

}// End of astrohelion namespace
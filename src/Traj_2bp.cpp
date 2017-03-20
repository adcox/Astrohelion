/**
 *  \file Traj_2bp.cpp
 *	\brief Derivative of Traj, specific to 2BP
 *
 *	\author Andrew Cox
 *	\version August 25, 2016
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

#include "Traj_2bp.hpp"

#include "Exceptions.hpp"
#include "Node.hpp"
#include "SysData_2bp.hpp"
#include "Utilities.hpp"


namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------
 
/**
 *	\brief Create a trajectory for a specific system
 *	\param sys a pointer to a system data object
 */
Traj_2bp::Traj_2bp(const SysData_2bp *sys) : Traj(sys){}

/**
 *	\brief Create a trajectory from another trajectory
 *	\param t a trajectory reference
 */
Traj_2bp::Traj_2bp(const Traj_2bp &t) : Traj(t){}

/**
 *	\brief Create a trajectory from its base class
 *	\param a an arc data reference
 */
Traj_2bp::Traj_2bp(const BaseArcset &a) : Traj(a){}

/**
 *  \brief Create a new trajectory object on the stack
 *  \details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  \param sys pointer to a system data object; should be a 
 *  CR3BP system as the pointer will be cast to that derived class
 *  \return a pointer to the newly created trajectory
 */
baseArcsetPtr Traj_2bp::create( const SysData *sys) const{
	const SysData_2bp *crSys = static_cast<const SysData_2bp*>(sys);
	return baseArcsetPtr(new Traj_2bp(crSys));
}//====================================================

/**
 *  \brief Create a new trajectory object on the stack that is a 
 *  duplicate of this object
 *  \details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 *  
 *  \return a pointer to the newly cloned trajectory
 */
baseArcsetPtr Traj_2bp::clone() const{
	return baseArcsetPtr(new Traj_2bp(*this));
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *	\brief Concatenate two trajectory objects
 *
 * 	When adding A + B, if the final state of A and initial state
 *	of B are the same, this algorithm will skip the initial state
 *	of B in the concatenation to avoid duplicating a state. This
 *	method also overrides the base class behavior and forces time to be
 *	continuous along the concatentated trajectory regardless of whether
 *	the final state of A and in itial state of B are the same
 *
 *	\param rhs the right-hand-side of the addition operation
 *	\return a reference to the concatenated arcset object
 */
Traj& Traj_2bp::operator +=(const Traj &rhs){
	// Create a copy of rhs (it is const)
	Traj temp(rhs);

	// Shift the time in temp by the final time in this trajectory
	double tf = getTimeByIx(-1);
	for(int s = 0; s < temp.getNumNodes(); s++){
		double t = tf + temp.getTimeByIx(s);
		temp.setTimeByIx(s, t);
	}

	// throw Exception("Traj_2bp::operator +=: Not currently implemented!");
	Traj::operator +=(temp);

	return *this;
}//====================================================

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
void Traj_2bp::saveCmds(mat_t* pMatFile) const{
	Traj::saveCmds(pMatFile);

	saveExtraParam(pMatFile, PARAMKEY_SMA, VARNAME_SMA);
	saveExtraParam(pMatFile, PARAMKEY_ECC, VARNAME_ECC);
	saveExtraParam(pMatFile, PARAMKEY_ANGMOM, VARNAME_ANGMOM);
	saveExtraParam(pMatFile, PARAMKEY_2BP_ENERGY, VARNAME_2BP_ENERGY);
}//====================================================

/**
 *  \brief Execute commands to read data from a Matlab file
 *  \param pMatFile pointer to an open Matlab file
 */
void Traj_2bp::readCmds(mat_t *pMatFile){
	Traj::readCmds(pMatFile);

	readExtraParamFromMat(pMatFile, PARAMKEY_SMA, VARNAME_SMA);
	readExtraParamFromMat(pMatFile, PARAMKEY_ECC, VARNAME_ECC);
	readExtraParamFromMat(pMatFile, PARAMKEY_ANGMOM, VARNAME_ANGMOM);
	readExtraParamFromMat(pMatFile, PARAMKEY_2BP_ENERGY, VARNAME_2BP_ENERGY);
}//====================================================


}// END of Astrohelion namespace
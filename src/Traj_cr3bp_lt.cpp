/**
 *  @file Traj_cr3bp_lt.cpp
 *	@brief 
 *
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
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

#include "Traj_cr3bp_lt.hpp"


#include "SysData_cr3bp_lt.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"
 

namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Create a trajectory for a specific system
 *	@param sys a pointer to a system data object
 */
Traj_cr3bp_lt::Traj_cr3bp_lt(const SysData_cr3bp_lt* sys) : Traj_cr3bp(sys){}

/**
 *	@brief Create a trajectory from another trajectory
 *	@param t a trajectory reference
 */
Traj_cr3bp_lt::Traj_cr3bp_lt(const Traj_cr3bp_lt &t) : Traj_cr3bp(t){}

/**
 *	@brief Create a trajectory from its base class
 *	@param a an arc data reference
 */
Traj_cr3bp_lt::Traj_cr3bp_lt(const BaseArcset &a) : Traj_cr3bp(a){}

/**
 *  @brief Create a new trajectory object on the stack
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @param sys pointer to a system data object; should be a 
 *  CR3BP LTVP system as the pointer will be cast to that derived class
 *  @return a pointer to the newly created trajectory
 */
baseArcsetPtr Traj_cr3bp_lt::create( const SysData *sys) const{
	const SysData_cr3bp_lt *crSys = static_cast<const SysData_cr3bp_lt*>(sys);
	return baseArcsetPtr(new Traj_cr3bp_lt(crSys));
}//====================================================

/**
 *  @brief Create a new trajectory object on the stack that is a 
 *  duplicate of this object
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @return a pointer to the newly cloned trajectory
 */
baseArcsetPtr Traj_cr3bp_lt::clone() const{
	return baseArcsetPtr(new Traj_cr3bp_lt(*this));
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

}
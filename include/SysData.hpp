/**
 *  @file SysData.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2016, Andrew Cox; Protected under the GNU GPL v3.0
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
#pragma once

#include "Core.hpp"
 
#include "matio.h"
#include <string>
#include <vector>


namespace astrohelion{

// Forward declarations
class DynamicsModel;

/**
 *	@brief Specifies the type of dynamical system
 *	
 *	Describes the type of system described by this data object. This
 *	is useful when a function or object specifies only the super-class
 *	SysData and the user needs to determine which derivative
 *	data type to use. Each derivative type will set its type variable to
 *	one of the type values.
 */
enum class SysData_tp {
	UNDEF_SYS,		//!< System type is undefined
	R2BP_SYS,		//!< Relative Motion 2-Body Problem
	CR3BP_SYS,		//!< Circular Restricted 3-Body Problem
	CR3BP_LTVP_SYS,	//!< Circular Restricted 3-Body Problem with Low Thrust, Velocity-Pointing;
	BCR4BPR_SYS 	//!< Bi-Circular Restricted 4-Body Problem, rotating frame
};

/**
 *	@ingroup model
 *	@brief Contains information about a system, like mass ratio, primary names, etc.
 *
 *	This object is a super-class for system data objects, "system" being
 *	a dynamical system or model. This class is an abstract base class
 *	with pure virtual functions and cannot be instantiated as a standalone
 *	object. Other derivative classes can be instantiated as objects and have
 *	variables and functions specific to different dynamical systems.
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 *	@copyright GNU GPL v3.0
 */
class SysData : public Core{

	public:
		SysData();
		SysData(const SysData&);
		virtual ~SysData();
		

		SysData& operator =(const SysData&);

		friend bool operator ==(const SysData&, const SysData&);
		friend bool operator !=(const SysData&, const SysData&);

		/**
		 *  @name Set and Get Functions
		 *  @{
		 */
		double getCharL() const;
		double getCharM() const;
		double getCharT() const;

		/**
		 *	@brief Retrieve the model that governs the motion for this system type
		 *	@return the model that governs the motion for this system type
		 */
		virtual const DynamicsModel* getDynamicsModel() const = 0;
		
		int getNumPrimaries() const;
		std::string getPrimary(int n) const;
		int getPrimID(int n) const;
		SysData_tp getType() const;
		std::string getTypeStr() const;
		//@}
		
		virtual void saveToMat(const char*) const;

		/**
		 *	@brief Save the system data object to a file
		 *	@param matFile a pointer to an open mat file
		 */
		virtual void saveToMat(mat_t *matFile) const = 0;
		
	protected:
		/** Number of primaries that exist in this system */
		int numPrimaries = 0;

		/** Vector containing names of primaries */
		std::vector<std::string> primaries {};

		/** Vector containing IDs for primaries */
		std::vector<int> primIDs {};

		/** Characteristic length (km) */
		double charL = 0;

		/** Characteristic time (sec) */
		double charT = 0;

		/** Characteristic mass (kg) */
		double charM = 0;

		/** Vector containing other parameters, like mu, k, nu, etc. */
		std::vector<double> otherParams {};

		/** The type of system this data object describes */
		SysData_tp type = SysData_tp::UNDEF_SYS;

		void copyData(const SysData&);

		/**
		 *  @brief Read the system data object from a Matlab data file
		 *  @param matfile A pointer to the open matfile
		 */
		virtual void readFromMat(mat_t *matfile) = 0;
};

}// END of Astrohelion namespace
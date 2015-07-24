/*
 *	Trajectory Propagation and Analysis Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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
#ifndef H_SYSDATA
#define H_SYSDATA

#include <string>
#include <vector>

//forward declarations
class tpat_model;

/**
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
class tpat_sys_data{

	public:
		/**
		 *	@brief Specifies the type of dynamical system
		 *	
		 *	Describes the type of system described by this data object. This
		 *	is useful when a function or object specifies only the super-class
		 *	tpat_sys_data and the user needs to determine which derivative
		 *	data type to use. Each derivative type will set its type variable to
		 *	one of the type values.
		 */
		enum system_t {
			UNDEF_SYS,		//!< System type is undefined
			CR3BP_SYS,		//!< Circular Restricted 3-Body Problem
			CR3BP_LTVP_SYS,	//!< Circular Restricted 3-Body Problem with Low Thrust, Velocity-Pointing;
			BCR4BPR_SYS 	//!< Bi-Circular Restricted 4-Body Problem, rotating frame
		};

		tpat_sys_data();
		tpat_sys_data(const tpat_sys_data&);
		virtual ~tpat_sys_data();
		

		tpat_sys_data& operator =(const tpat_sys_data&);

		friend bool operator ==(const tpat_sys_data&, const tpat_sys_data&);
		friend bool operator !=(const tpat_sys_data&, const tpat_sys_data&);

		double getCharL() const;
		double getCharM() const;
		double getCharT() const;
		virtual tpat_model* getModel() = 0;
		int getNumPrimaries() const;
		std::string getPrimary(int n) const;
		int getPrimID(int n) const;
		int getStateWidth() const;
		system_t getType() const;
		std::string getTypeStr() const;
		
	protected:
		/** Number of primaries that exist in this system */
		int numPrimaries = 0;

		/** Vector containing names of primaries */
		std::vector<std::string> primaries;

		/** Vector containing IDs for primaries */
		std::vector<int> primIDs;

		/** Characteristic length (km) */
		double charL = 0;

		/** Characteristic time (sec) */
		double charT = 0;

		/** Characteristic mass (kg) */
		double charM = 0;

		/** Vector containing other parameters, like mu, k, nu, etc. */
		std::vector<double> otherParams;

		/** The type of system this data object describes */
		system_t type = UNDEF_SYS;

		void copyData(const tpat_sys_data&);
};
#endif
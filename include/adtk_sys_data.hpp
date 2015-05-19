/**
 *	This object is a super-class for system data objects, "system" being
 *	a dynamical system or model. This class is an abstract base class
 *	with pure virtual functions and cannot be instantiated as a standalone
 *	object. Other derivative classes can be instantiated as objects and have
 *	variables and functions specific to different dynamical systems.
 *
 *	Author: Andrew Cox
 *
 *	Version: May 15, 2015
 */

/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (ADTK).
 *
 *  ADTK is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ADTK is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ATDK.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __H_SYSDATA_
#define __H_SYSDATA_

#include <string>

class adtk_sys_data{

	public:
		/**
		 *	System Type - Enumerated Type
		 *	
		 *	Describes the type of system described by this data object. This
		 *	is useful when a function or object specifies only the super-class
		 *	adtk_sys_data and the user needs to determine which derivative
		 *	data type to use. Each derivative type will set its type variable to
		 *	one of the type values.
		 *
		 *	Values:
		 *		UNDEF_SYS 	-	System type is undefined; no derivative type defined.
		 *		CR3BP_SYS 	- 	Circular Restricted 3-Body Problem; applies to adtk_cr3bp_sys_data
		 *		BCR4BPR_SYS - 	Bi-Circular Restricted 4-Body Problem, rotating frame; applies
		 *						to adtk_bcr4bpr_sys
		 */
		enum system_t {UNDEF_SYS, CR3BP_SYS, BCR4BPR_SYS};

		adtk_sys_data();
		virtual ~adtk_sys_data() {}

		adtk_sys_data& operator= (const adtk_sys_data&);

		double getCharL();
		double getCharT();
		double getCharM();
		system_t getType();
		std::string getTypeStr();

		/**
		 *	Each derivative class must define this function to provide the name of each primary
		 *	@param n the "index" of the primary, starts at 0
		 *	@return the name of the n'th primary
		 */
		virtual std::string getPrimary(int n) = 0;

	protected:

		/** Number of primaries that exist in this system */
		int numPrimaries;

		/** Characteristic length (km) */
		double charL;

		/** Characteristic time (sec) */
		double charT;

		/** Characteristic mass (kg) */
		double charM;

		/** The type of system this data object describes */
		system_t type;
};
#endif
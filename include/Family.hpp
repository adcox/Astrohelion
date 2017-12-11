/**
 * \file Family.hpp
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrohelion.
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

#include <string>
#include <vector>

#include "Core.hpp"

#include "SysData.hpp"

namespace astrohelion{

// Forward Declarations
class ControlLaw;

/**
 *	\ingroup fam
 *	\brief How to sort members of this family
 *	
 *	Families typically evolve naturally through one variable, and
 *	it makes sense to analyze trends using this parameter as the 
 *	independent variable.
 *
 *	The integer value of each type is set up to map to the index
 *	of the independent variable in the continuation process
 */
enum class FamSort_tp : int{
	SORT_X 		= 0,	//!< Sort by the x-coordinate in the IC
	SORT_Y 		= 1,	//!< Sort by the y-coordinate in the IC
	SORT_Z 		= 2,	//!< Sort by the z-coordinate in the IC
	SORT_VX 	= 3,	//!< Sort by the x velocity component in the IC
	SORT_VY  	= 4,	//!< Sort by the y velocity component in the IC
	SORT_VZ 	= 5,	//!< Sort by the z velocity component in the IC
	SORT_TOF 	= 6,	//!< Sort by Time-Of-Flight
	SORT_JC 	= 7,	//!< Sort by Jacobi Constant
	SORT_NONE 	= 999	//!< Do not adjust sorting; sortMembers() will do nothing
};

/**
 *  \ingroup fam
 *  \brief [brief description]
 *  \details [long description]
 *  \return [description]
 */
class Family : public Core{
	public:
		/**
		 *  \name *structors
		 *  \{
		 */
		Family(const SysData*);
		Family(const Family&);
		virtual ~Family();
		//\}

		/**
		 *  \name Operators
		 *  \{
		 */
		Family& operator= (const Family&);
		//\}

		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		std::string getName() const;
		double getMatchTol() const;
		virtual unsigned int getNumMembers() const = 0;
		FamSort_tp getSortType() const;
		const char* getSortTypeStr() const;
		const SysData* getSysData() const;

		void setMatchTol(double);
		void setName(std::string);
		void setSortType(FamSort_tp);
		//\}
		
		/**
		 *  \name Utility Functions
		 *  \{
		 */
		static const char* sortTypeToStr(FamSort_tp);
		//\}
		
		/**
		 *  \name File I/O
		 *  \{
		 */
		virtual void readFromMat(const char*, std::vector<ControlLaw*>&, bool bReconstruct = false) = 0;
		virtual void saveToMat(const char*) const = 0;
		//\}

		/**
		 *  \name Analysis Functions
		 *  \{
		 */
		virtual void sortMembers() = 0;
		virtual void reverseOrder() = 0;
		//\}

	protected:
		/** A pointer to the system data object that the describes the system this family exists in */
		const SysData *pSysData;

		std::string name = "UNNAMED";						//!< Descriptive name of family
		FamSort_tp sortType = FamSort_tp::SORT_NONE;		//!< Describes the variable to sort family members by

		double matchTol = 1e-8;								//!< Acceptable tolerance (non-dim units) when locating a member by a specific attribute

		const char* VARNAME_SORTTYPE = "SortType";		//!< Variable name for the sort type
		const char* VARNAME_NAME = "Name";					//!< Variable name for the name of the family

		/**
		 *  \name Utility Functions
		 *  \{
		 */
		virtual void copyMe(const Family&);
		//\}

		/**
		 *  \name Analysis Functions
		 *  \{
		 */
		std::vector<unsigned int> findMatches(double, std::vector<double>*) const;
		//\}
};

}// End of astrohelion namespace
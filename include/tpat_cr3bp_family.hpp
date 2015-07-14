/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (TPAT).
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

#ifndef __H_FAMILIY_
#define __H_FAMILIY_

#include "tpat_cr3bp_family_member.hpp"
#include "tpat_cr3bp_sys_data.hpp"

#include <cmath>
#include <string>
#include <vector>

// Forward declarations
class tpat_constraint;

/**
 *	@brief A data structure to represent a family of orbits
 *
 *	Currently, this structure assumes family members are periodic.
 */
class tpat_cr3bp_family{
	public:
		/**
		 *	@brief How to sort members of this family
		 *	
		 *	Families typically evolve naturally through one variable, and
		 *	it makes sense to analyze trends using this parameter as the 
		 *	independent variable.
		 */
		enum sortVar_t{
			SORT_X,		//!< Sort by the x-coordinate in the IC
			SORT_Y,		//!< Sort by the y-coordinate in the IC
			SORT_Z,		//!< Sort by the z-coordinate in the IC
			SORT_VX,	//!< Sort by the x velocity component in the IC
			SORT_VY,	//!< Sort by the y velocity component in the IC
			SORT_VZ,	//!< Sort by the z velocity component in the IC
			SORT_JC,	//!< Sort by Jacobi Constant
			SORT_TOF	//!< Sort by Time-Of-Flight
		};

		tpat_cr3bp_family(tpat_cr3bp_sys_data);
		tpat_cr3bp_family(std::string);
		tpat_cr3bp_family(const tpat_cr3bp_family&);
		~tpat_cr3bp_family();

		tpat_cr3bp_family& operator= (const tpat_cr3bp_family&);

		void addMember(tpat_cr3bp_family_member);
		
		tpat_cr3bp_family_member getMember(int) const;
		std::vector<tpat_cr3bp_family_member> getMemberByStateVar(double, int) const;
		std::vector<tpat_cr3bp_family_member> getMemberByTOF(double) const;
		std::vector<tpat_cr3bp_family_member> getMemberByJacobi(double) const;
		std::string getName() const;
		sortVar_t getSortType() const;
		tpat_cr3bp_sys_data getSysData() const;

		void setName(std::string);
		void setSortType(sortVar_t);

		void saveToFile(std::string filepath) const;
		void sort();

	protected:
		std::string name = "NULL";						//!< Descriptive name of the family
		std::vector<tpat_cr3bp_family_member> members;	//!< Contains all family members
		tpat_cr3bp_sys_data sysData;					//!< Describes the system this family exists in
		sortVar_t sortType;								//!< Describes the most natural variable to sort family members by

		double matchTol = 1e-8;		//!< Acceptable tolerance (non-dim units) when locating a member by a specific attribute
		double numNodes = 4;		//!< Number of nodes to use when representing a family member

		void copyMe(const tpat_cr3bp_family&);
		std::vector<int> findMatches(double, std::vector<double>*) const;
		std::vector<tpat_cr3bp_family_member> getMatchingMember(double, std::vector<double>*, tpat_constraint) const;
		void getCoord(int, std::vector<double>*) const;
};

#endif
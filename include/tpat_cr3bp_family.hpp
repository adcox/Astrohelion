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

#ifndef H_FAMILIY
#define H_FAMILIY

#include "tpat_cr3bp_family_member.hpp"
#include "tpat_cr3bp_sys_data.hpp"

#include "matio.h"

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
		 *
		 *	The integer value of each type is set up to map to the index
		 *	of the independent variable in the continuation process
		 */
		enum sortVar_t{
			SORT_X 		= 0,	//!< Sort by the x-coordinate in the IC
			SORT_Y 		= 1,	//!< Sort by the y-coordinate in the IC
			SORT_Z 		= 2,	//!< Sort by the z-coordinate in the IC
			SORT_VX 	= 3,	//!< Sort by the x velocity component in the IC
			SORT_VY  	= 4,	//!< Sort by the y velocity component in the IC
			SORT_VZ 	= 5,	//!< Sort by the z velocity component in the IC
			SORT_JC 	= 7,	//!< Sort by Jacobi Constant
			SORT_TOF 	= 6		//!< Sort by Time-Of-Flight
		};

		tpat_cr3bp_family(tpat_cr3bp_sys_data);
		tpat_cr3bp_family(const char*);
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
		const char* getSortTypeStr() const;
		tpat_cr3bp_sys_data getSysData() const;

		void setName(std::string);
		void setSortType(sortVar_t);

		void saveToMat(const char*);
		void sortEigs();
		void sortMembers();


	protected:
		std::string name = "NULL";						//!< Descriptive name of the family
		std::vector<tpat_cr3bp_family_member> members;	//!< Contains all family members
		tpat_cr3bp_sys_data sysData;					//!< Describes the system this family exists in
		sortVar_t sortType = SORT_X;								//!< Describes the most natural variable to sort family members by

		double matchTol = 1e-8;		//!< Acceptable tolerance (non-dim units) when locating a member by a specific attribute
		double numNodes = 4;		//!< Number of nodes to use when representing a family member

		const char* DATA_VAR_NAME = "MemberData";
		const char* SORTTYPE_VAR_NAME = "SortType";
		const char* NAME_VAR_NAME = "Name";
		const char* EIG_VAR_NAME = "Eigenvalues";
		const size_t DATA_WIDTH = 11;

		void copyMe(const tpat_cr3bp_family&);
		std::vector<int> findMatches(double, std::vector<double>*) const;
		std::vector<tpat_cr3bp_family_member> getMatchingMember(double, std::vector<double>*, tpat_constraint) const;
		void getCoord(int, std::vector<double>*) const;

		void loadMemberData(mat_t*);
		void loadEigVals(mat_t*);
		void saveMembers(mat_t*);
		void saveMiscData(mat_t*);
		void saveEigVals(mat_t*);
};

#endif
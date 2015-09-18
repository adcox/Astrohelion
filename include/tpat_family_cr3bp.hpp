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

#ifndef H_FAMILIY
#define H_FAMILIY

#include "tpat_family_member_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"

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
class tpat_family_cr3bp{
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

		/**
		 *	@brief Eigenvalue pair types
		 *
		 *	Eigenvalues come in three different types of pairs. They will
		 *	either be complex, real, or exactly equal to 1.0. Since the 
		 *	matrices are real, all complex eigenvalues will come in conjucate
		 *	pairs. In this type of problem (trajectory design), the other STM
		 *	eigenvalues also come in pairs: 2 1.0 eigenvalues, or two real 
		 *	eigenvalues that are reciprocals.
		 */
		enum eigValSet_t{
			EIGSET_COMP_CONJ,	//!< Complex conjugate pair
			EIGSET_ONES,		//!< Exactly equal to 1.0
			EIGSET_REAL_RECIP	//!< Real, reciprocal pair
		};

		// *structors
		tpat_family_cr3bp(tpat_sys_data_cr3bp);
		tpat_family_cr3bp(const char*);
		tpat_family_cr3bp(const tpat_family_cr3bp&);
		~tpat_family_cr3bp();

		// Operators
		tpat_family_cr3bp& operator= (const tpat_family_cr3bp&);

		// Set and Get Methods
		void addMember(tpat_family_member_cr3bp);
		
		tpat_family_member_cr3bp getMember(int) const;
		std::vector<tpat_family_member_cr3bp> getMemberByStateVar(double, int) const;
		std::vector<tpat_family_member_cr3bp> getMemberByTOF(double) const;
		std::vector<tpat_family_member_cr3bp> getMemberByJacobi(double) const;
		std::string getName() const;
		sortVar_t getSortType() const;
		const char* getSortTypeStr() const;
		tpat_sys_data_cr3bp getSysData() const;
		tpat_sys_data_cr3bp* getSysDataPtr();

		void setName(std::string);
		void setSortType(sortVar_t);

		// Utility and Operational Functions
		std::vector<int> findBifurcations();
		void saveToMat(const char*);
		void sortEigs();
		void sortMembers();


	protected:
		std::string name = "NULL";						//!< Descriptive name of the family
		std::vector<tpat_family_member_cr3bp> members;	//!< Contains all family members
		tpat_sys_data_cr3bp sysData;					//!< Describes the system this family exists in
		sortVar_t sortType = SORT_X;					//!< Describes the most natural variable to sort family members by

		double matchTol = 1e-8;		//!< Acceptable tolerance (non-dim units) when locating a member by a specific attribute
		double numNodes = 4;		//!< Number of nodes to use when representing a family member

		const char* DATA_VAR_NAME = "MemberData";		//!< Variable name for family member data
		const char* SORTTYPE_VAR_NAME = "SortType";		//!< Variable name for the sort type
		const char* NAME_VAR_NAME = "Name";				//!< Variable name for the name of the family
		const char* EIG_VAR_NAME = "Eigenvalues";		//!< Variable name for the eigenvalues
		const size_t DATA_WIDTH = 11;					//!< Number of elements in one row of member data

		void copyMe(const tpat_family_cr3bp&);
		std::vector<int> findMatches(double, std::vector<double>*) const;
		std::vector<tpat_family_member_cr3bp> getMatchingMember(double, std::vector<double>*, tpat_constraint) const;
		void getCoord(int, std::vector<double>*) const;

		void loadMemberData(mat_t*);
		void loadEigVals(mat_t*);
		void saveMembers(mat_t*);
		void saveMiscData(mat_t*);
		void saveEigVals(mat_t*);
};

#endif
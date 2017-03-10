/**
 *  \file Fam_cr3bp.hpp
 *	\brief 
 *	
 *	\author Andrew Cox
 *	\version May 25, 2016
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
#pragma once

#include "Core.hpp"
 
#include "FamMember_cr3bp.hpp"
#include "SysData_cr3bp.hpp"

#include "matio.h"

#include <cmath>
#include <string>
#include <vector>


namespace astrohelion{

// Forward declarations
class Constraint;


/**
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
	SORT_JC 	= 7,	//!< Sort by Jacobi Constant
	SORT_TOF 	= 6,	//!< Sort by Time-Of-Flight
	SORT_NONE 	= 999	//!< Do not adjust sorting; sortMembers() will do nothing
};
		
/**
 *	\ingroup traj cr3bp
 *	\brief A data structure to represent a family of orbits
 *
 *	Currently, this structure assumes family members are periodic.
 */
class Fam_cr3bp : public Core{
	public:
		/**
		 *  \name *structors
		 *  \{
		 */
		Fam_cr3bp(SysData_cr3bp);
		Fam_cr3bp(const char*);
		Fam_cr3bp(const Fam_cr3bp&);
		~Fam_cr3bp();
		//\}

		// Operators
		Fam_cr3bp& operator= (const Fam_cr3bp&);

		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		void addMember(FamMember_cr3bp);
		
		FamMember_cr3bp getMember(int) const;
		std::vector<FamMember_cr3bp> getMemberByStateVar(double, int) const;
		std::vector<FamMember_cr3bp> getMemberByTOF(double) const;
		std::vector<FamMember_cr3bp> getMemberByJacobi(double) const;
		std::string getName() const;
		int getNumMembers() const;
		FamSort_tp getSortType() const;
		const char* getSortTypeStr() const;
		SysData_cr3bp getSysData() const;
		SysData_cr3bp* getSysDataPtr();

		void setName(std::string);
		void setSortType(FamSort_tp);
		//\}

		std::vector<int> findBifurcations();
		void saveToMat(const char*);
		void sortEigs();
		void sortMembers();

	protected:
		std::string name = "NULL";								//!< Descriptive name of the family
		std::vector<FamMember_cr3bp> members {};				//!< Contains all family members
		SysData_cr3bp sysData = SysData_cr3bp();				//!< Describes the system this family exists in
		FamSort_tp sortType = FamSort_tp::SORT_NONE;			//!< Describes the most natural variable to sort family members by

		double matchTol = 1e-8;		//!< Acceptable tolerance (non-dim units) when locating a member by a specific attribute

		const char* DATA_VAR_NAME = "MemberData";		//!< Variable name for family member data
		const char* SORT_TYPE_VAR_NAME = "SortType";		//!< Variable name for the sort type
		const char* NAME_VAR_NAME = "Name";				//!< Variable name for the name of the family
		const char* EIG_VAR_NAME = "Eigenvalues";		//!< Variable name for the eigenvalues
		const char* EIGVEC_VAR_NAME = "Eigenvectors";	//!< Variable name for the eigenvectors
		const size_t DATA_WIDTH = 11;					//!< Number of elements in one row of member data

		void copyMe(const Fam_cr3bp&);
		std::vector<int> findMatches(double, std::vector<double>*) const;
		std::vector<FamMember_cr3bp> getMatchingMember(double, std::vector<double>*, Constraint) const;
		void getCoord(int, std::vector<double>*) const;

		void loadMemberData(mat_t*);
		void loadEigVals(mat_t*);
		void loadEigVecs(mat_t*);
		void loadSTMs(mat_t*);
		void saveMembers(mat_t*);
		void saveMiscData(mat_t*);
		void saveEigVals(mat_t*);
		void saveEigVecs(mat_t*);
		void saveSTMs(mat_t*);
};

}// END of Astrohelion namespace
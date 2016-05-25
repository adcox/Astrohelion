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

#include "tpat.hpp"
 
#include "tpat_famMember_cr3bp.hpp"
#include "tpat_sys_data_cr3bp.hpp"

#include "matio.h"

#include <cmath>
#include <string>
#include <vector>

// Forward declarations
class TPAT_Constraint;


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
enum class TPAT_SortFam_Tp : int{
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
 *	@brief A data structure to represent a family of orbits
 *
 *	Currently, this structure assumes family members are periodic.
 */
class TPAT_Fam_CR3BP : public TPAT{
	public:
		// *structors
		TPAT_Fam_CR3BP(TPAT_Sys_Data_CR3BP);
		TPAT_Fam_CR3BP(const char*);
		TPAT_Fam_CR3BP(const TPAT_Fam_CR3BP&);
		~TPAT_Fam_CR3BP();

		// Operators
		TPAT_Fam_CR3BP& operator= (const TPAT_Fam_CR3BP&);

		// Set and Get Methods
		void addMember(TPAT_FamMember_CR3BP);
		
		TPAT_FamMember_CR3BP getMember(int) const;
		std::vector<TPAT_FamMember_CR3BP> getMemberByStateVar(double, int) const;
		std::vector<TPAT_FamMember_CR3BP> getMemberByTOF(double) const;
		std::vector<TPAT_FamMember_CR3BP> getMemberByJacobi(double) const;
		std::string getName() const;
		int getNumMembers() const;
		TPAT_SortFam_Tp getSortType() const;
		const char* getSortTypeStr() const;
		TPAT_Sys_Data_CR3BP getSysData() const;
		TPAT_Sys_Data_CR3BP* getSysDataPtr();

		void setName(std::string);
		void setSortType(TPAT_SortFam_Tp);

		// Utility and Operational Functions
		std::vector<int> findBifurcations();
		void saveToMat(const char*);
		void sortEigs();
		void sortMembers();


	protected:
		std::string name = "NULL";								//!< Descriptive name of the family
		std::vector<TPAT_FamMember_CR3BP> members {};		//!< Contains all family members
		TPAT_Sys_Data_CR3BP sysData = TPAT_Sys_Data_CR3BP();	//!< Describes the system this family exists in
		TPAT_SortFam_Tp sortType = TPAT_SortFam_Tp::SORT_NONE;					//!< Describes the most natural variable to sort family members by

		double matchTol = 1e-8;		//!< Acceptable tolerance (non-dim units) when locating a member by a specific attribute

		const char* DATA_VAR_NAME = "MemberData";		//!< Variable name for family member data
		const char* SORT_TYPE_VAR_NAME = "SortType";		//!< Variable name for the sort type
		const char* NAME_VAR_NAME = "Name";				//!< Variable name for the name of the family
		const char* EIG_VAR_NAME = "Eigenvalues";		//!< Variable name for the eigenvalues
		const size_t DATA_WIDTH = 11;					//!< Number of elements in one row of member data

		void copyMe(const TPAT_Fam_CR3BP&);
		std::vector<int> findMatches(double, std::vector<double>*) const;
		std::vector<TPAT_FamMember_CR3BP> getMatchingMember(double, std::vector<double>*, TPAT_Constraint) const;
		void getCoord(int, std::vector<double>*) const;

		void loadMemberData(mat_t*);
		void loadEigVals(mat_t*);
		void saveMembers(mat_t*);
		void saveMiscData(mat_t*);
		void saveEigVals(mat_t*);
};

#endif
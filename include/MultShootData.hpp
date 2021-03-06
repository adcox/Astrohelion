/**
 *  @file MultShootData.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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

#include <iostream>
#include <map>
#include <vector>

#include "Constraint.hpp"
#include "Exceptions.hpp"
#include "Arcset.hpp"
#include "SysData.hpp"
#include "Arcset.hpp"



namespace astrohelion{

/**
 *  @brief Describe how time-of-flight data is stored in the free variable vector
 */
enum class MSTOF_tp : int{
	FIXED = 0,			//!< Times-of-flight are fixed and not part of the free variable vector
	VAR_FREE = 1,		//!< Times-of-flight are variable and part of the free variable vector
	VAR_EQUALARC = 2,	/*!< Times-of-flight are all equal; one TOF is stored in the free variable vector.
						 *	All times-of-flight must be of the same sign for this to work
						 */
	VAR_FIXSIGN = 3		/*!< Times-of-flight have same sign as initial guess; sqrt(TOF) is stored for 
						 * each arc
						 */
};

const char* MSTOF_tp_cStr(const MSTOF_tp&);
std::ostream& operator<<(std::ostream&, const MSTOF_tp&);

/**
 *  @brief The type of free variable being represented in the free-variable vector
 */
enum class MSVar_tp : int {EPOCH = 0,		//!< Epoch variable
							SLACK = 1,		//!< Slack variable (e.g., used in an inequality constraint)
							STATE = 2,		//!< State variable
							TOF = 3,		//!< Time-of-flight
							TOF_TOTAL = 4,	//!< Total time-of-flight for an arcset
							CTRL = 5		//!< Control variable
							};	

/**
 *  @brief The parent object of the variable
 *  @details "Parent" means that the variable is part of a larger object. For example,
 *  a node is the parent of the six state variables and the epoch associated with the node.
 *  Similarly, a segment is the parent of the time-of-flight variable along the segment.
 *  Constraints may own slack variables, and the entire arc may own quantities like total
 *  delta-V or total time-of-flight.
 */
enum class MSVarParent_tp : int {ARC = 0,	//!< The entire arc is the parent
								CON = 1,	//!< A constraint is the parent (for slack variables)
								NODE = 2,	//!< A node is the parent
								SEG = 3};	//!< A segment is the parent

/**
 *  @brief Represent a free variable in the free variable map
 *  @details A more complex object is required here to store
 *  information such as the variable type and ID.
 */
struct MSVarMap_Key{
	public:
		/**
		 * \name *structors
		 * \{
		 */
		MSVarMap_Key();
		MSVarMap_Key(MSVar_tp, int);
		MSVarMap_Key(const MSVarMap_Key&);
		//\}

		/**
		 * \name Operators
		 * \{
		 */
		MSVarMap_Key& operator =(const MSVarMap_Key&);
		bool friend operator <(const MSVarMap_Key&, const MSVarMap_Key&);
		//\}

		/**
		 * \name Utility Functions
		 * \{
		 */
		static const char* type2str(MSVar_tp);
		//\}

		MSVar_tp type = MSVar_tp::STATE;	//!< Type of variable
		int id = -1;						//!< ID of the parent object
	private:
		/**
		 * \name Utility Functions
		 * \{
		 */
		void copyMe(const MSVarMap_Key&);
		//\}
};

/**
 *  @brief Represent a free variable in the free variable vector
 *  @details This object stores information such as the position of
 *  the variable within the free-variable vector, the number of elements
 *  (i.e., rows) required by the variable (many variables are vector quantities).
 *  A parent object and map key are also stored.
 */
struct MSVarMap_Obj{
	public:
		/**
		 * \name *structors
		 * \{
		 */
		MSVarMap_Obj();
		MSVarMap_Obj(MSVar_tp);
		MSVarMap_Obj(MSVar_tp, int, int, int nRows = 1);
		MSVarMap_Obj(MSVarMap_Key, int, int nRows = 1);
		MSVarMap_Obj(const MSVarMap_Obj&);
		//\}

		/**
		 * \name Operators
		 * \{
		 */
		MSVarMap_Obj& operator =(const MSVarMap_Obj&);
		//\}

		/**
		 * \name Utility Functions
		 * \{
		 */
		bool matches(MSVar_tp, int) const;
		static const char* parent2str(MSVarParent_tp);
		//\}

		MSVarMap_Key key;					//!< Identifies this object by variable type and parent ID
		MSVarParent_tp parent = MSVarParent_tp::NODE;		//!< Object type that owns the represented variable
		int row0 = -1;		//!< Index of the first row of the free variable vector this variable occupies
		int nRows = -1;		//!< Number of rows of the free variable vector this variable occupies
	private:
		/**
		 * \name Utility Functions
		 * \{
		 */
		void copyMe(const MSVarMap_Obj&);
		//\}

		/**
		 * \name *structors
		 * \{
		 */
		void init();
		//\}
};

/**		
 *	@brief a custom data class to encapsulate data used in each iteration
 *	of the corrections process.
 *
 *	This data object can be passed to other functions, allowing us to break 
 *	the master corrections loop into smaller functions without requiring an
 *	obscene amount of arguments to be passed in.
 *	
 *	@author Andrew Cox
 *	@version May 16, 2016
 *	@copyright GNU GPL v3.0
 */
class MultShootData{
	public:

		/**
		 * \name *structors
		 * \{
		 */
		MultShootData(const Arcset*);
		MultShootData(const MultShootData&);
		//\}

		/**
		 * \name Operators
		 * \{
		 */
		MultShootData& operator =(const MultShootData&);
		//\}

		/**
		 * \name Set and Get Functions
		 * \{
		 */
		MSVarMap_Obj getVarMap_obj(MSVar_tp, int) const;
		//\}

		// Utilities

		// Variables
		const Arcset *pArcIn;						//!< A pointer to the arcset input for this corrections process
		Arcset *pArcOut = nullptr;					//!< A pointer to arcset that will be output
		std::vector<Arcset> propSegs {};			//!< A collection of all propagated segments; index here matches index within the input arcset
		std::vector<double> X0 {};					//!< Initial, uncorrected free-variable vector
		std::vector<double> X {};					//!< Free-Variable Vector
		std::vector<double> FX {};					//!< Constraint Function Vector
		std::vector<Tripletd> DF_elements {};		//!< Triplet vector that contains all nonzero elements of the DF matrix
		std::vector<double> deltaVs {};				//!< nx3 vector of non-dim delta-Vs
		std::vector<Constraint> allCons {};			//!< A list of all constraints
		std::map<MSVarMap_Key, MSVarMap_Obj> freeVarMap {};	//!< Structure that maps free variables to their rows in the free variable vector
		std::vector<int> slackAssignCon {};			//!< Indices of constraints, index of entry corresponds to a slack variable
		std::vector<int> conRows {};				//!< Each entry holds the row # for the constraint; i.e. 0th element holds row # for 0th constraint

		unsigned int numNodes = 0;							//!< Number of nodes in the entire arcset
		int count = 0;								//!< Count of number of iterations through corrections process

		int numSlack = 0;							//!< number of slack variables
		int totalCons = 0;							//!< Total # constraints -> # rows of DF
		int totalFree = 0;							//!< Total # free var. -> # cols of DF

		MSTOF_tp tofTp = MSTOF_tp::VAR_FREE;		//!< Describes how times-of-flight are stored in the free-variable vector
	protected:
		/**
		 * \name Utility Functions
		 * \{
		 */
		void copyMe(const MultShootData&);
		//\}
};

}// END of Astrohelion namespace
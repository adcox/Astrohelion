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
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "Constraint.hpp"
#include "Exceptions.hpp"
#include "Nodeset.hpp"
#include "SysData.hpp"
#include "Traj.hpp"

#include <map>
#include <vector>

namespace astrohelion{
// Forward Declarations

enum class MSVarType : int {EPOCH = 0,
							SLACK = 1,
							STATE = 2,
							TOF = 3,
							TOF_TOTAL = 4};

enum class MSVarParent : int {	ARC = 0,
								CON = 1,
								NODE = 2,
								SEG = 3};
struct MSVarMap_Key{
	public:

		MSVarType type = MSVarType::STATE;	//!< Type of variable
		int id = -1;				//!< ID of the parent object

		MSVarMap_Key();
		MSVarMap_Key(MSVarType, int);
		MSVarMap_Key(const MSVarMap_Key&);
		
		MSVarMap_Key& operator =(const MSVarMap_Key&);
		bool friend operator <(const MSVarMap_Key&, const MSVarMap_Key&);
		static const char* type2str(MSVarType);

	private:
		void copyMe(const MSVarMap_Key&);
};

struct MSVarMap_Obj{
	public:

		MSVarMap_Key key;					//!< Identifies this object by variable type and parent ID
		MSVarParent parent = MSVarParent::NODE;		//!< Object type that owns the represented variable
		int row0 = -1;		//!< Index of the first row of the free variable vector this variable occupies
		int nRows = -1;		//!< Number of rows of the free variable vector this variable occupies

		MSVarMap_Obj();
		MSVarMap_Obj(MSVarType);
		MSVarMap_Obj(MSVarType, int, int);
		MSVarMap_Obj(MSVarMap_Key, int);
		MSVarMap_Obj(const MSVarMap_Obj&);

		MSVarMap_Obj& operator =(const MSVarMap_Obj&);

		bool matches(MSVarType, int) const;
		
	private:
		void copyMe(const MSVarMap_Obj&);
		void init();
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

		// *structors
		MultShootData(const Nodeset*);
		MultShootData(const MultShootData&);

		// Operators
		MultShootData& operator =(const MultShootData&);

		// Set and Get
		MSVarMap_Obj getVarMap_obj(MSVarType, int) const;

		// Utilities

		// Variables
		const SysData *sysData;				//!< A pointer to the system data object used for this corrections process
		const Nodeset *nodeset;				//!< A pointer to the nodeset input for this corrections process
		std::vector<double> X0 {};					//!< Initial, uncorrected free-variable vector
		std::vector<double> X {};					//!< Free-Variable Vector
		std::vector<double> FX {};					//!< Constraint Function Vector
		std::vector<double> DF {};					//!< Jacobian Matrix
		std::vector<double> deltaVs {};				//!< nx3 vector of non-dim delta-Vs
		std::vector<Traj> propSegs {};			//!< A collection of all propagated segments
		std::vector<Constraint> allCons {};	//!< A list of all constraints
		std::map<MSVarMap_Key, MSVarMap_Obj> freeVarMap {};	//!< Structure that maps free variables to their rows in the free variable vector
		std::vector<int> slackAssignCon {};			//!< Indices of constraints, index of entry corresponds to a slack variable
		std::vector<int> conRows {};				//!< Each entry holds the row # for the constraint; i.e. 0th element holds row # for 0th constraint
		
		/**
		 * A scalar coefficient for each free variable type to scale them to the appropriate magnitude.
		 * These scalings should make numerical processes more successful but must be reversed before
		 * the free variables are passed into various computation functions (e.g. simulations). Thus,
		 * every constraint computation function (stored in the models) is responsible for reversing any
		 * scaling necessary to accurately compute constraints and partial derivatives.
		 */
		std::vector<double> freeVarScale {};

		int numNodes = 0;			//!< Number of nodes in the entire nodeset
		int count = 0;				//!< Count of number of iterations through corrections process

		int numSlack = 0;			//!< number of slack variables
		int totalCons = 0;			//!< Total # constraints -> # rows of DF
		int totalFree = 0;			//!< Total # free var. -> # cols of DF

		bool varTime = true;		//!< Whether or not the corrector is using variable time
		bool equalArcTime = false;	//!< Whether or not each arc must have an equal duration
	
	protected:
		void copyMe(const MultShootData&);
};

}// END of Astrohelion namespace
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

#ifndef H_TPAT_MULTSHOOT_DATA
#define H_TPAT_MULTSHOOT_DATA

#include "tpat_constraint.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_nodeset.hpp"
#include "tpat_sys_data.hpp"
#include "tpat_traj.hpp"

#include <vector>

// Forward Declarations

struct ms_varMap_obj{
public:

	enum ms_varMap_tp {EPOCH, SLACK, STATE, TOF, TOF_TOTAL};
	enum ms_varMap_parent {ARC, CON, NODE, SEG};

	ms_varMap_tp type = STATE;			//!< Type of variable
	ms_varMap_parent parent = NODE;		//!< Object type that owns the represented variable
	int row0 = -1;		//!< Index of the first row of the free variable vector this variable occupies
	int nRows = -1;		//!< Number of rows of the free variable vector this variable occupies
	int id = -1;		//!< ID of the parent object

	ms_varMap_obj(ms_varMap_tp t) : type(t){
		init();
	}//================================================
	
	ms_varMap_obj(ms_varMap_tp t, int row0, int id) : type(t){
		this->row0 = row0;
		this->id = id;
		init();
	}//================================================

	ms_varMap_obj(const ms_varMap_obj &obj){
		copyMe(obj);
	}//================================================

	ms_varMap_obj& operator =(const ms_varMap_obj &obj){
		copyMe(obj);
		return *this;
	}//================================================

	bool matches(ms_varMap_tp t, int id) const{
		return type == t && this->id == id;
	}//================================================

	static const char* type2str(ms_varMap_tp tp){
		switch(tp){
			case EPOCH: return "EPOCH"; break;
			case SLACK: return "SLACK"; break;
			case STATE: return "STATE"; break;
			case TOF: return "TOF"; break;
			case TOF_TOTAL: return "TOF_TOTAL"; break;
		}
		return "Unrecognized Type";
	}//===============================================
	
private:
	void init(){
		nRows = 1;
		switch(type){
			case STATE:
				parent = NODE;
				nRows = 6;
				break;
			case EPOCH: parent = NODE; break;
			case TOF: parent = SEG; break;
			case TOF_TOTAL: parent = ARC; break;
			case SLACK: parent = CON; break;
			default: throw tpat_exception("ms_varMap_obj constructor: Unrecognized type");
		}
	}//================================================

	void copyMe(const ms_varMap_obj &obj){
		type = obj.type;
		parent = obj.parent;
		row0 = obj.row0;
		nRows = obj.nRows;
		id = obj.id;
	}//================================================
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
class tpat_multShoot_data{
	public:

		// *structors
		tpat_multShoot_data(const tpat_nodeset*);
		tpat_multShoot_data(const tpat_multShoot_data&);

		// Operators
		tpat_multShoot_data& operator =(const tpat_multShoot_data&);

		// Set and Get
		ms_varMap_obj getVarMap_obj(ms_varMap_obj::ms_varMap_tp, int) const;

		// Utilities

		// Variables
		const tpat_sys_data *sysData;				//!< A pointer to the system data object used for this corrections process
		const tpat_nodeset *nodeset;				//!< A pointer to the nodeset input for this corrections process
		std::vector<double> X0 {};					//!< Initial, uncorrected free-variable vector
		std::vector<double> X {};					//!< Free-Variable Vector
		std::vector<double> FX {};					//!< Constraint Function Vector
		std::vector<double> DF {};					//!< Jacobian Matrix
		std::vector<double> deltaVs {};				//!< nx3 vector of non-dim delta-Vs
		std::vector<tpat_traj> propSegs {};			//!< A collection of all propagated segments
		std::vector<tpat_constraint> allCons {};	//!< A list of all constraints
		std::vector<ms_varMap_obj> freeVarMap {};	//!< Structure that maps free variables to their rows in the free variable vector
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
		void copyMe(const tpat_multShoot_data&);
};

#endif
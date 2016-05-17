/**
 *  @file 
 *	@brief 
 *
 *	@author Andrew Cox
 *	@version 
 *	@copyright GNU GPL v3.0
 */
 
/*
 *  Trajectory Propagation and Analysis Toolkit 
 *  Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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

#include "tpat_multShoot_data.hpp"

#include "tpat_utilities.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *  @brief Construct a new tpat_multShoot_data object
 *  @param set pointer to the nodeset being corrected
 */
tpat_multShoot_data::tpat_multShoot_data(const tpat_nodeset *set) : sysData(set->getSysData()), nodeset(set){}

/**
 *  @brief Copy constructor
 *  @param it reference to another tpat_multShoot_data object
 */
tpat_multShoot_data::tpat_multShoot_data(const tpat_multShoot_data &it) : sysData(it.nodeset->getSysData()), nodeset(it.nodeset){
	copyMe(it);
}//============================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *  @brief Assignment operator
 *  @param it reference to another tpat_multShoot_data object
 */
tpat_multShoot_data& tpat_multShoot_data::operator =(const tpat_multShoot_data &it){
	copyMe(it);
	sysData = it.sysData;	// Copying ADDRESS!
	nodeset = it.nodeset;	// Copying ADDRESS!
	return *this;
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *  @brief Retrieve a free variable map object by type and reference ID
 * 
 *  @param type the type of object this map object represents (i.e., EPOCH, SLACK, STATE, TOF, TOF_TOTAL)
 *  @param refID the ID of the parent object
 * 
 *  @return the variable map object that represents the desired variable
 *  @see ms_varMap_type
 *  @throws tpat_exception if the object cannot be located
 *  @todo This function could be greatly sped up by leveraging a hash table
 */
ms_varMap_obj tpat_multShoot_data::getVarMap_obj(ms_varMap_obj::ms_varMap_tp type, int refID) const{
	for(size_t i = 0; i < freeVarMap.size(); i++){
		if(freeVarMap[i].matches(type, refID))
			return freeVarMap[i];
	}
	printErr("Trying to locate %s object with id %d\n", ms_varMap_obj::type2str(type), refID);
	printf("freeVarMap = {\n");
	for(size_t i = 0; i < freeVarMap.size(); i++){
		printf("  %02zu: %s obj, ID = %d\n", i, ms_varMap_obj::type2str(freeVarMap[i].type), freeVarMap[i].id);
	}
	printf("}\n");
	throw tpat_exception("tpat_multShoot_data::getVarMap_obj: Could not locate object");
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *  @brief Copy all parameters from one tpat_multShoot_data object to another
 *  @param it reference to source tpat_multShoot_data object
 */
void tpat_multShoot_data::copyMe(const tpat_multShoot_data &it){
	
	allCons = it.allCons;
	conRows = it.conRows;
	count = it.count;
	deltaVs = it.deltaVs;
	DF = it.DF;
	equalArcTime = it.equalArcTime;
	freeVarMap = it.freeVarMap;
	freeVarScale = it.freeVarScale;
	FX = it.FX;
	numNodes = it.numNodes;
	numSlack = it.numSlack;
	primPos = it.primPos;
	primVel = it.primVel;
	propSegs = it.propSegs;
	slackAssignCon = it.slackAssignCon;
	totalCons = it.totalCons;
	totalFree = it.totalFree;
	varTime = it.varTime;
	X = it.X;
	X0 = it.X0;
}//============================================
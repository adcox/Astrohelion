/**
 *  @file tpat_multShoot_data.cpp
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


//-----------------------------------------------------------------------------
//  ** MSVarMap_Key Functions
//-----------------------------------------------------------------------------


MSVarMap_Key::MSVarMap_Key(){}

MSVarMap_Key::MSVarMap_Key(MSVarType type, int id){
	this->type = type;
	this->id = id;
}
MSVarMap_Key::MSVarMap_Key(const MSVarMap_Key &k){ copyMe(k); }

MSVarMap_Key& MSVarMap_Key::operator =(const MSVarMap_Key &k){
	copyMe(k);
	return *this;
}//================================================

bool operator <(const MSVarMap_Key &lhs, const MSVarMap_Key &rhs){
	bool isLess =  lhs.type < rhs.type || (lhs.type == rhs.type && lhs.id < rhs.id);

	// printf("Comparing %s (%d) < %s (%d) : %s\n", type2str(lhs.type),
	// 	lhs.id, type2str(rhs.type), rhs.id, isLess ? "TRUE" : "FALSE");
	
	return isLess;
}//================================================

void MSVarMap_Key::copyMe(const MSVarMap_Key &k){
	type = k.type;
	id = k.id;
}//================================================

const char* MSVarMap_Key::type2str(MSVarType tp){
	switch(tp){
		case MSVarType::EPOCH: return "EPOCH"; break;
		case MSVarType::SLACK: return "SLACK"; break;
		case MSVarType::STATE: return "TPAT_Constraint_Tp::STATE"; break;
		case MSVarType::TOF: return "TOF"; break;
		case MSVarType::TOF_TOTAL: return "TOF_TOTAL"; break;
	}
	return "Unrecognized Type";
}//===============================================



//-----------------------------------------------------------------------------
//  ** MSVarMap_Obj Functions
//-----------------------------------------------------------------------------


MSVarMap_Obj::MSVarMap_Obj() : key() {}

MSVarMap_Obj::MSVarMap_Obj(MSVarType t) : key(){
	key.type = t;
	init();
}//================================================

MSVarMap_Obj::MSVarMap_Obj(MSVarType t, int row0, int id) : key(){
	key.type = t;
	key.id = id;
	this->row0 = row0;
	init();
}//================================================

MSVarMap_Obj::MSVarMap_Obj(MSVarMap_Key k, int rowNum) : key(k){
	this->row0 = rowNum;
	init();
}//================================================

MSVarMap_Obj::MSVarMap_Obj(const MSVarMap_Obj &obj) : key(){
	copyMe(obj);
}//================================================

MSVarMap_Obj& MSVarMap_Obj::operator =(const MSVarMap_Obj &obj){
	copyMe(obj);
	return *this;
}//================================================

bool MSVarMap_Obj::matches(MSVarType t, int id) const{
	return key.type == t && this->key.id == id;
}//================================================

void MSVarMap_Obj::init(){
	nRows = 1;
	switch(key.type){
		case MSVarType::STATE:
			parent = MSVarParent::NODE;
			nRows = 6;
			break;
		case MSVarType::EPOCH: parent = MSVarParent::NODE; break;
		case MSVarType::TOF: parent = MSVarParent::SEG; break;
		case MSVarType::TOF_TOTAL: parent = MSVarParent::ARC; break;
		case MSVarType::SLACK: parent = MSVarParent::CON; break;
		default: throw TPAT_Exception("MSVarMap_Obj constructor: Unrecognized type");
	}
}//================================================

void MSVarMap_Obj::copyMe(const MSVarMap_Obj &obj){
	key = obj.key;
	parent = obj.parent;
	row0 = obj.row0;
	nRows = obj.nRows;
}//================================================

//-----------------------------------------------------------------------------
//  ** TPAT_MULTSHOOT_DATA Functions
//-----------------------------------------------------------------------------


//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *  @brief Construct a new TPAT_MultShoot_Data object
 *  @param set pointer to the nodeset being corrected
 */
TPAT_MultShoot_Data::TPAT_MultShoot_Data(const TPAT_Nodeset *set) : sysData(set->getSysData()), nodeset(set){
	numNodes = set->getNumNodes();
}//====================================================

/**
 *  @brief Copy constructor
 *  @param it reference to another TPAT_MultShoot_Data object
 */
TPAT_MultShoot_Data::TPAT_MultShoot_Data(const TPAT_MultShoot_Data &it) : sysData(it.sysData), nodeset(it.nodeset){
	copyMe(it);
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *  @brief Assignment operator
 *  @param it reference to another TPAT_MultShoot_Data object
 */
TPAT_MultShoot_Data& TPAT_MultShoot_Data::operator =(const TPAT_MultShoot_Data &it){
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
 *  @param type the type of object this map object represents (i.e., EPOCH, SLACK, TPAT_Constraint_Tp::STATE, TOF, TOF_TOTAL)
 *  @param refID the ID of the parent object
 * 
 *  @return the variable map object that represents the desired variable
 *  @see ms_varMap_type
 *  @throws TPAT_Exception if the object cannot be located
 *  @todo This function could be greatly sped up by leveraging a hash table
 */
MSVarMap_Obj TPAT_MultShoot_Data::getVarMap_obj(MSVarType type, int refID) const{
	return freeVarMap.at(MSVarMap_Key(type, refID));
	// for(size_t i = 0; i < freeVarMap.size(); i++){
	// 	if(freeVarMap[i].matches(type, refID))
	// 		return freeVarMap[i];
	// }
	// printErr("Trying to locate %s object with id %d\n", MSVarMap_Obj::type2str(type), refID);
	// printf("freeVarMap = {\n");
	// for(size_t i = 0; i < freeVarMap.size(); i++){
	// 	printf("  %02zu: %s obj, ID = %d\n", i, MSVarMap_Obj::type2str(freeVarMap[i].key.type), freeVarMap[i].key.id);
	// }
	// printf("}\n");
	// throw TPAT_Exception("TPAT_MultShoot_Data::getVarMap_obj: Could not locate object");
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *  @brief Copy all parameters from one TPAT_MultShoot_Data object to another
 *  @param it reference to source TPAT_MultShoot_Data object
 */
void TPAT_MultShoot_Data::copyMe(const TPAT_MultShoot_Data &it){
	
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
	propSegs = it.propSegs;
	slackAssignCon = it.slackAssignCon;
	totalCons = it.totalCons;
	totalFree = it.totalFree;
	varTime = it.varTime;
	X = it.X;
	X0 = it.X0;				
}//============================================
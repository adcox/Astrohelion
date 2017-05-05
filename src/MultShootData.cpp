/**
 *  \file MultShootData.cpp
 *	\brief Contains member functions for several classes that support
 *	multiple shooting data manipulation and transfer
 *
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of Astrohelion
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

#include "MultShootData.hpp"

#include "Utilities.hpp"


namespace astrohelion{
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------	
//-----------------------------------------------------------------------------
//  ** MSVarMap_Key Functions
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
/**
 *  \brief Construct a default Multiple Shooting Variable Map Key object
 */
MSVarMap_Key::MSVarMap_Key(){}

/**
 *  \brief Construct a Multiple Shooting Variable Map Key object
 *  \details [long description]
 * 
 *  \param type Variable type
 *  \param id ID associated with the variable (e.g., Node ID for 
 *  epochs and states, Segment ID for TOF, etc.)
 */
MSVarMap_Key::MSVarMap_Key(MSVar_tp type, int id){
	this->type = type;
	this->id = id;
}//================================================

/**
 *  \brief Copy constructory
 *  \param k another MSVarMap_Key object
 */
MSVarMap_Key::MSVarMap_Key(const MSVarMap_Key &k){ copyMe(k); }

/**
 *  \brief Assignment operator
 *  \details [long description]
 * 
 *  \param k Another MSVarMap_Key object
 *  \return reference to this object after assignment to the other object
 */
MSVarMap_Key& MSVarMap_Key::operator =(const MSVarMap_Key &k){
	copyMe(k);
	return *this;
}//================================================

/**
 *  \brief Comparator for map operations
 *  \details Results are arbitrary; used only in internal
 *  map sorting processes
 * 
 *  \param lhs left-hand-side object reference
 *  \param rhs right-hand-side object reference
 * 
 *  \return TRUE if lhs.type < rhs.type or if lhs.type == rhs.type and lhs.id < rhs.id
 */
bool operator <(const MSVarMap_Key &lhs, const MSVarMap_Key &rhs){
	bool isLess =  lhs.type < rhs.type || (lhs.type == rhs.type && lhs.id < rhs.id);

	// printf("Comparing %s (%d) < %s (%d) : %s\n", type2str(lhs.type),
	// 	lhs.id, type2str(rhs.type), rhs.id, isLess ? "TRUE" : "FALSE");
	
	return isLess;
}//================================================

/**
 *  \brief Utility function to copy all attributes of one object to another
 *  \param k Reference to another MSVarMap_Key object
 */
void MSVarMap_Key::copyMe(const MSVarMap_Key &k){
	type = k.type;
	id = k.id;
}//================================================

/**
 *  \brief Retrieve a string to describe the type
 *  \details [long description]
 * 
 *  \param tp Type associated with a MSVarMap_Key object
 *  \return a string to describe the type
 */
const char* MSVarMap_Key::type2str(MSVar_tp tp){
	switch(tp){
		case MSVar_tp::EPOCH: return "EPOCH"; break;
		case MSVar_tp::SLACK: return "SLACK"; break;
		case MSVar_tp::STATE: return "STATE"; break;
		case MSVar_tp::TOF: return "TOF"; break;
		case MSVar_tp::TOF_TOTAL: return "TOF_TOTAL"; break;
	}
	return "Unrecognized Type";
}//===============================================


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//  ** MSVarMap_Obj Functions
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


/**
 *  \brief Construct a default Multiple Shooting Variable Map object
 */
MSVarMap_Obj::MSVarMap_Obj() : key() {}

/**
 *  \brief Construct a Multiple Shooting Variable Map object
 *  \param t Variable type
 */
MSVarMap_Obj::MSVarMap_Obj(MSVar_tp t) : key(){
	key.type = t;
	init();
}//================================================

/**
 *  \brief Construct a Multiple Shooting Variable Map object
 *  \details [long description]
 * 
 *  \param t Variable type
 *  \param row0 index of the first row this object occupies in the free variable vector.
 *  If row0 is set equal to -1, the variable is not stored in the free variable vector and
 *  is instead stored in the input nodeset (and not varied).
 *  \param id ID associated with the object that stores this object information (e.g.,
 *  a Node ID for state and epoch variables, or a segment ID for a time-of-flight)
 *  \param nRows Number of rows this object occupies in the free variable vector
 */
MSVarMap_Obj::MSVarMap_Obj(MSVar_tp t, int row0, int id, int nRows) : key(){
	key.type = t;
	key.id = id;
	this->row0 = row0;
	this->nRows = nRows;
	init();
}//====================================================

/**
 *  \brief Construct a Multiple Shooting Variable Map object
 *  \details [long description]
 * 
 *  \param k MSVarMap_Key object that stores information about the variable type and ID
 *  \param row0 index of the first row this object occupies in the free variable vector.
 *  If row0 is set equal to -1, the variable is not stored in the free variable vector and
 *  is instead stored in the input nodeset (and not varied).
 *  \param nRows Number of rows this object occupies in the free variable vector
 */
MSVarMap_Obj::MSVarMap_Obj(MSVarMap_Key k, int row0, int nRows) : key(k){
	this->row0 = row0;
	this->nRows = nRows;
	init();
}//================================================

/**
 *  \brief Copy constructor
 *  \param obj reference to another MSVarMap_Obj object
 */
MSVarMap_Obj::MSVarMap_Obj(const MSVarMap_Obj &obj) : key(){
	copyMe(obj);
}//================================================

/**
 *  \brief Assignment operator
 *  \param obj reference to another object
 *  \return a reference to this object after the assignment is completed
 */
MSVarMap_Obj& MSVarMap_Obj::operator =(const MSVarMap_Obj &obj){
	copyMe(obj);
	return *this;
}//================================================

/**
 *  \brief Determine whether or not this object matches a specified
 *  type and ID
 * 
 *  \param t Variable type
 *  \param id Variable ID
 * 
 *  \return TRUE if this object is associated with a MSVarMap_Key that matches the
 *  specified type and ID
 */
bool MSVarMap_Obj::matches(MSVar_tp t, int id) const{
	return key.type == t && this->key.id == id;
}//================================================

/**
 *  \brief Initialize the object after construction
 *  \details Set the parent based on key/variable type
 */
void MSVarMap_Obj::init(){
	switch(key.type){
		case MSVar_tp::STATE:
			parent = MSVarParent_tp::NODE;
			break;
		case MSVar_tp::EPOCH:
			parent = MSVarParent_tp::NODE;
			break;
		case MSVar_tp::TOF:
			parent = MSVarParent_tp::SEG;
			break;
		case MSVar_tp::TOF_TOTAL:
			parent = MSVarParent_tp::ARC;
			break;
		case MSVar_tp::SLACK:
			parent = MSVarParent_tp::CON;
			break;
		default: throw Exception("MSVarMap_Obj constructor: Unrecognized type");
	}
}//================================================

/**
 *  \brief Copy all attributes from one object to this one
 *  \param obj reference to another MSVarMap_Obj object
 */
void MSVarMap_Obj::copyMe(const MSVarMap_Obj &obj){
	key = obj.key;
	parent = obj.parent;
	row0 = obj.row0;
	nRows = obj.nRows;
}//================================================

/**
 *  \brief Retrieve a string that describes the parent type
 *  \param p Parent type
 *  \return a string that describes the parent type
 */
const char* MSVarMap_Obj::parent2str(MSVarParent_tp p){
	switch(p){
		case MSVarParent_tp::ARC: return "ENTIRE ARC";
		case MSVarParent_tp::CON: return "CONSTRAINT";
		case MSVarParent_tp::NODE: return "NODE";
		case MSVarParent_tp::SEG: return "SEG";
		default: return "Unrecognized Type";
	}
}//================================================

//-----------------------------------------------------------------------------
//  ** MultShootData Functions
//-----------------------------------------------------------------------------


//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *  \brief Construct a new MultShootData object
 *  \param set pointer to the nodeset being corrected
 */
MultShootData::MultShootData(const Arcset *set) : nodesIn(set){
	numNodes = set->getNumNodes();
}//====================================================

/**
 *  \brief Copy constructor
 *  \param it reference to another MultShootData object
 */
MultShootData::MultShootData(const MultShootData &it) : nodesIn(it.nodesIn){
	copyMe(it);
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *  \brief Assignment operator
 *  \param it reference to another MultShootData object
 */
MultShootData& MultShootData::operator =(const MultShootData &it){
	copyMe(it);
	return *this;
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *  \brief Retrieve a free variable map object by type and reference ID
 * 
 *  \param type the type of object this map object represents (i.e., EPOCH, SLACK, Constraint_tp::STATE, TOF, TOF_TOTAL)
 *  \param refID the ID of the parent object
 * 
 *  \return the variable map object that represents the desired variable
 *  @see ms_varMap_type
 *  \throws Exception if the object cannot be located
 *  @todo This function could be greatly sped up by leveraging a hash table
 */
MSVarMap_Obj MultShootData::getVarMap_obj(MSVar_tp type, int refID) const{
	if(freeVarMap.count(MSVarMap_Key(type, refID)) == 0){
		printErr("Attempted to access MSVar_tp %s, refID %d\n", MSVarMap_Key::type2str(type), refID);
		// printf("freeVarMap = {\n");
		// for(auto& obj : freeVarMap){
		// 	printf("  %s obj, ID = %d\n", MSVarMap_Key::type2str(obj.first.type), obj.first.id);
		// }
		// printf("}\n");
		throw Exception("MultShootData::getVarMap_obj: Variable map object does not exist\n");
	}
	return freeVarMap.at(MSVarMap_Key(type, refID));
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *  \brief Copy all parameters from one MultShootData object to another
 *  \param it reference to source MultShootData object
 */
void MultShootData::copyMe(const MultShootData &it){
	
	allCons = it.allCons;
	conRows = it.conRows;
	count = it.count;
	deltaVs = it.deltaVs;
	DF_elements = it.DF_elements;
	bEqualArcTime = it.bEqualArcTime;
	freeVarMap = it.freeVarMap;
	FX = it.FX;
	numNodes = it.numNodes;
	numSlack = it.numSlack;
	slackAssignCon = it.slackAssignCon;
	totalCons = it.totalCons;
	totalFree = it.totalFree;
	bVarTime = it.bVarTime;
	X = it.X;
	X0 = it.X0;
	propSegs = it.propSegs;

	nodesIn = it.nodesIn;		// Copying address
	nodesOut = it.nodesOut;		// Copying address			
}//============================================


}// END of Astrohelion namespace
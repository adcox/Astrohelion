/** 
 *	@file tpat_sys_data.cpp
 *
 *	System Data (Abstract Base Class)
 *
 * 	Stores information like the number and names of primaries, and the value 
 * 	of different characteristic quantities.
 *
 *	Note that because this class contains pure virtual functions, it cannot be 
 *	instantiated as an object. It CAN be used as a pointer, but all instances must 
 *	be one of the derived system data classes.
 */
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
#include "tpat.hpp"

#include "tpat_sys_data.hpp"
#include "tpat_utilities.hpp"
 
/**
 *	@brief Create a new system data object and initialize all values to zero.
 */
tpat_sys_data::tpat_sys_data(){}

/**
 *	@brief Copy constructor
 *
 *	@param d a system data object
 */
tpat_sys_data::tpat_sys_data(const tpat_sys_data& d){
	copyData(d);
}//============================================

/**
 *	@brief Destructor
 */
tpat_sys_data::~tpat_sys_data(){}

/**
 *	@brief Copy the system data object
 *
 *	@param d a system data object
 *	@return this data object, set to equal the input object
 */
tpat_sys_data& tpat_sys_data::operator =(const tpat_sys_data &d){
	copyData(d);
	return *this;
}//==========================================

/**
 *	@brief compare two system data objects. They are the same if they 
 *	represent the same system (exact same primaries)
 *
 *	@param lhs
 *	@param rhs
 *	@return whether or not they are the same
 */
bool operator ==(const tpat_sys_data &lhs, const tpat_sys_data &rhs){
	if(lhs.type != rhs.type)
		return false;

	if(lhs.numPrimaries != rhs.numPrimaries)
		return false;

	for(int p = 0; p < lhs.numPrimaries; p++){
		if(lhs.primIDs[p] != rhs.primIDs[p])
			return false;
	}
	
	return true;
}//===========================================

/**
 *	@brief compare two system data objects. They are the same if they 
 *	represent the same system (exact same primaries)
 *
 *	@param lhs
 *	@param rhs
 *	@return whether or not they are different
 */
bool operator !=(const tpat_sys_data &lhs, const tpat_sys_data &rhs){
	return ! operator==(lhs, rhs);
}//===========================================

void tpat_sys_data::copyData(const tpat_sys_data &d){
	charL = d.charL;
	charT = d.charT;
	charM = d.charM;
	type = d.type;
	primaries = d.primaries;
	primIDs = d.primIDs;
	otherParams = d.otherParams;
}//==========================================

/**
 *	@return the number of primaries this system models
 */
int tpat_sys_data::getNumPrimaries() const { return primaries.size(); }

/**
 *	@brief Retrieve the name of one of the system primaries
 *	@param n the "index" of the primary, starts at 0
 *	@return the name of the n'th primary
 */
std::string tpat_sys_data::getPrimary(int n) const{ 
	try{
		return primaries.at(n);
	}catch(std::out_of_range &err){
		printErr("Could not find primary at index %d; out of range\n", n);
		return "NULL";
	}
}//=============================

/**
 *	@brief Retrieve a unique numerical ID for one of the system primaries
 *	@param n the index of the primary, starts at 0
 *	@return a unique numerical ID for this primary; useful for comparing systems
 */
int tpat_sys_data::getPrimID(int n) const{ return primIDs.at(n); }

/**
 *	@return the characteristic length (km) associated with this system
 */
double tpat_sys_data::getCharL() const { return charL; }

/**
 *	@return the charactersitic mass (kg) associated with this system
 */
double tpat_sys_data::getCharM() const { return charM; }

/**
 *	@return the characteristic time (s) associated with this system
 */
double tpat_sys_data::getCharT() const { return charT; }

/**
 *	@return the system_t type associated with this system
 */
tpat_sys_data::system_t tpat_sys_data::getType() const { return type; }

/**
 *	@return a string (human-readable) version of the system type
 */
std::string tpat_sys_data::getTypeStr() const{
	switch (type){
		case tpat_sys_data::UNDEF_SYS:
			return "Undefined System Type";
		case tpat_sys_data::CR3BP_SYS:
			return "CR3BP";
		case tpat_sys_data::CR3BP_LTVP_SYS:
			return "CR3BP, Low Thrust, Velocity-Pointing";
		case tpat_sys_data::BCR4BPR_SYS:
			return "BCR4BP, Rotating Coord.";
		default:
			return "Unrecognized value... May be a bug";
	}
}//=================================================


//
/** 
 *	@file tpat_sys_data.cpp
 *
 *	@brief System Data (Abstract Base Class)
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

#include "tpat_sys_data.hpp"
#include "tpat_utilities.hpp"
 
/**
 *	@brief Create a new system data object and initialize all values to zero.
 */
TPAT_Sys_Data::TPAT_Sys_Data(){}

/**
 *	@brief Copy constructor
 *
 *	@param d a system data object
 */
TPAT_Sys_Data::TPAT_Sys_Data(const TPAT_Sys_Data& d){
	copyData(d);
}//============================================

/**
 *	@brief Destructor
 */
TPAT_Sys_Data::~TPAT_Sys_Data(){}

/**
 *	@brief Copy the system data object
 *
 *	@param d a system data object
 *	@return this data object, set to equal the input object
 */
TPAT_Sys_Data& TPAT_Sys_Data::operator =(const TPAT_Sys_Data &d){
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
bool operator ==(const TPAT_Sys_Data &lhs, const TPAT_Sys_Data &rhs){
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
bool operator !=(const TPAT_Sys_Data &lhs, const TPAT_Sys_Data &rhs){
	return ! operator==(lhs, rhs);
}//===========================================

/**
 *	@brief Copy this system data object
 *	@param d a system data reference
 */
void TPAT_Sys_Data::copyData(const TPAT_Sys_Data &d){
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
int TPAT_Sys_Data::getNumPrimaries() const { return primaries.size(); }

/**
 *	@brief Retrieve the name of one of the system primaries
 *	@param n the "index" of the primary, starts at 0
 *	@return the name of the n'th primary
 */
std::string TPAT_Sys_Data::getPrimary(int n) const{ 
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
int TPAT_Sys_Data::getPrimID(int n) const{ return primIDs.at(n); }

/**
 *	@return the characteristic length (km) associated with this system
 */
double TPAT_Sys_Data::getCharL() const { return charL; }

/**
 *	@return the charactersitic mass (kg) associated with this system
 */
double TPAT_Sys_Data::getCharM() const { return charM; }

/**
 *	@return the characteristic time (s) associated with this system
 */
double TPAT_Sys_Data::getCharT() const { return charT; }

/**
 *	@return the tpat_system_tp type associated with this system
 */
TPAT_System_Tp TPAT_Sys_Data::getType() const { return type; }

/**
 *	@return a string (human-readable) version of the system type
 */
std::string TPAT_Sys_Data::getTypeStr() const{
	switch (type){
		case TPAT_System_Tp::UNDEF_SYS:
			return "Undefined System Type";
		case TPAT_System_Tp::CR3BP_SYS:
			return "CR3BP";
		case TPAT_System_Tp::CR3BP_LTVP_SYS:
			return "CR3BP, Low Thrust, Velocity-Pointing";
		case TPAT_System_Tp::BCR4BPR_SYS:
			return "BCR4BP, Rotating Coord.";
		default:
			return "Unrecognized value... May be a bug";
	}
}//=================================================

/**
 *  @brief Save the system data object to a file
 * 
 *  @param filepath relative or absolute path to the file
 */
void TPAT_Sys_Data::saveToMat(const char *filepath) const{
	/*	Create a new Matlab MAT file with the given name and optional
	 *	header string. If no header string is given, the default string 
	 *	used containing the software, version, and date in it. If a header
	 *	string is specified, at most the first 116 characters are written to
	 *	the file. Arguments are:
	 *	const char *matname 	- 	the name of the file
	 *	const char *hdr_str 	- 	the 116 byte header string
	 *	enum mat_ft 			- 	matlab file version MAT_FT_MAT5 or MAT_FT_MAT4
	 */
	mat_t *matfp = Mat_CreateVer(filepath, NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		printErr("TPAT_Fam_CR3BP::saveToMat: Error creating MAT file\n");
	}else{
		// save things
		saveToMat(matfp);
	}

	Mat_Close(matfp);
}//======================================================
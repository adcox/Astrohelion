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

#include "tpat_sys_data.hpp"

using namespace std;

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
 *	@brief Copy the system data object
 *
 *	@param d a system data object
 *	@return this data object, set to equal the input object
 */
tpat_sys_data& tpat_sys_data::operator =(const tpat_sys_data &d){
	copyData(d);
	return *this;
}//==========================================

void tpat_sys_data::copyData(const tpat_sys_data &d){
	charL = d.charL;
	charT = d.charT;
	charM = d.charM;
	type = d.type;
	primaries = d.primaries;
	primIDs = d.primIDs;
}//==========================================

/**
 *	@return the number of primaries this system models
 */
int tpat_sys_data::getNumPrimaries() const { return numPrimaries; }

/**
 *	@brief Retrieve the name of one of the system primaries
 *	@param n the "index" of the primary, starts at 0
 *	@return the name of the n'th primary
 */
std::string tpat_sys_data::getPrimary(int n) const{ return primaries.at(n); }

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
		case tpat_sys_data::BCR4BPR_SYS:
			return "BCR4BP, Rotating Coord.";
		default:
			return "Unrecognized value... May be a bug";
	}
}//=================================================
/** 
 *	@file adtk_sys_data.cpp
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

#include "adtk_sys_data.hpp"

using namespace std;

/**
 *	Create a new system data object and initialize all values to zero.
 */
adtk_sys_data::adtk_sys_data(){
	numPrimaries = 0;
	charL = 0;
	charT = 0;
	charM = 0;
	type = UNDEF_SYS;
}

/**
 *	Copy the system data object
 *
 *	@param d a system data object
 *	@return this data object, set to equal the input object
 */
adtk_sys_data& adtk_sys_data::operator= (const adtk_sys_data &d){
	charL = d.charL;
	charT = d.charT;
	charM = d.charM;
	type = d.type;
	return *this;
}//==========================================

/**
 *	@return the characteristic length (km) associated with this system
 */
double adtk_sys_data::getCharL(){ return charL; }

/**
 *	@return the charactersitic mass (kg) associated with this system
 */
double adtk_sys_data::getCharM(){ return charM; }

/**
 *	@return the characteristic time (s) associated with this system
 */
double adtk_sys_data::getCharT(){ return charT; }

/**
 *	@return the system_t type associated with this system
 */
adtk_sys_data::system_t adtk_sys_data::getType(){ return type; }

/**
 *	@return a string (human-readable) version of the system type
 */
string adtk_sys_data::getTypeStr(){
	switch (type){
		case adtk_sys_data::UNDEF_SYS:
			return "Undefined System Type";
		case adtk_sys_data::CR3BP_SYS:
			return "CR3BP";
		case adtk_sys_data::BCR4BPR_SYS:
			return "BCR4BP, Rotating Coord.";
		default:
			return "Unrecognized value... May be a bug";
	}
}//=================================================
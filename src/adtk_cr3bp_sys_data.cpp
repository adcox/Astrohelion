/**
 *	adtk_cr3bp_sys_data.cpp
 *
 * 	System Data object specifically for CR3BP
 */

#include "adtk_cr3bp_sys_data.hpp"
#include "adtk_body_data.hpp"
#include "adtk_constants.hpp"

#include <exception>
#include <iostream>
#include <math.h>

using namespace std;

/**
 *	Default constructor
 */
adtk_cr3bp_sys_data::adtk_cr3bp_sys_data() : adtk_sys_data(){
	numPrimaries = 2;
	type = adtk_sys_data::CR3BP_SYS;
	P1 = "P1";
	P2 = "P2";
	mu = 0;
}//========================================

/**
 *	Create a system data object using data from the two primaries
 *	@param P1 the name of the larger primary
 *	@param P2 the name of the smaller primary; P2 must orbit P1
 */
adtk_cr3bp_sys_data::adtk_cr3bp_sys_data(std::string P1, std::string P2){
	numPrimaries = 2;
	type = adtk_sys_data::CR3BP_SYS;
	
	adtk_body_data p1Data(P1);
	adtk_body_data p2Data(P2);

	this->P1 = p1Data.getName();
	this->P2 = p2Data.getName();

	// Check to make sure P1 is P2's parent
	if(p2Data.getName().compare(p1Data.getName())){
		charL = p2Data.getOrbitRad();
		charM = p1Data.getMass() + p2Data.getMass();
		charT = sqrt(pow(charL, 3)/(G*charM));

		mu = p2Data.getMass()/charM;
	}else{
		cout << "adtk_cr3bp_sys_Data constructor :: P1 must be the parent of P2" << endl;
		throw;
	}
}//===================================================

/**
 *	Copy operator; makes a clean copy of a data object into this one
 *	@param d a CR3BP system data object
 *	@return this system data object
 */
adtk_cr3bp_sys_data& adtk_cr3bp_sys_data::operator= (const adtk_cr3bp_sys_data &d){
	adtk_sys_data::operator= (d);
	P1 = d.P1;
	P2 = d.P2;
	mu = d.mu;
	return *this;
}

/**
 *	@return the non-dimensional mass ratio for the system
 */
double adtk_cr3bp_sys_data::getMu(){ return mu; }

/**
 *	@param n the index of the primary (0 for P1, 1 for P2)
 *	@return the name of the primary
 */
std::string adtk_cr3bp_sys_data::getPrimary(int n){
	return n == 0 ? P1 : P2;
}
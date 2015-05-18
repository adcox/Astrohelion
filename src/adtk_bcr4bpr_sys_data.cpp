/**
 *	@file adtk_bcr4bpr_sys_data.cpp
 *
 * 	System Data object specifically for BCR4BP, rotating coordinates
 */

#include "adtk_bcr4bpr_sys_data.hpp"
 
#include "adtk_body_data.hpp"
#include "adtk_constants.hpp"

#include <exception>
#include <iostream>
#include <math.h>

using namespace std;

/**
 *	Default constructor
 */
adtk_bcr4bpr_sys_data::adtk_bcr4bpr_sys_data() : adtk_sys_data(){
	numPrimaries = 2;
	type = adtk_sys_data::BCR4BPR_SYS;
	P1 = "P1";
	P2 = "P2";
	P3 = "P3";
	mu = 0;
	nu = 0;
	charLRatio = 0;
	k = 1;
}//========================================

/**
 *	Create a system data object using data from the two primaries
 *	@param P1 the name of the larger primary
 *	@param P2 the name of the medium primary; P2 must orbit P1
 *	@param P3 the name of the smallest primary; P3 must orbit P2
 */
adtk_bcr4bpr_sys_data::adtk_bcr4bpr_sys_data(std::string P1, std::string P2, std::string P3){
	numPrimaries = 3;
	type = adtk_sys_data::BCR4BPR_SYS;
	
	adtk_body_data p1Data(P1);
	adtk_body_data p2Data(P2);
	adtk_body_data p3Data(P3);

	this->P1 = p1Data.getName();
	this->P2 = p2Data.getName();
	this->P3 = p3Data.getName();

	// Check to make sure P1 is P2's parent
	if(p2Data.getParent().compare(p1Data.getName()) == 0 && p2Data.getParent().compare(p1Data.getName()) == 0){
		k = 1.0/100.0;	// Scaling factor to make non-dimensional numbers nicer
		charL = k * p2Data.getOrbitRad();
		charM = k * (p1Data.getMass() + p2Data.getMass() + p3Data.getMass());
		charLRatio = p3Data.getOrbitRad()/charL;
		charT = sqrt(pow(charL, 3)/(G*charM));

		mu = (p2Data.getMass() + p3Data.getMass())/charM;
		nu = p3Data.getMass()/charM;
	}else{
		cout << "adtk_bcr4bpr_sys_data constructor :: P1 must be the parent of P2 and P2 must be the parent of P3" << endl;
		throw;
	}
}//===================================================

/**
 *	Copy operator; makes a clean copy of a data object into this one
 *	@param d a BCR4BPR system data object
 *	@return this system data object
 */
adtk_bcr4bpr_sys_data& adtk_bcr4bpr_sys_data::operator= (const adtk_bcr4bpr_sys_data &d){
	adtk_sys_data::operator= (d);
	P1 = d.P1;
	P2 = d.P2;
	P3 = d.P3;
	k = d.k;
	mu = d.mu;
	nu = d.nu;
	charLRatio = d.charLRatio;
	return *this;
}//=====================================

/**
 *	@return the non-dimensional mass ratio for the secondary system (P2 + P3)
 */
double adtk_bcr4bpr_sys_data::getMu(){ return mu; }

/**
 *	@return the non-dimensional mass ratio for P3
 */
double adtk_bcr4bpr_sys_data::getNu(){ return nu; }

/**
 *	@return the ratio between P3's orbital radius and P2's orbital radius, non-dimensional units
 */
double adtk_bcr4bpr_sys_data::getCharLRatio(){ return charLRatio; }

/**
 *	@return the scaling constant for this system
 */
double adtk_bcr4bpr_sys_data::getK(){ return k; }

/**
 *	@param n the index of the primary (0 for P1, 1 for P2, 2 for P3)
 *	@return the name of the primary
 */
std::string adtk_bcr4bpr_sys_data::getPrimary(int n){
	switch(n){
		case 0:
			return P1;
		case 1:
			return P2;
		case 2:
			return P3;
		default:
			cout << "adtk_bcr4bpr_sys_data :: Cannot get name of primary #" << n << endl;
			throw;
	}
}
/**
 *	@file tpat_bcr4bpr_sys_data.cpp
 *
 * 	System Data object specifically for BCR4BP, rotating coordinates
 */

#include "tpat_bcr4bpr_sys_data.hpp"
 
#include "tpat_body_data.hpp"
#include "tpat_constants.hpp"
#include "tpat_utilities.hpp"
 
#include <cmath>
#include <exception>
#include <iostream>

using namespace std;

/**
 *	@brief Default constructor
 */
tpat_bcr4bpr_sys_data::tpat_bcr4bpr_sys_data() : tpat_sys_data(){
	numPrimaries = 3;
	type = tpat_sys_data::BCR4BPR_SYS;
}//========================================

/**
 *	@brief Create a system data object using data from the two primaries
 *	@param P1 the name of the larger primary
 *	@param P2 the name of the medium primary; P2 must orbit P1
 *	@param P3 the name of the smallest primary; P3 must orbit P2
 */
tpat_bcr4bpr_sys_data::tpat_bcr4bpr_sys_data(std::string P1, std::string P2, std::string P3){
	numPrimaries = 3;
	type = tpat_sys_data::BCR4BPR_SYS;
	
	tpat_body_data p1Data(P1);
	tpat_body_data p2Data(P2);
	tpat_body_data p3Data(P3);

	primaries.push_back(p1Data.getName());
	primaries.push_back(p2Data.getName());
	primaries.push_back(p3Data.getName());

	primIDs.push_back(p1Data.getID());
	primIDs.push_back(p2Data.getID());
	primIDs.push_back(p3Data.getID());

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
		cout << "tpat_bcr4bpr_sys_data constructor :: P1 must be the parent of P2 and P2 must be the parent of P3" << endl;
		throw;
	}
}//===================================================

/**
 *	@brief Copy constructor
 *	@param d
 */
tpat_bcr4bpr_sys_data::tpat_bcr4bpr_sys_data(const tpat_bcr4bpr_sys_data &d) : tpat_sys_data(d){
	copyData(d);
}//==============================================

/**
 *	@brief Copy data from a system data object into this one
 *	@param d
 */
void tpat_bcr4bpr_sys_data::copyData(const tpat_bcr4bpr_sys_data &d){
	k = d.k;
	mu = d.mu;
	nu = d.nu;
	charLRatio = d.charLRatio;
}//==============================================

/**
 *	@brief Copy operator; makes a clean copy of a data object into this one
 *	@param d a BCR4BPR system data object
 *	@return this system data object
 */
tpat_bcr4bpr_sys_data& tpat_bcr4bpr_sys_data::operator= (const tpat_bcr4bpr_sys_data &d){
	tpat_sys_data::operator= (d);
	copyData(d);
	return *this;
}//=====================================

/**
 *	@return the non-dimensional mass ratio for the secondary system (P2 + P3)
 */
double tpat_bcr4bpr_sys_data::getMu() const { return mu; }

/**
 *	@return the non-dimensional mass ratio for P3
 */
double tpat_bcr4bpr_sys_data::getNu() const { return nu; }

/**
 *	@return the ratio between P3's orbital radius and P2's orbital radius, non-dimensional units
 */
double tpat_bcr4bpr_sys_data::getCharLRatio() const { return charLRatio; }

/**
 *	@return the scaling constant for this system
 */
double tpat_bcr4bpr_sys_data::getK() const { return k; }

/**
 *	@return the angle between the P1/P2 line and the inertial x-axis at time t = 0,
 *	angle in radians
 */
double tpat_bcr4bpr_sys_data::getTheta0() const { return theta0; }

/**
 *	@return the angle between the P2/P3 line (projected onto inertial XY plane) and the 
 *	x-axis at time t = 0, angle in radians
 */
double tpat_bcr4bpr_sys_data::getPhi0() const { return phi0; }

/**
 *	@return the inclination of the P2/P3 orbital plane relative to the P1/P2 orbital plane, radians.
 */
double tpat_bcr4bpr_sys_data::getGamma() const { return gamma; }

/**
 *	@brief Set the angle theta0
 *	@param t angle in radians
 */
void tpat_bcr4bpr_sys_data::setTheta0(double t){ theta0 = t; }

/**
 *	@brief Set the angle phi0
 *	@param t angle in radians
 */
void tpat_bcr4bpr_sys_data::setPhi0(double t){ phi0 = t; }

/**
 *	@brief Set the angle gamma
 *	@param t angle in radians
 */
void tpat_bcr4bpr_sys_data::setGamma(double t){ gamma = t; }

/**
 *	@brief Save system data, like the names of the primaries and the system mass ratio, to a .mat file
 *	@param matFile a pointer to the .mat file
 */
void tpat_bcr4bpr_sys_data::saveToMat(mat_t *matFile){
	size_t dims[2] = {1,1};

	// Initialize character array (larger than needed), copy in the name of the primary, then create a var.
	char p1_str[64];
	strcpy(p1_str, primaries.at(0).c_str());
	dims[1] = primaries.at(0).length();
	matvar_t *p1_var = Mat_VarCreate("P1", MAT_C_CHAR, MAT_T_UTF8, 2, dims, p1_str, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p1_var, "P1", MAT_COMPRESSION_NONE);

	char p2_str[64];
	strcpy(p2_str, primaries.at(1).c_str());
	dims[1] = primaries.at(1).length();
	matvar_t *p2_var = Mat_VarCreate("P2", MAT_C_CHAR, MAT_T_UTF8, 2, dims, &(p2_str[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p2_var, "P2", MAT_COMPRESSION_NONE);

	char p3_str[64];
	strcpy(p3_str, primaries.at(2).c_str());
	dims[1] = primaries.at(2).length();
	matvar_t *p3_var = Mat_VarCreate("P3", MAT_C_CHAR, MAT_T_UTF8, 2, dims, &(p3_str[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p3_var, "P3", MAT_COMPRESSION_NONE);

	dims[1] = 1;
	matvar_t *theta0_var = Mat_VarCreate("Theta0", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &theta0, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, theta0_var, "Theta0", MAT_COMPRESSION_NONE);

	matvar_t *phi0_var = Mat_VarCreate("Phi0", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &phi0, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, phi0_var, "Phi0", MAT_COMPRESSION_NONE);

	matvar_t *gamma_var = Mat_VarCreate("Gamma", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &gamma, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, gamma_var, "Gamma", MAT_COMPRESSION_NONE);
}//===================================================





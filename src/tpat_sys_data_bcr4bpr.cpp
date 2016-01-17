/**
 *	@file tpat_sys_data_bcr4bpr.cpp
 *	@brief Derivative of tpat_sys_data, specific to BCR4BPR
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

#include "tpat_sys_data_bcr4bpr.hpp"
 
#include "tpat_body_data.hpp"
#include "tpat_constants.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_utilities.hpp"
 
#include <cmath>
#include <cstring>
#include <exception>

// Static variable initialization
double tpat_sys_data_bcr4bpr::REF_EPOCH = 172650160;	// 2005/06/21 18:21:35

/**
 *	@brief Default constructor
 */
tpat_sys_data_bcr4bpr::tpat_sys_data_bcr4bpr() : tpat_sys_data(){
	numPrimaries = 3;
	type = tpat_sys_data::BCR4BPR_SYS;
	otherParams.assign(7,0);
}//========================================

/**
 *	@brief Create a system data object using data from the two primaries
 *	@param P1 the name of the larger primary
 *	@param P2 the name of the medium primary; P2 must orbit P1
 *	@param P3 the name of the smallest primary; P3 must orbit P2
 */
tpat_sys_data_bcr4bpr::tpat_sys_data_bcr4bpr(std::string P1, std::string P2, std::string P3){
	numPrimaries = 3;
	type = tpat_sys_data::BCR4BPR_SYS;
	otherParams.assign(7,0);

	initFromPrimNames(P1, P2, P3);
}//===================================================

/**
 *	@brief Initialize all data fiels using the primaries' names
 *	@param P1 name of the largest primary in the entire system
 *	@param P2 name of the larger primary in the secondary system
 *	@param P3 name of the smaller primary int he secondary system
 */
void tpat_sys_data_bcr4bpr::initFromPrimNames(std::string P1, std::string P2, std::string P3){
	tpat_body_data p1Data(P1);
	tpat_body_data p2Data(P2);
	tpat_body_data p3Data(P3);

	primaries.clear();
	primaries.push_back(p1Data.getName());
	primaries.push_back(p2Data.getName());
	primaries.push_back(p3Data.getName());

	primIDs.clear();
	primIDs.push_back(p1Data.getID());
	primIDs.push_back(p2Data.getID());
	primIDs.push_back(p3Data.getID());

	// Check to make sure P1 is P2's parent
	if(p2Data.getParent().compare(p1Data.getName()) == 0 && p2Data.getParent().compare(p1Data.getName()) == 0){
		double k = 1.0/100.0;	// Scaling factor to make non-dimensional numbers nicer
		// double k = 1.0;
		
		charL = k * p2Data.getOrbitRad();
		charM = k * (p1Data.getMass() + p2Data.getMass() + p3Data.getMass());
		double charLRatio = p3Data.getOrbitRad()/charL;
		charT = sqrt(pow(charL, 3)/(G*charM));

		double mu = (p2Data.getMass() + p3Data.getMass())/charM;
		double nu = p3Data.getMass()/charM;
		double theta = 0;
		double phi = 0;
		double gamma = 5.14*PI/180;

		otherParams[0] = mu;
		otherParams[1] = nu;
		otherParams[2] = k;
		otherParams[3] = charLRatio;
		otherParams[4] = theta;
		otherParams[5] = phi;
		otherParams[6] = gamma;
	}else{
		throw tpat_exception("P1 must be the parent of P2 and P2 must be the parent of P3");
	}
}//===================================================

/**
 *	@brief Copy constructor
 *	@param d
 */
tpat_sys_data_bcr4bpr::tpat_sys_data_bcr4bpr(const tpat_sys_data_bcr4bpr &d) : tpat_sys_data(d){}

/**
 *	@brief Copy operator; makes a clean copy of a data object into this one
 *	@param d a BCR4BPR system data object
 *	@return this system data object
 */
tpat_sys_data_bcr4bpr& tpat_sys_data_bcr4bpr::operator= (const tpat_sys_data_bcr4bpr &d){
	tpat_sys_data::operator= (d);
	return *this;
}//=====================================

/**
 *	@brief Retrieve the model that governs the motion for this system type
 *	@return the model that governs the motion for this system type
 */
const tpat_model* tpat_sys_data_bcr4bpr::getModel() const { return &model; }

/**
 *	@return the non-dimensional mass ratio for the secondary system (P2 + P3)
 */
double tpat_sys_data_bcr4bpr::getMu() const { return otherParams.at(0); }

/**
 *	@return the non-dimensional mass ratio for P3
 */
double tpat_sys_data_bcr4bpr::getNu() const { return otherParams.at(1); }

/**
 *	@return the ratio between P3's orbital radius and P2's orbital radius, non-dimensional units
 */
double tpat_sys_data_bcr4bpr::getCharLRatio() const { return otherParams.at(3); }

/**
 *	@return the scaling constant for this system
 */
double tpat_sys_data_bcr4bpr::getK() const {
	if(otherParams.size() != 7)
		printf("");

	return otherParams.at(2); }

/**
 *	@return the angle between the P1/P2 line and the inertial x-axis at time t = 0,
 *	angle in radians
 */
double tpat_sys_data_bcr4bpr::getTheta0() const { return otherParams.at(4); }

/**
 *	@return the angle between the P2/P3 line (projected onto inertial XY plane) and the 
 *	x-axis at time t = 0, angle in radians
 */
double tpat_sys_data_bcr4bpr::getPhi0() const { return otherParams.at(5); }

/**
 *	@return the inclination of the P2/P3 orbital plane relative to the P1/P2 orbital plane, radians.
 */
double tpat_sys_data_bcr4bpr::getGamma() const { return otherParams.at(6); }

/**
 *	@brief Set the angle theta0
 *	@param t angle in radians
 */
void tpat_sys_data_bcr4bpr::setTheta0(double t){ otherParams.at(4) = t; }

/**
 *	@brief Set the angle phi0
 *	@param t angle in radians
 */
void tpat_sys_data_bcr4bpr::setPhi0(double t){ otherParams.at(5) = t; }

/**
 *	@brief Set the angle gamma
 *	@param t angle in radians
 */
void tpat_sys_data_bcr4bpr::setGamma(double t){ otherParams.at(6) = t; }

/**
 *	@brief Save system data, like the names of the primaries and the system mass ratio, to a .mat file
 *	@param matFile a pointer to the .mat file
 */
void tpat_sys_data_bcr4bpr::saveToMat(mat_t *matFile) const{
	size_t dims[2] = {1,1};

	// Initialize character array (larger than needed), copy in the name of the primary, then create a var.
	char p1_str[64];
	strcpy(p1_str, primaries.at(0).c_str());
	dims[1] = primaries.at(0).length();
	matvar_t *p1_var = Mat_VarCreate("P1", MAT_C_CHAR, MAT_T_UINT8, 2, dims, p1_str, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p1_var, "P1", MAT_COMPRESSION_NONE);

	char p2_str[64];
	strcpy(p2_str, primaries.at(1).c_str());
	dims[1] = primaries.at(1).length();
	matvar_t *p2_var = Mat_VarCreate("P2", MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(p2_str[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p2_var, "P2", MAT_COMPRESSION_NONE);

	char p3_str[64];
	strcpy(p3_str, primaries.at(2).c_str());
	dims[1] = primaries.at(2).length();
	matvar_t *p3_var = Mat_VarCreate("P3", MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(p3_str[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p3_var, "P3", MAT_COMPRESSION_NONE);

	dims[1] = 1;
	double mu = otherParams[0];
	matvar_t *mu_var = Mat_VarCreate("Mu", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &mu, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, mu_var, "Mu", MAT_COMPRESSION_NONE);

	double nu = otherParams[1];
	matvar_t *nu_var = Mat_VarCreate("Nu", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &nu, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, nu_var, "Nu", MAT_COMPRESSION_NONE);

	double k = otherParams[2];
	matvar_t *k_var = Mat_VarCreate("K", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &k, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, k_var, "K", MAT_COMPRESSION_NONE);

	double r = otherParams[3];
	matvar_t *ratio_var = Mat_VarCreate("CharLRatio", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &r, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, ratio_var, "CharLRatio", MAT_COMPRESSION_NONE);

	double t = otherParams[4];
	matvar_t *theta0_var = Mat_VarCreate("Theta0", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &t, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, theta0_var, "Theta0", MAT_COMPRESSION_NONE);

	double p = otherParams[5];
	matvar_t *phi0_var = Mat_VarCreate("Phi0", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &p, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, phi0_var, "Phi0", MAT_COMPRESSION_NONE);

	double g = otherParams[6];
	matvar_t *gamma_var = Mat_VarCreate("Gamma", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &g, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, gamma_var, "Gamma", MAT_COMPRESSION_NONE);


}//===================================================

/**
 *	@brief Load system data from a mat file
 *	@param matFile a pointer to the mat file in question
 */
void tpat_sys_data_bcr4bpr::readFromMat(mat_t *matFile){
	std::string P1 = readStringFromMat(matFile, "P1", MAT_T_UINT8, MAT_C_CHAR);
	std::string P2 = readStringFromMat(matFile, "P2", MAT_T_UINT8, MAT_C_CHAR);
	std::string P3 = readStringFromMat(matFile, "P3", MAT_T_UINT8, MAT_C_CHAR);
	initFromPrimNames(P1, P2, P3);

	otherParams[4] = readDoubleFromMat(matFile, "Theta0");
	otherParams[5] = readDoubleFromMat(matFile, "Phi0");
	otherParams[6] = readDoubleFromMat(matFile, "Gamma");
}//====================================================





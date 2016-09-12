/**
 *	@file SysData_bc4bp.cpp
 *	@brief Derivative of SysData, specific to BCR4BPR
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2016, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of Astrohelion
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

#include "SysData_bc4bp.hpp"
 
#include "BodyData.hpp"
#include "Common.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"
 
#include <cmath>
#include <cstring>
#include <exception>

namespace astrohelion{

// Static variable initialization
double SysData_bc4bp::REF_EPOCH = 172650160;	// 2005/06/21 18:21:35

/**
 *	@brief Default constructor
 */
SysData_bc4bp::SysData_bc4bp() : SysData(){
	numPrimaries = 3;
	type = SysData_tp::BCR4BPR_SYS;
	otherParams.assign(8,0);
}//========================================

/**
 *	@brief Create a system data object using data from the two primaries
 *	@param P1 the name of the larger primary
 *	@param P2 the name of the medium primary; P2 must orbit P1
 *	@param P3 the name of the smallest primary; P3 must orbit P2
 */
SysData_bc4bp::SysData_bc4bp(std::string P1, std::string P2, std::string P3){
	numPrimaries = 3;
	type = SysData_tp::BCR4BPR_SYS;
	otherParams.assign(8,0);

	initFromPrimNames(P1, P2, P3);
}//===================================================

/**
 *  @brief Load the system data object from a Matlab data file
 * 
 *  @param filepath path to the data file
 *  @throws Exception if the data file cannot be loaded
 */
SysData_bc4bp::SysData_bc4bp(const char *filepath){
	numPrimaries = 3;
	type = SysData_tp::BCR4BPR_SYS;
	otherParams.assign(8,0);

	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw Exception("SysData_bc4bp: Could not open data file");
	}
	readFromMat(matfp);
	Mat_Close(matfp);
}//===================================================

/**
 *	@brief Initialize all data fiels using the primaries' names
 *	@param P1 name of the largest primary in the entire system
 *	@param P2 name of the larger primary in the secondary system
 *	@param P3 name of the smaller primary int he secondary system
 *	@throws Exception if the system architecture is incorrect, i.e.
 *	P1 must be the parent of P2 and P2 must be the parent of P3
 */
void SysData_bc4bp::initFromPrimNames(std::string P1, std::string P2, std::string P3){
	BodyData p1Data(P1);
	BodyData p2Data(P2);
	BodyData p3Data(P3);

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
		otherParams[7] = REF_EPOCH;
	}else{
		throw Exception("P1 must be the parent of P2 and P2 must be the parent of P3");
	}
}//===================================================

/**
 *	@brief Copy constructor
 *	@param d
 */
SysData_bc4bp::SysData_bc4bp(const SysData_bc4bp &d) : SysData(d){}

/**
 *	@brief Copy operator; makes a clean copy of a data object into this one
 *	@param d a BCR4BPR system data object
 *	@return this system data object
 */
SysData_bc4bp& SysData_bc4bp::operator= (const SysData_bc4bp &d){
	SysData::operator= (d);
	return *this;
}//=====================================

/**
 *	@brief Retrieve the model that governs the motion for this system type
 *	@return the model that governs the motion for this system type
 */
const DynamicsModel* SysData_bc4bp::getDynamicsModel() const { return &model; }

/**
 *	@return the non-dimensional mass ratio for the secondary system (P2 + P3)
 */
double SysData_bc4bp::getMu() const { return otherParams.at(0); }

/**
 *	@return the non-dimensional mass ratio for P3
 */
double SysData_bc4bp::getNu() const { return otherParams.at(1); }

/**
 *	@return the ratio between P3's orbital radius and P2's orbital radius, non-dimensional units
 */
double SysData_bc4bp::getCharLRatio() const { return otherParams.at(3); }

/**
 *	@return the scaling constant for this system
 */
double SysData_bc4bp::getK() const {
	return otherParams.at(2);
}//====================================================

/**
 *  @brief Retrieve the epoch associated with T = 0 (seconds, J2000, UTC)
 *  @details The angles Theta0 and Phi0 coincide with time T = 0, but this is
 *  NOT the J2000 epoch T = 0; Rather, T = 0 is the epoch relative to some reference
 *  epoch T0, which is returned by this function.
 *  @return the epoch associated with T = 0 (seconds, J2000, UTC)
 */
double SysData_bc4bp::getEpoch0() const { return otherParams.at(7); }

/**
 *	@return the angle between the P1/P2 line and the inertial x-axis at time t = 0,
 *	angle in radians
 */
double SysData_bc4bp::getTheta0() const { return otherParams.at(4); }

/**
 *	@return the angle between the P2/P3 line (projected onto inertial XY plane) and the 
 *	x-axis at time t = 0, angle in radians
 */
double SysData_bc4bp::getPhi0() const { return otherParams.at(5); }

/**
 *	@return the inclination of the P2/P3 orbital plane relative to the P1/P2 orbital plane, radians.
 */
double SysData_bc4bp::getGamma() const { return otherParams.at(6); }

/**
 *  @brief Set the reference epoch associated with T = 0 for this system (seconds, J2000, UTC)
 *  @details The angles Theta0 and Phi0 coincide with time T = 0, but this is
 *  NOT the J2000 epoch T = 0; Rather, T = 0 is the epoch relative to some reference
 *  epoch T0, which is set by this function.
 * 
 *  @param T the reference epoch associated with T = 0 for this system (seconds, J2000, UTC)
 */
void SysData_bc4bp::setEpoch0(double T){ otherParams.at(7) = T; }

/**
 *	@brief Set the angle theta0
 *	@param t angle in radians
 */
void SysData_bc4bp::setTheta0(double t){ otherParams.at(4) = t; }

/**
 *	@brief Set the angle phi0
 *	@param t angle in radians
 */
void SysData_bc4bp::setPhi0(double t){ otherParams.at(5) = t; }

/**
 *	@brief Set the angle gamma
 *	@param t angle in radians
 */
void SysData_bc4bp::setGamma(double t){ otherParams.at(6) = t; }

/**
 *  @brief Save the system data to a matlab file
 * 
 *  @param filepath path to the data file
 */
void SysData_bc4bp::saveToMat(const char *filepath) const{
	SysData::saveToMat(filepath);
}//==================================================

/**
 *	@brief Save system data, like the names of the primaries and the system mass ratio, to a .mat file
 *	@param matFile a pointer to the .mat file
 */
void SysData_bc4bp::saveToMat(mat_t *matFile) const{
	size_t dims[2] = {1,1};

	// Initialize character array (larger than needed), copy in the name of the primary, then create a var.
	char p1_str[64];
	strcpy(p1_str, primaries.at(0).c_str());
	dims[1] = primaries.at(0).length();
	matvar_t *p1_var = Mat_VarCreate("P1", MAT_C_CHAR, MAT_T_UINT8, 2, dims, p1_str, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, p1_var, "P1", MAT_COMPRESSION_NONE);

	char p2_str[64];
	strcpy(p2_str, primaries.at(1).c_str());
	dims[1] = primaries.at(1).length();
	matvar_t *p2_var = Mat_VarCreate("P2", MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(p2_str[0]), MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, p2_var, "P2", MAT_COMPRESSION_NONE);

	char p3_str[64];
	strcpy(p3_str, primaries.at(2).c_str());
	dims[1] = primaries.at(2).length();
	matvar_t *p3_var = Mat_VarCreate("P3", MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(p3_str[0]), MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, p3_var, "P3", MAT_COMPRESSION_NONE);

	dims[1] = 1;
	double mu = otherParams[0];
	matvar_t *mu_var = Mat_VarCreate("Mu", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &mu, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, mu_var, "Mu", MAT_COMPRESSION_NONE);

	double nu = otherParams[1];
	matvar_t *nu_var = Mat_VarCreate("Nu", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &nu, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, nu_var, "Nu", MAT_COMPRESSION_NONE);

	double k = otherParams[2];
	matvar_t *k_var = Mat_VarCreate("K", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &k, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, k_var, "K", MAT_COMPRESSION_NONE);

	double r = otherParams[3];
	matvar_t *ratio_var = Mat_VarCreate("CharLRatio", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &r, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, ratio_var, "CharLRatio", MAT_COMPRESSION_NONE);

	double t = otherParams[4];
	matvar_t *theta0_var = Mat_VarCreate("Theta0", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &t, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, theta0_var, "Theta0", MAT_COMPRESSION_NONE);

	double p = otherParams[5];
	matvar_t *phi0_var = Mat_VarCreate("Phi0", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &p, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, phi0_var, "Phi0", MAT_COMPRESSION_NONE);

	double g = otherParams[6];
	matvar_t *gamma_var = Mat_VarCreate("Gamma", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &g, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, gamma_var, "Gamma", MAT_COMPRESSION_NONE);

	double T = otherParams[7];
	matvar_t *epoch_var = Mat_VarCreate("Epoch0", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &T, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, epoch_var, "Epoch0", MAT_COMPRESSION_NONE);
}//===================================================

/**
 *	@brief Load system data from a mat file
 *	@param matFile a pointer to the mat file in question
 */
void SysData_bc4bp::readFromMat(mat_t *matFile){
	std::string P1 = astrohelion::readStringFromMat(matFile, "P1", MAT_T_UINT8, MAT_C_CHAR);
	std::string P2 = astrohelion::readStringFromMat(matFile, "P2", MAT_T_UINT8, MAT_C_CHAR);
	std::string P3 = astrohelion::readStringFromMat(matFile, "P3", MAT_T_UINT8, MAT_C_CHAR);
	
	numPrimaries = 3;
	type = SysData_tp::BCR4BPR_SYS;
	otherParams.assign(8,0);
	initFromPrimNames(P1, P2, P3);

	otherParams[4] = astrohelion::readDoubleFromMat(matFile, "Theta0");
	otherParams[5] = astrohelion::readDoubleFromMat(matFile, "Phi0");
	otherParams[6] = astrohelion::readDoubleFromMat(matFile, "Gamma");
	otherParams[7] = astrohelion::readDoubleFromMat(matFile, "Epoch0");
}//====================================================




}// END of Astrohelion namespace
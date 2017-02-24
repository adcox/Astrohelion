/**
 *	@file SysData_2bp.cpp
 *	@brief Derivative of SysData, specific to 2BP
 *	
 *	@author Andrew Cox
 *	@version August 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "SysData_2bp.hpp"
 
#include "BodyData.hpp"
#include "Common.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"

#include <cstring>
#include <cmath>
#include <exception>

namespace astrohelion{
/**
 *	@brief Default constructor
 */
SysData_2bp::SysData_2bp() : SysData(){
	numPrimaries = 1;
	type = SysData_tp::R2BP_SYS;
	otherParams.assign(1,0);	// make mu = 0
}//========================================

/**
 *	@brief Create a system data object using data from the primary
 *	@param P1 the name of the primary
 */
SysData_2bp::SysData_2bp(std::string P1){
	numPrimaries = 1;
	type = SysData_tp::R2BP_SYS;
	otherParams.assign(1,0);	// make mu = 0
	
	initFromPrimNames(P1);
}//===================================================

/**
 *  @brief Load the system data object from a Matlab data file
 * 
 *  @param filepath path to the data file
 */
void SysData_2bp::initFromFile(const char *filepath){
	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw Exception("SysData_2bp: Could not load data from file");
	}
	readFromMat(matfp);
	Mat_Close(matfp);
}//===================================================

/**
 *  @brief Destructor
 */
SysData_2bp::~SysData_2bp(){}

/**
 *	@brief Initialize all data fields using the name of the primary
 *	@param P1 the name of the primary
 */
void SysData_2bp::initFromPrimNames(std::string P1){
	BodyData p1Data(P1);

	primaries.clear();
	primIDs.clear();

	primaries.push_back(p1Data.getName());
	primIDs.push_back(p1Data.getID());

	// Characteristic quantities are all unity; system uses dimensional quantities
	charL = 1;
	charM = 1;
	charT = 1;

	otherParams.at(0) = G*p1Data.getMass();	// Mu for the primary body
}//===================================================

/**
 *	@brief Retrieve the model that governs the motion for this system type
 *	@return the model that governs the motion for this system type
 */
const DynamicsModel* SysData_2bp::getDynamicsModel() const { return &model; }

/**
 *	@brief Copy constructor
 *	@param d
 */
SysData_2bp::SysData_2bp(const SysData_2bp &d) : SysData(d){}

/**
 *	@brief Copy operator; makes a clean copy of a data object into this one
 *	@param d a CR3BP system data object
 *	@return this system data object
 */
SysData_2bp& SysData_2bp::operator= (const SysData_2bp &d){
	SysData::operator= (d);
	return *this;
}//===================================================

/**
 *	@return the mass parameter (km^3/s^2)
 */
double SysData_2bp::getMu() const { return otherParams.at(0); }

/**
 *  @brief Save the system data to a matlab file
 * 
 *  @param filepath path to the data file
 */
void SysData_2bp::saveToMat(const char *filepath) const{
	SysData::saveToMat(filepath);
}//==================================================

/**
 *	@brief Save system data, like the names of the primaries and the system mass ratio, to a .mat file
 *	@param matFile a pointer to the .mat file
 */
void SysData_2bp::saveToMat(mat_t *matFile) const{
	size_t dims[2] = {1,1};

	if(primaries.size() < 1){
		astrohelion::printErr("Primaries size is %zu\n", primaries.size());
		throw Exception("SysData_2bp::saveToMat: There are no primaries?");
	}

	// Initialize character array (larger than needed), copy in the name of the primary, then create a var.
	char p1_str[64];
	strcpy(p1_str, primaries.at(0).c_str());
	dims[1] = primaries.at(0).length();
	matvar_t *p1_var = Mat_VarCreate("P1", MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(p1_str[0]), MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, p1_var, "P1", MAT_COMPRESSION_NONE);

	dims[1] = 1;
	double mu = otherParams[0];
	matvar_t *mu_var = Mat_VarCreate("Mu", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &mu, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, mu_var, "Mu", MAT_COMPRESSION_NONE);
}//===================================================

/**
 *	@brief Populate data fiels for this data object by reading the primaries'
 *	names from a Mat file
 *	@param matFile a pointer to the Mat file in question
 */
void SysData_2bp::readFromMat(mat_t *matFile){
	std::string P1 = astrohelion::readStringFromMat(matFile, "P1", MAT_T_UINT8, MAT_C_CHAR);

	numPrimaries = 1;
	type = SysData_tp::R2BP_SYS;
	otherParams.assign(1,0);	// make mu = 0
	initFromPrimNames(P1);
}//===================================================


}// END of Astrohelion namespace
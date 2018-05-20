/**
 *	@file SysData_cr3bp_lt.cpp
 *	@brief Derivative of SysData, specific to CR3BP-LTVP
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "SysData_cr3bp_lt.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"

#include <cstring>


namespace astrohelion{

//------------------------------------------------------------------------------------------------------
//      Constructors
//------------------------------------------------------------------------------------------------------

/**
 *	@brief Default constructor
 */
SysData_cr3bp_lt::SysData_cr3bp_lt() : SysData_cr3bp(){
	numPrimaries = 2;
	type = SysData_tp::CR3BP_LT_SYS;
	otherParams.assign(2,0);	// make mu = 0 and M0 = 0
	otherParams[1] = 1;		// Assign M0 = 1 to avoid singularities
}//========================================

/**
 *	@brief Create a system data object using data from the two primaries
 *	@param P1 the name of the larger primary
 *	@param P2 the name of the smaller primary; P2 must orbit P1
 *	@param refMass reference mass, kilograms
 */
SysData_cr3bp_lt::SysData_cr3bp_lt(std::string P1, std::string P2, double refMass){
	numPrimaries = 2;
	type = SysData_tp::CR3BP_LT_SYS;
	otherParams.assign(2,0);	// make mu = 0 and M0 = 0
	
	initFromPrimNames(P1, P2);

	if(refMass <= 0)
		throw Exception("SysData_cr3bp_lt: Cannot use negative or zero reference mass");
	
	otherParams[1] = refMass;
}//===================================================

/**
 *  @brief Load the system data object from a Matlab data file
 * 
 *  @param filepath path to the data file
 */
SysData_cr3bp_lt::SysData_cr3bp_lt(const char *filepath){
	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(nullptr == matfp){
		char msg[128];
		sprintf(msg, "SysData_cr3bp_lt: Could not load data from %s", 
			filepath);
		throw Exception(msg);
	}
	readFromMat(matfp);
	Mat_Close(matfp);
}//===================================================

/**
 *	@brief Copy constructor
 *	@param d 
 */
SysData_cr3bp_lt::SysData_cr3bp_lt(const SysData_cr3bp_lt &d) : SysData_cr3bp(d){}

//------------------------------------------------------------------------------------------------------
//      Operators
//------------------------------------------------------------------------------------------------------

/**
 *	@brief Copy operator; makes a clean copy of a data object into this one
 *	@param d a CR3BP system data object
 *	@return this system data object
 */
SysData_cr3bp_lt& SysData_cr3bp_lt::operator =(const SysData_cr3bp_lt &d){
	SysData_cr3bp::operator= (d);
	return *this;
}//===================================================

//------------------------------------------------------------------------------------------------------
//      Set and Get Functions
//------------------------------------------------------------------------------------------------------

/**
 *	@brief Retrieve the model that governs the motion for this system type
 *	@return the model that governs the motion for this system type
 */
const DynamicsModel* SysData_cr3bp_lt::getDynamicsModel() const { return &model; }

/**
 *	@brief Get the reference mass for the spacecraft
 *	@return the reference mass of the spacecraft, kg
 */
double SysData_cr3bp_lt::getRefMass() const { return otherParams[1]; }

/**
 *	@brief Set the reference mass for the spacecraft
 *	@param mass the mass, kg
 */
void SysData_cr3bp_lt::setRefMass(double mass){ otherParams[1] = mass; }

//------------------------------------------------------------------------------------------------------
//      Utility Functions
//------------------------------------------------------------------------------------------------------

/**
 *  @brief Save the system data to a matlab file
 * 
 *  @param filepath path to the data file
 */
void SysData_cr3bp_lt::saveToMat(const char *filepath) const{
	SysData::saveToMat(filepath);
}//==================================================

/**
 *	@brief Save system data, like the names of the primaries and the system mass ratio, to a .mat file
 *	@param matFile a pointer to the .mat file
 *	@throws Exception if the number of primaries is incorrect
 */
void SysData_cr3bp_lt::saveToMat(mat_t *matFile) const{
	size_t dims[2] = {1,1};

	if(primaries.size() < 2){
		astrohelion::printErr("Primaries size is %zu\n", primaries.size());
		throw Exception("SysData_cr3bp::saveToMat: There are no primaries?");
	}

	// Initialize character array (larger than needed), copy in the name of the primary, then create a var.
	char p1_str[64];
	std::strcpy(p1_str, primaries.at(0).c_str());
	dims[1] = primaries.at(0).length();
	matvar_t *p1_var = Mat_VarCreate("P1", MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(p1_str[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p1_var, "P1", MAT_COMPRESSION_NONE);

	char p2_str[64];
	std::strcpy(p2_str, primaries.at(1).c_str());
	dims[1] = primaries.at(1).length();
	matvar_t *p2_var = Mat_VarCreate("P2", MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(p2_str[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p2_var, "P2", MAT_COMPRESSION_NONE);

	dims[1] = 1;	
	double mu = otherParams[0];
	matvar_t *mu_var = Mat_VarCreate("Mu", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &mu, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, mu_var, "Mu", MAT_COMPRESSION_NONE);

	// Save mass in dimensional units
	double m = getRefMass();
	matvar_t *mass_var = Mat_VarCreate("Mass0", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &m, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, mass_var, "Mass0", MAT_COMPRESSION_NONE);
}//===================================================

/**
 *	@brief Populate data fiels for this data object by reading the primaries'
 *	names from a Mat file
 *	@param matFile a pointer to the Mat file in question
 */
void SysData_cr3bp_lt::readFromMat(mat_t *matFile){
	std::string P1 = astrohelion::readStringFromMat(matFile, "P1", MAT_T_UINT8, MAT_C_CHAR);
	std::string P2 = astrohelion::readStringFromMat(matFile, "P2", MAT_T_UINT8, MAT_C_CHAR);

	type = SysData_tp::CR3BP_LT_SYS;
	otherParams.assign(4,0);
	initFromPrimNames(P1, P2);

	setRefMass(readDoubleFromMat(matFile, "Mass0"));
}//===================================================




}// END of Astrohelion namespace
/**
 *	@file tpat_sys_data_cr3bp.cpp
 *	@brief Derivative of TPAT_Sys_Data, specific to CR3BP
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

#include "tpat_sys_data_cr3bp.hpp"
 
#include "tpat_body_data.hpp"
#include "tpat_constants.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_utilities.hpp"

#include <cstring>
#include <cmath>
#include <exception>

/**
 *	@brief Default constructor
 */
TPAT_Sys_Data_CR3BP::TPAT_Sys_Data_CR3BP() : TPAT_Sys_Data(){
	numPrimaries = 2;
	type = TPAT_System_Tp::CR3BP_SYS;
	otherParams.assign(1,0);	// make mu = 0
}//========================================

/**
 *	@brief Create a system data object using data from the two primaries
 *	@param P1 the name of the larger primary
 *	@param P2 the name of the smaller primary; P2 must orbit P1
 */
TPAT_Sys_Data_CR3BP::TPAT_Sys_Data_CR3BP(std::string P1, std::string P2){
	numPrimaries = 2;
	type = TPAT_System_Tp::CR3BP_SYS;
	otherParams.assign(1,0);	// make mu = 0
	
	initFromPrimNames(P1, P2);
}//===================================================

/**
 *  @brief Load the system data object from a Matlab data file
 * 
 *  @param filepath path to the data file
 */
TPAT_Sys_Data_CR3BP::TPAT_Sys_Data_CR3BP(const char *filepath){
	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw TPAT_Exception("TPAT_Sys_Data_CR3BP: Could not load data from file");
	}
	readFromMat(matfp);
	Mat_Close(matfp);
}//===================================================

/**
 *  @brief Destructor
 */
TPAT_Sys_Data_CR3BP::~TPAT_Sys_Data_CR3BP(){}

/**
 *	@brief Initialize all data fields using the names of the primaries
 *	@param P1 the name of the larger primary
 *	@param P2 the name of the smaller primary
 */
void TPAT_Sys_Data_CR3BP::initFromPrimNames(std::string P1, std::string P2){
	TPAT_Body_Data p1Data(P1);
	TPAT_Body_Data p2Data(P2);

	primaries.clear();
	primIDs.clear();

	primaries.push_back(p1Data.getName());
	primIDs.push_back(p1Data.getID());
	primaries.push_back(p2Data.getName());
	primIDs.push_back(p2Data.getID());

	// Check to make sure P1 is P2's parent
	if(p2Data.getParent().compare(p1Data.getName()) == 0){
		charL = p2Data.getOrbitRad();
		charM = p1Data.getMass() + p2Data.getMass();
		charT = sqrt(pow(charL, 3)/(G*charM));

		otherParams.at(0) = p2Data.getMass()/charM;	// Non-dimensional mass ratio mu
	}else{
		throw TPAT_Exception("P1 must be the parent of P2");
	}
}//===================================================

/**
 *	@brief Retrieve the model that governs the motion for this system type
 *	@return the model that governs the motion for this system type
 */
const TPAT_Model* TPAT_Sys_Data_CR3BP::getModel() const { return &model; }

/**
 *	@brief Copy constructor
 *	@param d
 */
TPAT_Sys_Data_CR3BP::TPAT_Sys_Data_CR3BP(const TPAT_Sys_Data_CR3BP &d) : TPAT_Sys_Data(d){}

/**
 *	@brief Copy operator; makes a clean copy of a data object into this one
 *	@param d a CR3BP system data object
 *	@return this system data object
 */
TPAT_Sys_Data_CR3BP& TPAT_Sys_Data_CR3BP::operator= (const TPAT_Sys_Data_CR3BP &d){
	TPAT_Sys_Data::operator= (d);
	return *this;
}//===================================================

/**
 *	@return the non-dimensional mass ratio for the system
 */
double TPAT_Sys_Data_CR3BP::getMu() const { return otherParams.at(0); }

/**
 *  @brief Save the system data to a matlab file
 * 
 *  @param filepath path to the data file
 */
void TPAT_Sys_Data_CR3BP::saveToMat(const char *filepath) const{
	TPAT_Sys_Data::saveToMat(filepath);
}//==================================================

/**
 *	@brief Save system data, like the names of the primaries and the system mass ratio, to a .mat file
 *	@param matFile a pointer to the .mat file
 */
void TPAT_Sys_Data_CR3BP::saveToMat(mat_t *matFile) const{
	size_t dims[2] = {1,1};

	if(primaries.size() < 2){
		printErr("Primaries size is %zu\n", primaries.size());
		throw TPAT_Exception("TPAT_Sys_Data_CR3BP::saveToMat: There are no primaries?");
	}

	// Initialize character array (larger than needed), copy in the name of the primary, then create a var.
	char p1_str[64];
	strcpy(p1_str, primaries.at(0).c_str());
	dims[1] = primaries.at(0).length();
	matvar_t *p1_var = Mat_VarCreate("P1", MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(p1_str[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p1_var, "P1", MAT_COMPRESSION_NONE);

	char p2_str[64];
	strcpy(p2_str, primaries.at(1).c_str());
	dims[1] = primaries.at(1).length();
	matvar_t *p2_var = Mat_VarCreate("P2", MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(p2_str[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p2_var, "P2", MAT_COMPRESSION_NONE);

	dims[1] = 1;
	double mu = otherParams[0];
	matvar_t *mu_var = Mat_VarCreate("Mu", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &mu, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, mu_var, "Mu", MAT_COMPRESSION_NONE);
}//===================================================

/**
 *	@brief Populate data fiels for this data object by reading the primaries'
 *	names from a Mat file
 *	@param matFile a pointer to the Mat file in question
 */
void TPAT_Sys_Data_CR3BP::readFromMat(mat_t *matFile){
	std::string P1 = readStringFromMat(matFile, "P1", MAT_T_UINT8, MAT_C_CHAR);
	std::string P2 = readStringFromMat(matFile, "P2", MAT_T_UINT8, MAT_C_CHAR);

	numPrimaries = 2;
	type = TPAT_System_Tp::CR3BP_SYS;
	otherParams.assign(1,0);	// make mu = 0
	initFromPrimNames(P1, P2);
}//===================================================
/**
 *	\file SysData_cr3bp.cpp
 *	\brief Derivative of SysData, specific to CR3BP
 *	
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
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

#include "SysData_cr3bp.hpp"
 
#include "BodyData.hpp"
#include "Common.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"

#include <cstring>
#include <cmath>
#include <exception>

namespace astrohelion{
/**
 *	\brief Default constructor
 */
SysData_cr3bp::SysData_cr3bp() : SysData(){
	numPrimaries = 2;
	type = SysData_tp::CR3BP_SYS;
	otherParams.assign(1,0);	// make mu = 0
}//========================================

/**
 *	\brief Create a system data object using data from the two primaries
 *	\param P1 the name of the larger primary
 *	\param P2 the name of the smaller primary; P2 must orbit P1
 */
SysData_cr3bp::SysData_cr3bp(std::string P1, std::string P2){
	numPrimaries = 2;
	type = SysData_tp::CR3BP_SYS;
	otherParams.assign(1,0);	// make mu = 0
	
	initFromPrimNames(P1, P2);
}//===================================================

/**
 *  \brief Load the system data object from a Matlab data file
 * 
 *  \param filepath path to the data file
 */
SysData_cr3bp::SysData_cr3bp(const char *filepath){
	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw Exception("SysData_cr3bp: Could not load data from file");
	}
	readFromMat(matfp);
	Mat_Close(matfp);
}//===================================================

/**
 *  \brief Destructor
 */
SysData_cr3bp::~SysData_cr3bp(){}

/**
 *	\brief Initialize all data fields using the names of the primaries
 *	\param P1 the name of the larger primary
 *	\param P2 the name of the smaller primary
 */
void SysData_cr3bp::initFromPrimNames(std::string P1, std::string P2){
	BodyData p1Data(P1);
	BodyData p2Data(P2);

	primaries.clear();
	primIDs.clear();

	primaries.push_back(p1Data.getName());
	primIDs.push_back(p1Data.getID());
	primaries.push_back(p2Data.getName());
	primIDs.push_back(p2Data.getID());

	// Check to make sure P1 is P2's parent
	if(p2Data.getParent().compare(p1Data.getName()) == 0){
		double totalGM = p1Data.getGravParam() + p2Data.getGravParam();

		charL = p2Data.getOrbitRad();
		charM = totalGM/G;
		charT = sqrt(pow(charL, 3)/totalGM);

		otherParams.at(0) = p2Data.getGravParam()/totalGM;	// Non-dimensional mass ratio mu
	}else{
		throw Exception("P1 must be the parent of P2");
	}
}//===================================================

/**
 *	\brief Retrieve the model that governs the motion for this system type
 *	\return the model that governs the motion for this system type
 */
const DynamicsModel* SysData_cr3bp::getDynamicsModel() const { return &model; }

/**
 *  \brief Retrieve the object that serves up control law data for this system type
 *  \return the control law object
 */
const ControlLaw* SysData_cr3bp::getControlLaw() const { return &control; }

/**
 *	\brief Copy constructor
 *	\param d
 */
SysData_cr3bp::SysData_cr3bp(const SysData_cr3bp &d) : SysData(d){}

/**
 *	\brief Copy operator; makes a clean copy of a data object into this one
 *	\param d a CR3BP system data object
 *	\return this system data object
 */
SysData_cr3bp& SysData_cr3bp::operator= (const SysData_cr3bp &d){
	SysData::operator= (d);
	return *this;
}//===================================================

/**
 *	\return the non-dimensional mass ratio for the system
 */
double SysData_cr3bp::getMu() const { return otherParams.at(0); }

/**
 *  \brief Save the system data to a matlab file
 * 
 *  \param filepath path to the data file
 */
void SysData_cr3bp::saveToMat(const char *filepath) const{
	SysData::saveToMat(filepath);
}//==================================================

/**
 *	\brief Save system data, like the names of the primaries and the system mass ratio, to a .mat file
 *	\param matFile a pointer to the .mat file
 */
void SysData_cr3bp::saveToMat(mat_t *matFile) const{
	size_t dims[2] = {1,1};

	if(primaries.size() < 2){
		astrohelion::printErr("Primaries size is %zu\n", primaries.size());
		throw Exception("SysData_cr3bp::saveToMat: There are no primaries?");
	}

	// Initialize character array (larger than needed), copy in the name of the primary, then create a var.
	char p1_str[64];
	strcpy(p1_str, primaries.at(0).c_str());
	dims[1] = primaries.at(0).length();
	matvar_t *p1_var = Mat_VarCreate("P1", MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(p1_str[0]), MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, p1_var, "P1", MAT_COMPRESSION_NONE);

	char p2_str[64];
	strcpy(p2_str, primaries.at(1).c_str());
	dims[1] = primaries.at(1).length();
	matvar_t *p2_var = Mat_VarCreate("P2", MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(p2_str[0]), MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, p2_var, "P2", MAT_COMPRESSION_NONE);

	dims[1] = 1;
	double mu = otherParams[0];
	matvar_t *mu_var = Mat_VarCreate("Mu", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &mu, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, mu_var, "Mu", MAT_COMPRESSION_NONE);
}//===================================================

/**
 *	\brief Populate data fiels for this data object by reading the primaries'
 *	names from a Mat file
 *	\param matFile a pointer to the Mat file in question
 */
void SysData_cr3bp::readFromMat(mat_t *matFile){
	std::string P1 = astrohelion::readStringFromMat(matFile, "P1", MAT_T_UINT8, MAT_C_CHAR);
	std::string P2 = astrohelion::readStringFromMat(matFile, "P2", MAT_T_UINT8, MAT_C_CHAR);

	numPrimaries = 2;
	type = SysData_tp::CR3BP_SYS;
	otherParams.assign(1,0);	// make mu = 0
	initFromPrimNames(P1, P2);
}//===================================================


}// END of Astrohelion namespace
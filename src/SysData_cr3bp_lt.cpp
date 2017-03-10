/**
 *	\file SysData_cr3bp_lt.cpp
 *	\brief Derivative of SysData, specific to CR3BP-LTVP
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

#include "SysData_cr3bp_lt.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"

#include <cstring>


namespace astrohelion{
/**
 *	\brief Default constructor
 */
SysData_cr3bp_lt::SysData_cr3bp_lt() : SysData_cr3bp(){
	type = SysData_tp::CR3BP_LT_SYS;
	otherParams.assign(4,0);
}//========================================

/**
 *	\brief Create a system data object using data from the two primaries
 *	\param P1 the name of the larger primary
 *	\param P2 the name of the smaller primary; P2 must orbit P1
 *	\param T Thrust value, Newtons
 *	\param I Specific Impulse (Isp), seconds
 *	\param M0 reference mass, kilograms
 */
SysData_cr3bp_lt::SysData_cr3bp_lt(std::string P1, std::string P2, double T, double I, double M0){
	numPrimaries = 2;
	type = SysData_tp::CR3BP_LT_SYS;
	otherParams.assign(4,0);
	
	initFromPrimNames(P1, P2);	// use function from cr3bp_sys_data to initialize most everything
	otherParams[1] = (T/1000)*charT*charT/charL/M0;		// thrust, non-dimensionalized
	otherParams[2] = I;									// Isp, dimensional
	otherParams[3] = M0;								// Initial mass, dimensional
}//===================================================

/**
 *  \brief Load the system data object from a Matlab data file
 * 
 *  \param filepath path to the data file
 */
SysData_cr3bp_lt::SysData_cr3bp_lt(const char *filepath){
	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw Exception("SysData_cr3bp_lt: Could not load data from file");
	}
	readFromMat(matfp);
	Mat_Close(matfp);
}//===================================================

/**
 *	\brief Copy constructor
 *	\param d
 */
SysData_cr3bp_lt::SysData_cr3bp_lt(const SysData_cr3bp_lt &d) : SysData_cr3bp(d){}

/**
 *	\brief Copy operator; makes a clean copy of a data object into this one
 *	\param d a CR3BP system data object
 *	\return this system data object
 */
SysData_cr3bp_lt& SysData_cr3bp_lt::operator= (const SysData_cr3bp_lt &d){
	SysData_cr3bp::operator= (d);
	return *this;
}//===================================================

/**
 *	\brief Retrieve the model that governs the motion for this system type
 *	\return the model that governs the motion for this system type
 */
const DynamicsModel* SysData_cr3bp_lt::getDynamicsModel() const { return &model; }

/**
 *  \brief Retrieve the object that serves up control law data for this system type
 *  \return the control law object
 */
const ControlLaw* SysData_cr3bp_lt::getControlLaw() const { return &control; }

/**
 *	\brief Get the nondimensional thrust for P3 in this system
 *	\return the thrust, nondimensional (nondimensionalized in mass by spacecraft reference mass) 
 */
double SysData_cr3bp_lt::getThrust() const { return otherParams[1]; }

/**
 *	\brief Get the specific impulse for P3 in this system
 *	\return the specific impulse, seconds
 */
double SysData_cr3bp_lt::getIsp() const { return otherParams[2]; }

/**
 *	\brief Get the reference mass for P3
 *	\return the reference mass of P3, kg
 */
double SysData_cr3bp_lt::getMass() const { return otherParams[3]; }

/**
 *	\brief Set the thrust for P3 for this system
 *	\param d the thrust, non-dimensional units (nondimensionalized in mass by spacecraft reference mass)
 */
void SysData_cr3bp_lt::setThrust(double d){ otherParams[1] = d; }

/**
 *	\brief Set the specific impulse for P3 for this system
 *	\param d the specific impulse, seconds
 */
void SysData_cr3bp_lt::setIsp(double d){ otherParams[2] = d; }

/**
 *	\brief Set the reference mass for P3
 *	\param d the mass, kg
 */
void SysData_cr3bp_lt::setMass(double d){ otherParams[3] = d; }

/**
 *  \brief Save the system data to a matlab file
 * 
 *  \param filepath path to the data file
 */
void SysData_cr3bp_lt::saveToMat(const char *filepath) const{
	SysData::saveToMat(filepath);
}//==================================================

/**
 *	\brief Save system data, like the names of the primaries and the system mass ratio, to a .mat file
 *	\param matFile a pointer to the .mat file
 *	\throws Exception if the number of primaries is incorrect
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
	astrohelion::saveVar(matFile, p1_var, "P1", MAT_COMPRESSION_NONE);

	char p2_str[64];
	std::strcpy(p2_str, primaries.at(1).c_str());
	dims[1] = primaries.at(1).length();
	matvar_t *p2_var = Mat_VarCreate("P2", MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(p2_str[0]), MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, p2_var, "P2", MAT_COMPRESSION_NONE);

	dims[1] = 1;	
	double mu = otherParams[0];
	matvar_t *mu_var = Mat_VarCreate("Mu", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &mu, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, mu_var, "Mu", MAT_COMPRESSION_NONE);

	// Save thrust in nondimensional units
	double T = otherParams[1];//*1000*charL*otherParams[3]/charT/charT;
	matvar_t *thrust_var = Mat_VarCreate("Thrust", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &T, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, thrust_var, "Thrust", MAT_COMPRESSION_NONE);

	// Save Isp in dimensional units
	double I = otherParams[2];
	matvar_t *imp_var = Mat_VarCreate("Isp", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &I, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, imp_var, "Isp", MAT_COMPRESSION_NONE);

	// Save mass in dimensional units
	double m = otherParams[3];
	matvar_t *mass_var = Mat_VarCreate("Mass0", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &m, MAT_F_DONT_COPY_DATA);
	astrohelion::saveVar(matFile, mass_var, "Mass0", MAT_COMPRESSION_NONE);
}//===================================================

/**
 *	\brief Populate data fiels for this data object by reading the primaries'
 *	names from a Mat file
 *	\param matFile a pointer to the Mat file in question
 */
void SysData_cr3bp_lt::readFromMat(mat_t *matFile){
	std::string P1 = astrohelion::readStringFromMat(matFile, "P1", MAT_T_UINT8, MAT_C_CHAR);
	std::string P2 = astrohelion::readStringFromMat(matFile, "P2", MAT_T_UINT8, MAT_C_CHAR);

	type = SysData_tp::CR3BP_LT_SYS;
	otherParams.assign(4,0);
	initFromPrimNames(P1, P2);
	otherParams[1] = astrohelion::readDoubleFromMat(matFile, "Thrust");
	otherParams[2] = astrohelion::readDoubleFromMat(matFile, "Isp");
	otherParams[3] = astrohelion::readDoubleFromMat(matFile, "Mass0");
}//===================================================




}// END of Astrohelion namespace
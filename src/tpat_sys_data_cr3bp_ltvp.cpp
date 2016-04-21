/**
 *	@file tpat_sys_data_cr3bp_ltvp.cpp
 *	@brief Derivative of tpat_sys_data, specific to CR3BP-LTVP
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

#include "tpat_sys_data_cr3bp_ltvp.hpp"
#include "tpat_utilities.hpp"

#include <cstring>
/**
 *	@brief Default constructor
 */
tpat_sys_data_cr3bp_ltvp::tpat_sys_data_cr3bp_ltvp() : tpat_sys_data_cr3bp(){
	type = tpat_sys_data::CR3BP_LTVP_SYS;
	otherParams.assign(4,0);
}//========================================

/**
 *	@brief Create a system data object using data from the two primaries
 *	@param P1 the name of the larger primary
 *	@param P2 the name of the smaller primary; P2 must orbit P1
 *	@param T Thrust value, Newtons
 *	@param I Specific Impulse (Isp), seconds
 *	@param M0 initial mass (at t = 0), kilograms
 */
tpat_sys_data_cr3bp_ltvp::tpat_sys_data_cr3bp_ltvp(std::string P1, std::string P2, double T, double I, double M0){
	numPrimaries = 2;
	type = tpat_sys_data::CR3BP_LTVP_SYS;
	otherParams.assign(4,0);
	
	initFromPrimNames(P1, P2);	// use function from cr3bp_sys_data to initialize most everything
	otherParams[1] = (T/1000)*charT*charT/charL/charM;	// thrust, non-dimensionalized
	otherParams[2] = I/charT;	// Isp, non-dimensionalized
	otherParams[3] = M0/charM;	// Initial mass, non-dimensionalized
}//===================================================

/**
 *  @brief Load the system data object from a Matlab data file
 * 
 *  @param filepath path to the data file
 */
tpat_sys_data_cr3bp_ltvp::tpat_sys_data_cr3bp_ltvp(const char *filepath){
	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw tpat_exception("tpat_sys_data_cr3bp_ltvp: Could not load data from file");
	}
	readFromMat(matfp);
	Mat_Close(matfp);
}//===================================================

/**
 *	@brief Copy constructor
 *	@param d
 */
tpat_sys_data_cr3bp_ltvp::tpat_sys_data_cr3bp_ltvp(const tpat_sys_data_cr3bp_ltvp &d) : tpat_sys_data_cr3bp(d){}

/**
 *	@brief Copy operator; makes a clean copy of a data object into this one
 *	@param d a CR3BP system data object
 *	@return this system data object
 */
tpat_sys_data_cr3bp_ltvp& tpat_sys_data_cr3bp_ltvp::operator= (const tpat_sys_data_cr3bp_ltvp &d){
	tpat_sys_data_cr3bp::operator= (d);
	return *this;
}//===================================================

/**
 *	@brief Retrieve the model that governs the motion for this system type
 *	@return the model that governs the motion for this system type
 */
const tpat_model* tpat_sys_data_cr3bp_ltvp::getModel() const { return &model; }

/**
 *	@brief Get the non-dimensional thrust for P3 in this system
 *	@return the thrust (non-dimensional)
 */
double tpat_sys_data_cr3bp_ltvp::getThrust() const { return otherParams[1]; }

/**
 *	@brief Get the non-dimensional specific impulse for P3 in this system
 *	@return the specific impulse (non-dimensional)
 */
double tpat_sys_data_cr3bp_ltvp::getIsp() const { return otherParams[2]; }

/**
 *	@brief Get the non-dimensional initial mass for P3
 *	@return the non-dimensional initial mass for P3
 */
double tpat_sys_data_cr3bp_ltvp::getM0() const { return otherParams[3]; }

/**
 *	@brief Set the thrust for P3 for this system
 *	@param d the thrust, non-dimensional units
 */
void tpat_sys_data_cr3bp_ltvp::setThrust(double d){ otherParams[1] = d; }

/**
 *	@brief Set the specific impulse for P3 for this system
 *	@param d the specific impulse, non-dimensional units
 */
void tpat_sys_data_cr3bp_ltvp::setIsp(double d){ otherParams[2] = d; }

/**
 *	@brief Set the initial mass for P3 (non-dim)
 *	@param d the initial mass, non-dimensional units
 */
void tpat_sys_data_cr3bp_ltvp::setM0(double d){ otherParams[3] = d; }

/**
 *	@brief Set the thrust for P3 for this system using a dimensional quantity
 *	@param d the thrust, in Newtons
 */
void tpat_sys_data_cr3bp_ltvp::setThrustDim(double d){ otherParams[1] = (d/1000)*charT*charT/charM/charL; }

/**
 *	@brief Set the specific impulse for P3 for this system using a dimensional quantity
 *	@param d the specific impulse, in seconds
 */
void tpat_sys_data_cr3bp_ltvp::setIspDim(double d){ otherParams[2] = d/charT; }

/**
 *	@brief Set the initial mass for P3 using a dimensional quantity
 *	@param d the initial mass, in kilograms
 */
void tpat_sys_data_cr3bp_ltvp::setM0Dim(double d){ otherParams[3] = d/charM; }

/**
 *  @brief Save the system data to a matlab file
 * 
 *  @param filepath path to the data file
 */
void tpat_sys_data_cr3bp_ltvp::saveToMat(const char *filepath) const{
	tpat_sys_data::saveToMat(filepath);
}//==================================================

/**
 *	@brief Save system data, like the names of the primaries and the system mass ratio, to a .mat file
 *	@param matFile a pointer to the .mat file
 */
void tpat_sys_data_cr3bp_ltvp::saveToMat(mat_t *matFile) const{
	size_t dims[2] = {1,1};

	if(primaries.size() < 2){
		printErr("Primaries size is %zu\n", primaries.size());
		throw tpat_exception("tpat_sys_data_cr3bp::saveToMat: There are no primaries?");
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

	double T = otherParams[1];
	matvar_t *thrust_var = Mat_VarCreate("Thrust", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &T, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, thrust_var, "Thrust", MAT_COMPRESSION_NONE);

	double I = otherParams[2];
	matvar_t *imp_var = Mat_VarCreate("Isp", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &I, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, imp_var, "Isp", MAT_COMPRESSION_NONE);

	double m = otherParams[3];
	matvar_t *mass_var = Mat_VarCreate("Mass0", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &m, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, mass_var, "Mass0", MAT_COMPRESSION_NONE);
}//===================================================

/**
 *	@brief Populate data fiels for this data object by reading the primaries'
 *	names from a Mat file
 *	@param matFile a pointer to the Mat file in question
 */
void tpat_sys_data_cr3bp_ltvp::readFromMat(mat_t *matFile){
	std::string P1 = readStringFromMat(matFile, "P1", MAT_T_UINT8, MAT_C_CHAR);
	std::string P2 = readStringFromMat(matFile, "P2", MAT_T_UINT8, MAT_C_CHAR);

	type = tpat_sys_data::CR3BP_LTVP_SYS;
	otherParams.assign(4,0);
	initFromPrimNames(P1, P2);
	otherParams[2] = readDoubleFromMat(matFile, "Isp");
	otherParams[3] = readDoubleFromMat(matFile, "Mass0");
}//===================================================
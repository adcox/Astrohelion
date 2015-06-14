/**
 *	@tpat_cr3bp_sys_data.cpp
 *
 *	tpat_cr3bp_sys_data.cpp
 *
 * 	System Data object specifically for CR3BP
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

#include "tpat_cr3bp_sys_data.hpp"
 
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
tpat_cr3bp_sys_data::tpat_cr3bp_sys_data() : tpat_sys_data(){
	numPrimaries = 2;
	type = tpat_sys_data::CR3BP_SYS;
}//========================================

/**
 *	@brief Create a system data object using data from the two primaries
 *	@param P1 the name of the larger primary
 *	@param P2 the name of the smaller primary; P2 must orbit P1
 */
tpat_cr3bp_sys_data::tpat_cr3bp_sys_data(std::string P1, std::string P2){
	numPrimaries = 2;
	type = tpat_sys_data::CR3BP_SYS;
	
	tpat_body_data p1Data(P1);
	tpat_body_data p2Data(P2);

	primaries.push_back(p1Data.getName());
	primIDs.push_back(p1Data.getID());
	primaries.push_back(p2Data.getName());
	primIDs.push_back(p2Data.getID());

	// Check to make sure P1 is P2's parent
	if(p2Data.getName().compare(p1Data.getName())){
		charL = p2Data.getOrbitRad();
		charM = p1Data.getMass() + p2Data.getMass();
		charT = sqrt(pow(charL, 3)/(G*charM));

		mu = p2Data.getMass()/charM;
	}else{
		throw tpat_exception("P1 must be the parent of P2");
	}
}//===================================================

/**
 *	@brief Copy constructor
 *	@param d
 */
tpat_cr3bp_sys_data::tpat_cr3bp_sys_data(const tpat_cr3bp_sys_data &d) : tpat_sys_data(d){
	mu = d.mu;
}

/**
 *	@brief Copy operator; makes a clean copy of a data object into this one
 *	@param d a CR3BP system data object
 *	@return this system data object
 */
tpat_cr3bp_sys_data& tpat_cr3bp_sys_data::operator= (const tpat_cr3bp_sys_data &d){
	tpat_sys_data::operator= (d);
	mu = d.mu;
	return *this;
}//===================================================

/**
 *	@return the non-dimensional mass ratio for the system
 */
double tpat_cr3bp_sys_data::getMu() const { return mu; }

/**
 *	@brief Save system data, like the names of the primaries and the system mass ratio, to a .mat file
 *	@param matFile a pointer to the .mat file
 */
void tpat_cr3bp_sys_data::saveToMat(mat_t *matFile){
	size_t dims[2] = {1,1};

	// Initialize character array (larger than needed), copy in the name of the primary, then create a var.
	char p1_str[64];
	strcpy(p1_str, primaries.at(0).c_str());
	dims[1] = primaries.at(0).length();
	matvar_t *p1_var = Mat_VarCreate("P1", MAT_C_CHAR, MAT_T_UTF8, 2, dims, &(p1_str[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p1_var, "P1", MAT_COMPRESSION_NONE);

	char p2_str[64];
	strcpy(p2_str, primaries.at(1).c_str());
	dims[1] = primaries.at(1).length();
	matvar_t *p2_var = Mat_VarCreate("P2", MAT_C_CHAR, MAT_T_UTF8, 2, dims, &(p2_str[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p2_var, "P2", MAT_COMPRESSION_NONE);

	dims[1] = 1;	
	matvar_t *mu_var = Mat_VarCreate("Mu", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &mu, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, mu_var, "Mu", MAT_COMPRESSION_NONE);
}//===================================================
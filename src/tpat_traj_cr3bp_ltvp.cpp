/**
 *	@file tpat_traj_cr3bp_ltvp.cpp
 *	@brief Derivative of tpat_traj, specific to CR3BP-LTVP
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

#include "tpat_traj_cr3bp_ltvp.hpp"
#include "tpat_utilities.hpp"
 
//-----------------------------------------------------
// 		Constructor Functions
//-----------------------------------------------------

/**
 *	@brief Default constructor; calls constructor for super-class tpat_traj
 *	and additionall initializes the jacobi matrix
 */
tpat_traj_cr3bp_ltvp::tpat_traj_cr3bp_ltvp() : tpat_traj(){
	numExtraParam = 1;	// Jacobi
	extraParamRowSize.push_back(1);	// each jacobi value has one element
	extraParam.push_back(std::vector<double>(0));
}//====================================================

/**
 *	@brief Create a CR3BP LTVP trajectory object for the specified system
 *	@param data a system data object describing the system
 */
tpat_traj_cr3bp_ltvp::tpat_traj_cr3bp_ltvp(tpat_sys_data_cr3bp_ltvp data){
	numExtraParam = 1;	// Jacobi
	extraParamRowSize.push_back(1);	// each jacobi value has one element
	extraParam.push_back(std::vector<double>(0));
	sysData = data;
}//====================================================

/**
 *	@brief Initialize all vectors to have size n; fill each vector with zeros.
 */
tpat_traj_cr3bp_ltvp::tpat_traj_cr3bp_ltvp(int n) : tpat_traj(n){
	numExtraParam = 1;	// Jacobi
	extraParamRowSize.push_back(1);	// each jacobi value has one element
	extraParam.push_back(std::vector<double>(0));
	extraParam.at(0).reserve(n);
}//====================================================

/**
 *	@brief Copy the specified trajectory
 *	@param t trajectory
 */
tpat_traj_cr3bp_ltvp::tpat_traj_cr3bp_ltvp(const tpat_traj_cr3bp_ltvp &t) : tpat_traj(t){
	copyMe(t);
}//====================================================

//-----------------------------------------------------
// 		Operators
//-----------------------------------------------------

/**
 *	@brief Copy operator; copy a trajectory object into this one.
 *	@param t a trajectory object
 *	@return this trajectory object
 */
tpat_traj_cr3bp_ltvp& tpat_traj_cr3bp_ltvp::operator= (const tpat_traj_cr3bp_ltvp& t){
	tpat_traj::operator= (t);
	copyMe(t);
	return *this;
}//=====================================================

/**
 *	@brief Copy all data objects specific to the CR3BP Trajectory
 *	@param t the source trajectory; copy its attributes into this one
 */
void tpat_traj_cr3bp_ltvp::copyMe(const tpat_traj_cr3bp_ltvp &t){
	sysData = t.sysData;
}//=====================================================

/**
 *	@brief Retrieve data about the system this trajectory was propagated in
 *	@return the system data object
 */
tpat_sys_data_cr3bp_ltvp tpat_traj_cr3bp_ltvp::getSysData() const { return sysData; }


/**
 *	@brief Retrieve the system type for this trajectory
 *	@return the system type
 */
tpat_sys_data::system_t tpat_traj_cr3bp_ltvp::getType() const{
	return sysData.getType();
}//=====================================================

/**
 *	@brief Get a pointer to the vector of Jacobi values along the trajectory
 *	@return a pointer to the vector of Jacobi values along the trajectory
 */
std::vector<double>* tpat_traj_cr3bp_ltvp::getJacobi() { return &(extraParam.at(0)); }

/**
 *	@brief Retrieve a Jacobi value at one point on the trajectory
 *	@param n the index of the point; if n < 0, it will count backwards
 *	from the end of the trajectory
 *	@return a Jacobi value (non-dimensional units)
 */
double tpat_traj_cr3bp_ltvp::getJacobi(int n) const {
	if(n < 0)
		n += extraParam.at(0).size();

	return extraParam.at(0)[n];
}//=====================================================

/**
 *	@brief Set the system data for this trajectory
 *	@param d a system data object
 */
void tpat_traj_cr3bp_ltvp::setSysData(tpat_sys_data_cr3bp_ltvp d){ sysData = d; }

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_traj_cr3bp_ltvp::saveToMat(const char* filename){
	// TODO: Check for propper file extension, add if necessary

	/*	Create a new Matlab MAT file with the given name and optional
	 *	header string. If no header string is given, the default string 
	 *	used containing the software, version, and date in it. If a header
	 *	string is specified, at most the first 116 characters are written to
	 *	the file. Arguments are:
	 *	const char *matname 	- 	the name of the file
	 *	const char *hdr_str 	- 	the 116 byte header string
	 *	enum mat_ft 			- 	matlab file @version MAT_FT_MAT5 or MAT_FT_MAT4
	 */
	mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		printErr("Error creating MAT file\n");
	}else{
		saveState(matfp);
		saveTime(matfp);
		saveSTMs(matfp);
		saveJacobi(matfp);
		sysData.saveToMat(matfp);
	}

	Mat_Close(matfp);
}//========================================

/**
 *	@brief Save the Jacobi vector to a file
 * 	@param matFile a pointer to the destination matlab file 
 */
void tpat_traj_cr3bp_ltvp::saveJacobi(mat_t *matFile){
	size_t dims[2] = {extraParam.size(), 1};
	matvar_t *matvar = Mat_VarCreate("Jacobi", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(extraParam[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "Jacobi", MAT_COMPRESSION_NONE);
}//=================================================


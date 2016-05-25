/**
 *  @file tpat_traj_bc4bp.cpp
 *	@brief Derivative of TPAT_Traj, specific to BCR4BPR
 *
 *	@author Andrew Cox
 *	@version September 2, 2015
 *	@copyright GNU GPL v3.0
 */
 
/*
 *  Trajectory Propagation and Analysis Toolkit 
 *  Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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

#include "tpat_traj_bc4bp.hpp"

#include "tpat_exceptions.hpp"
#include "tpat_node.hpp"
#include "tpat_sim_engine.hpp"
#include "tpat_sys_data_bc4bp.hpp"
#include "tpat_utilities.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Create a trajectory for a specific system
 *	@param sys a pointer to a system data object
 */
TPAT_Traj_BC4BP::TPAT_Traj_BC4BP(const TPAT_Sys_Data_BC4BP *sys) : TPAT_Traj(sys){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from another trajectory
 *	@param t a trajectory reference
 */
TPAT_Traj_BC4BP::TPAT_Traj_BC4BP(const TPAT_Traj_BC4BP &t) : TPAT_Traj(t){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from its base class
 *	@param a an arc data reference
 */
TPAT_Traj_BC4BP::TPAT_Traj_BC4BP(const TPAT_Base_Arcset &a) : TPAT_Traj(a){
	initExtraParam();
}//====================================================

/**
 *  @brief Create a new trajectory object on the stack
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @param sys pointer to a system data object; should be a 
 *  BCR4BPR system as the pointer will be cast to that derived class
 *  @return a pointer to the newly created trajectory
 */
baseArcsetPtr TPAT_Traj_BC4BP::create( const TPAT_Sys_Data *sys) const{
	const TPAT_Sys_Data_BC4BP *bcSys = static_cast<const TPAT_Sys_Data_BC4BP*>(sys);
	return baseArcsetPtr(new TPAT_Traj_BC4BP(bcSys));
}//====================================================

/**
 *  @brief Create a new trajectory object on the stack that is a 
 *  duplicate of this object
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @return a pointer to the newly cloned trajectory
 */
baseArcsetPtr TPAT_Traj_BC4BP::clone() const{
	return baseArcsetPtr(new TPAT_Traj_BC4BP(*this));
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@return the angle between the P1/P2 line and the inertial x-axis, radians
 */
double TPAT_Traj_BC4BP::getTheta0(){
	const TPAT_Sys_Data_BC4BP *bcSys = static_cast<const TPAT_Sys_Data_BC4BP *>(sysData);
	return bcSys->getTheta0();
}//====================================================

/**
 *	@return the angle between the P2/P3 line (projected into the inertial XY plane)
 *	and the inertial x-axis, radians
 */
double TPAT_Traj_BC4BP::getPhi0(){
	const TPAT_Sys_Data_BC4BP *bcSys = static_cast<const TPAT_Sys_Data_BC4BP *>(sysData);
	return bcSys->getPhi0();
}//====================================================

/**
 *	@return the inclination of the P2/P3 orbital plane relative to the P1/P2 orbital
 *	plane, radians
 */
double TPAT_Traj_BC4BP::getGamma(){
	const TPAT_Sys_Data_BC4BP *bcSys = static_cast<const TPAT_Sys_Data_BC4BP *>(sysData);
	return bcSys->getGamma();
}//====================================================

/**
 *	@param ix the index of the dqdT vector to retrieve
 *	@return the i'th 6-element dqdT vector. If ix is negative, the count
 *	will proceed from the end of the vector, i.e. -1 will return the final time, 
 *	-2 will give the second to last value, etc.
 *	@throws TPAT_Exception if <tt>ix</tt> is out of bounds
 */
std::vector<double> TPAT_Traj_BC4BP::get_dqdTByIx(int ix){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > ((int)nodes.size()))
		throw TPAT_Exception("TPAT_Traj_BC4BP::getdqdT: invalid index");

	return getExtraParam(ix, 0);
}//====================================================

/**
 *	@brief Set the value of the dqdT vector for the specified step
 *	@param ix the index of the step; if < 0, it will count backwards from the end
 *	@param dqdT a pointer to the dqdT vector; this MUST have at least 6 elements,
 *	or the function will read unallocated memory.
 *	@throws TPAT_Exception if <tt>ix</tt> is out of bounds
 */
void TPAT_Traj_BC4BP::set_dqdTByIx(int ix, const double *dqdT){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > ((int)nodes.size()))
		throw TPAT_Exception("TPAT_Traj_BC4BP::setdqdT: invalid index");

	for(int i = 0; i < 6; i++)
		nodes[ix].setExtraParam(0+i, dqdT[i]);
}//====================================================

/**
 *	@brief Set the value of the dqdT vector for the specified step
 *	@param ix the index of the step; if < 0, it will count backwards from the end
 *	@param dqdT a vector (6 elements) representing the dqdT vector
 *	@throws TPAT_Exception if <tt>ix</tt> is out of bounds
 */
void TPAT_Traj_BC4BP::set_dqdTByIx(int ix, std::vector<double> dqdT){
	if(dqdT.size() != 6)
		throw TPAT_Exception("TPAT_Traj_BC4BP::set_dqdT: Cannot accept a dqdT with anything other than 6 elements");

	set_dqdTByIx(ix, &(dqdT[0]));
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Initialize the extra param vector for info specific to this trajectory
 */
void TPAT_Traj_BC4BP::initExtraParam(){
	// Add a variable for dqdT
	numExtraParam = 1;
	extraParamRowSize.push_back(6);
}//====================================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void TPAT_Traj_BC4BP::saveToMat(const char* filename) const{
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
		saveAccel(matfp);
		saveEpoch(matfp, "Time");
		saveSTMs(matfp);
		saveExtraParam(matfp, 0, "dqdT");
		sysData->saveToMat(matfp);
	}

	Mat_Close(matfp);
}//========================================

/**
 *  @brief Populate data in this nodeset from a matlab file
 * 
 *  @param filepath the path to the matlab data file
 *  @throws TPAT_Exception if the Matlab file cannot be opened
 */
void TPAT_Traj_BC4BP::readFromMat(const char *filepath){
	TPAT_Traj::readFromMat(filepath);
	
	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw TPAT_Exception("TPAT_Traj_BC4BP: Could not load data from file");
	}

	readExtraParamFromMat(matfp, 0, "dqdT");
	
	Mat_Close(matfp);
}//====================================================



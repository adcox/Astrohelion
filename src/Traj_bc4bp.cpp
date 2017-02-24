/**
 *  @file Traj_bc4bp.cpp
 *	@brief Derivative of Traj, specific to BCR4BPR
 *
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of Astrohelion
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

#include "Traj_bc4bp.hpp"

#include "Exceptions.hpp"
#include "Node.hpp"
#include "SimEngine.hpp"
#include "SysData_bc4bp.hpp"
#include "Utilities.hpp"

namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Create a trajectory for a specific system
 *	@param pSys a pointer to a system data object
 */
Traj_bc4bp::Traj_bc4bp(const SysData_bc4bp *pSys) : Traj(pSys){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from another trajectory
 *	@param t a trajectory reference
 */
Traj_bc4bp::Traj_bc4bp(const Traj_bc4bp &t) : Traj(t){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from its base class
 *	@param a an arc data reference
 */
Traj_bc4bp::Traj_bc4bp(const BaseArcset &a) : Traj(a){
	initExtraParam();
}//====================================================

/**
 *  @brief Load the trajectory from a saved data file
 * 
 *  @param filepath Absolute or relative path to the data file
 *  @param pSys pointer to the system data object. Load the system object
 *  from the same file using the filepath constructor of the SysData_bc4bp
 *  object
 */
Traj_bc4bp::Traj_bc4bp(const char* filepath, const SysData_bc4bp *pSys) : Traj(pSys){
	initExtraParam();
	readFromMat(filepath);
}//====================================================

/**
 *  @brief Create a new trajectory object on the stack
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @param pSys pointer to a system data object; should be a 
 *  BCR4BPR system as the pointer will be cast to that derived class
 *  @return a pointer to the newly created trajectory
 */
baseArcsetPtr Traj_bc4bp::create( const SysData *pSys) const{
	const SysData_bc4bp *bcSys = static_cast<const SysData_bc4bp*>(pSys);
	return baseArcsetPtr(new Traj_bc4bp(bcSys));
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
baseArcsetPtr Traj_bc4bp::clone() const{
	return baseArcsetPtr(new Traj_bc4bp(*this));
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
double Traj_bc4bp::getTheta0(){
	const SysData_bc4bp *bcSys = static_cast<const SysData_bc4bp *>(pSysData);
	return bcSys->getTheta0();
}//====================================================

/**
 *	@return the angle between the P2/P3 line (projected into the inertial XY plane)
 *	and the inertial x-axis, radians
 */
double Traj_bc4bp::getPhi0(){
	const SysData_bc4bp *bcSys = static_cast<const SysData_bc4bp *>(pSysData);
	return bcSys->getPhi0();
}//====================================================

/**
 *	@return the inclination of the P2/P3 orbital plane relative to the P1/P2 orbital
 *	plane, radians
 */
double Traj_bc4bp::getGamma(){
	const SysData_bc4bp *bcSys = static_cast<const SysData_bc4bp *>(pSysData);
	return bcSys->getGamma();
}//====================================================

/**
 *	@param ix the index of the dqdT vector to retrieve
 *	@return the i'th 6-element dqdT vector. If ix is negative, the count
 *	will proceed from the end of the vector, i.e. -1 will return the final time, 
 *	-2 will give the second to last value, etc.
 *	@throws Exception if <tt>ix</tt> is out of bounds
 */
std::vector<double> Traj_bc4bp::get_dqdTByIx(int ix){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > static_cast<int>(nodes.size()))
		throw Exception("Traj_bc4bp::getdqdT: invalid index");

	return getExtraParamVec(ix, "dqdT");
}//====================================================

/**
 *	@brief Set the value of the dqdT vector for the specified step
 *	@param ix the index of the step; if < 0, it will count backwards from the end
 *	@param dqdT a pointer to the dqdT vector; this MUST have at least 6 elements,
 *	or the function will read unallocated memory.
 *	@throws Exception if <tt>ix</tt> is out of bounds
 */
void Traj_bc4bp::set_dqdTByIx(int ix, const double *dqdT){
	if(ix < 0)
		ix += nodes.size();

	if(ix < 0 || ix > static_cast<int>(nodes.size()))
		throw Exception("Traj_bc4bp::setdqdT: invalid index");

	std::vector<double> dqdT_vec(dqdT, dqdT+6);
	nodes[ix].setExtraParamVec("dqdT", dqdT_vec);
}//====================================================

/**
 *	@brief Set the value of the dqdT vector for the specified step
 *	@param ix the index of the step; if < 0, it will count backwards from the end
 *	@param dqdT a vector (6 elements) representing the dqdT vector
 *	@throws Exception if <tt>ix</tt> is out of bounds
 */
void Traj_bc4bp::set_dqdTByIx(int ix, std::vector<double> dqdT){
	if(dqdT.size() != 6)
		throw Exception("Traj_bc4bp::set_dqdT: Cannot accept a dqdT with anything other than 6 elements");

	set_dqdTByIx(ix, &(dqdT[0]));
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Initialize the extra param vector for info specific to this trajectory
 */
void Traj_bc4bp::initExtraParam(){
}//====================================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void Traj_bc4bp::saveToMat(const char* filename) const{
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
		astrohelion::printErr("Error creating MAT file\n");
	}else{
		saveState(matfp);
		saveAccel(matfp);
		saveEpoch(matfp, "Time");
		saveTOF(matfp, "TOFs");
		saveSTMs(matfp);
		saveExtraParamVec(matfp, "dqdT", 6, "dqdT");
		pSysData->saveToMat(matfp);
	}

	Mat_Close(matfp);
}//========================================

/**
 *  @brief Populate data in this nodeset from a matlab file
 * 
 *  @param filepath the path to the matlab data file
 *  @throws Exception if the Matlab file cannot be opened
 */
void Traj_bc4bp::readFromMat(const char *filepath){
	Traj::readFromMat(filepath);
	
	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw Exception("Traj_bc4bp: Could not load data from file");
	}

	readExtraParamVecFromMat(matfp, "dqdT", 6, "dqdT");
	
	Mat_Close(matfp);
}//====================================================


}// END of Astrohelion namespace
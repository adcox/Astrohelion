/**
 *  @file Traj_2bp.cpp
 *	@brief Derivative of Traj, specific to 2BP
 *
 *	@author Andrew Cox
 *	@version August 25, 2016
 *	@copyright GNU GPL v3.0
 */
 
/*
 *  Astrohelion 
 *  Copyright 2016, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "Traj_2bp.hpp"

#include "Exceptions.hpp"
#include "Node.hpp"
#include "SysData_2bp.hpp"
#include "Utilities.hpp"


namespace astrohelion{
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------
 
/**
 *	@brief Create a trajectory for a specific system
 *	@param sys a pointer to a system data object
 */
Traj_2bp::Traj_2bp(const SysData_2bp *sys) : Traj(sys){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from another trajectory
 *	@param t a trajectory reference
 */
Traj_2bp::Traj_2bp(const Traj_2bp &t) : Traj(t){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from its base class
 *	@param a an arc data reference
 */
Traj_2bp::Traj_2bp(const BaseArcset &a) : Traj(a){
	initExtraParam();
}//====================================================

/**
 *  @brief Create a new trajectory object on the stack
 *  @details the <tt>delete</tt> function must be called to 
 *  free the memory allocated to this object to avoid 
 *  memory leaks
 * 
 *  @param sys pointer to a system data object; should be a 
 *  CR3BP system as the pointer will be cast to that derived class
 *  @return a pointer to the newly created trajectory
 */
baseArcsetPtr Traj_2bp::create( const SysData *sys) const{
	const SysData_2bp *crSys = static_cast<const SysData_2bp*>(sys);
	return baseArcsetPtr(new Traj_2bp(crSys));
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
baseArcsetPtr Traj_2bp::clone() const{
	return baseArcsetPtr(new Traj_2bp(*this));
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

/**
 *	@brief Concatenate two trajectory objects
 *
 * 	When adding A + B, if the final state of A and initial state
 *	of B are the same, this algorithm will skip the initial state
 *	of B in the concatenation to avoid duplicating a state. This
 *	method also overrides the base class behavior and forces time to be
 *	continuous along the concatentated trajectory regardless of whether
 *	the final state of A and in itial state of B are the same
 *
 *	@param rhs the right-hand-side of the addition operation
 *	@return a reference to the concatenated arcset object
 */
Traj& Traj_2bp::operator +=(const Traj &rhs){
	// Create a copy of rhs (it is const)
	Traj temp(rhs);

	// Shift the time in temp by the final time in this trajectory
	double tf = getTimeByIx(-1);
	for(int s = 0; s < temp.getNumNodes(); s++){
		double t = tf + temp.getTimeByIx(s);
		temp.setTimeByIx(s, t);
	}

	// throw Exception("Traj_2bp::operator +=: Not currently implemented!");
	Traj::operator +=(temp);

	return *this;
}//====================================================

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Initialize the extra param vector for info specific to this trajectory
 */
void Traj_2bp::initExtraParam(){
	// Nothing to do here, for now
}//====================================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void Traj_2bp::saveToMat(const char* filename) const{
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
		saveSTMs(matfp);
		pSysData->saveToMat(matfp);
	}

	Mat_Close(matfp);
}//========================================

/**
 *  @brief Populate data in this trajectory from a matlab file
 * 
 *  @param filepath the path to the matlab data file
 *  @throws Exception if the data file cannot be loaded
 */
void Traj_2bp::readFromMat(const char *filepath){
	Traj::readFromMat(filepath);

	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw Exception("Traj: Could not load data from file");
	}

	Mat_Close(matfp);
}//====================================================



}// END of Astrohelion namespace
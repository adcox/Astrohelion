/**
 *	@file tpat_bcr4bpr_traj.cpp
 *	
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

#include "tpat_bcr4bpr_traj.hpp"

#include "tpat_utilities.hpp"

//-----------------------------------------------------
// 		Constructor Functions
//-----------------------------------------------------

/**
 *	@brief Construct a basic BCR4BPR trajectory object
 */
tpat_bcr4bpr_traj::tpat_bcr4bpr_traj() : tpat_trajectory() {
	// dqdT.assign(6,0);
	dqdT.clear();
}

/**
 *	@brief Construct a BCR4BPR trajectory object for the specified system
 *	@param data a system data object describing the BCR4BPR system
 */
tpat_bcr4bpr_traj::tpat_bcr4bpr_traj(tpat_bcr4bpr_sys_data data){
	// dqdT.assign(6,0);
	dqdT.clear();
	sysData = data;
}//====================================================

/**
 *	@brief Construct a BCR4BPR trajectory object with room for a specified number of states
 *	@param n the number of states this trajectory will contain
 */
tpat_bcr4bpr_traj::tpat_bcr4bpr_traj(int n) : tpat_trajectory(n){
	// dqdT.assign(n*6,0);
	dqdT.reserve(n*6);
}//====================================================

/**
 *	@brief Copy the specified trajectory
 *	@param t a BCR4BPR trajectory object
 */
tpat_bcr4bpr_traj::tpat_bcr4bpr_traj(const tpat_bcr4bpr_traj &t) : tpat_trajectory(t){
	copyMe(t);
}//====================================================

//-----------------------------------------------------
// 		Operator Functions
//-----------------------------------------------------

/**
 *	@brief Copy operator; copy a trajectory object into this one.
 *	@param t a trajectory object
 *	@return this trajectory object
 */
tpat_bcr4bpr_traj& tpat_bcr4bpr_traj::operator= (const tpat_bcr4bpr_traj& t){
	tpat_trajectory::operator= (t);
	copyMe(t);
	return *this;
}//====================================================

void tpat_bcr4bpr_traj::copyMe(const tpat_bcr4bpr_traj &t){
	sysData = t.sysData;
	dqdT = t.dqdT;
}//====================================================

/**
 *	@brief Sum two BCR4BPR trajectories
 *
 *	Both trajectories must be propagated in the same system, or else an error will be thrown.
 *	The states from both trajectories are appended [lhs, rhs] without any modification, as is 
 *	the vector of time values because the system is non-autonomous.
 *  Only the STMs from <tt>lhs</tt> are copied to the new trajectory since the 
 *	<tt>rhs</tt> STMs will have no meaning in the new, combined trajectory. The same is true
 *	for the dqdT vector
 *
 *	@param lhs a trajectory
 *	@param rhs a trajectory
 *	@return a new trajectory
 */
tpat_bcr4bpr_traj operator +(const tpat_bcr4bpr_traj &lhs, const tpat_bcr4bpr_traj &rhs){
	if(lhs.sysData != rhs.sysData){
		throw tpat_exception("Cannot sum two BCR4BPR trajectories from different systems!");
	}

	// create a new trajectory object with space for both sets of data to be combined
	tpat_bcr4bpr_traj newTraj(lhs.numPoints + rhs.numPoints);
	
	// Copy the states and times from the LHS into the new guy
	std::copy(lhs.state.begin(), lhs.state.end(), newTraj.getState()->begin());
	std::copy(lhs.times.begin(), lhs.times.end(), newTraj.getTime()->begin());
	
	// Append the rhs state and timeto the end of the new guy's vectors; don't adjust the
	// time because the system is non-autonomous
	std::copy(rhs.state.begin(), rhs.state.end(), newTraj.getState()->begin() + lhs.numPoints);
	std::copy(rhs.times.begin(), rhs.times.end(), newTraj.getTime()->begin() + lhs.numPoints);

	// Copy only the lhs STMs, because there is no gaurantee the two summed trajectories
	// will be continuous and smooth in time and state; same for dqdT
	std::copy(lhs.allSTM.begin(), lhs.allSTM.end(), newTraj.getSTM()->begin());
	std::copy(lhs.dqdT.begin(), lhs.dqdT.end(), newTraj.get_dqdT()->begin());
	
	return newTraj;
}//========================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@return the angle between the P1/P2 line and the inertial x-axis, radians
 */
double tpat_bcr4bpr_traj::getTheta0(){ return sysData.getTheta0(); }

/**
 *	@return the angle between the P2/P3 line (projected into the inertial XY plane)
 *	and the inertial x-axis, radians
 */
double tpat_bcr4bpr_traj::getPhi0(){ return sysData.getPhi0(); }

/**
 *	@return the inclination of the P2/P3 orbital plane relative to the P1/P2 orbital
 *	plane, radians
 */
double tpat_bcr4bpr_traj::getGamma(){ return sysData.getGamma(); }

/**
 *	@return the system data object describing this system
 */
tpat_bcr4bpr_sys_data tpat_bcr4bpr_traj::getSysData(){ return sysData; }

/**
 *	@brief Retrieve a pointer to the dqdT array for in-place editing.
 *	
 *	@return a pointer to the vector of dqdT values;
 */
std::vector<double>* tpat_bcr4bpr_traj::get_dqdT(){ return &dqdT; }

/**
 *	@param i the index of the dqdT vector to retrieve
 *	@return the i'th 6-element dqdT vector. If i is negative, the count
 *	will proceed from the end of the vector, i.e. -1 will return the final time, 
 *	-2 will give the second to last value, etc.
 */
std::vector<double> tpat_bcr4bpr_traj::get_dqdT(int i){
	if(i < 0)
		i += dqdT.size()/6;
	
	std::vector<double> temp(dqdT.begin()+i*6, dqdT.begin()+(i+1)*6);
	return temp;
}//===============================================

tpat_sys_data::system_t tpat_bcr4bpr_traj::getType() const{
	return sysData.getType();
}

/**
 *	@brief Set the system data object
 *	@param data a data object describing the BCR4BP
 */
void tpat_bcr4bpr_traj::setSysData(tpat_bcr4bpr_sys_data data){ sysData = data; }

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	@brief Calls the basic trajectory setLength() method and implements extra catches
 *	specific to the BCR4BPR trajectory object
 */
void tpat_bcr4bpr_traj::setLength(){
	tpat_trajectory::setLength();

	if(dqdT.size()/6 != times.size()){
		printErr("Warning: dqdT vector has different length than time vector!\n");
	}
}//================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_bcr4bpr_traj::saveToMat(const char* filename){
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
		save_dqdT(matfp);
		sysData.saveToMat(matfp);
	}

	Mat_Close(matfp);
}//========================================

/**
 *	@brief Save the Jacobi vector to a file
 * 	@param matFile a pointer to the destination matlab file 
 */
void tpat_bcr4bpr_traj::save_dqdT(mat_t *matFile){
	size_t dims[2] = {static_cast<size_t>(numPoints), 6};
	matvar_t *matvar = Mat_VarCreate("dqdT", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(dqdT[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "dqdT", MAT_COMPRESSION_NONE);
}//=================================================
/**
 *	@file tpat_traj_bcr4bpr.cpp
 *	@brief Derivative of tpat_traj, specific to BCR4BPR
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

#include "tpat_traj_bcr4bpr.hpp"

#include "tpat_utilities.hpp"

//-----------------------------------------------------
// 		Constructor Functions
//-----------------------------------------------------

/**
 *	@brief Construct a basic BCR4BPR trajectory object
 */
tpat_traj_bcr4bpr::tpat_traj_bcr4bpr() : tpat_traj() {
	initExtraParam();
}

/**
 *	@brief Construct a BCR4BPR trajectory object for the specified system
 *	@param data a system data object describing the BCR4BPR system
 */
tpat_traj_bcr4bpr::tpat_traj_bcr4bpr(tpat_sys_data_bcr4bpr data){
	initExtraParam();
	sysData = data;
}//====================================================

/**
 *	@brief Construct a BCR4BPR trajectory object with room for a specified number of states
 *	@param n the number of states this trajectory will contain
 */
tpat_traj_bcr4bpr::tpat_traj_bcr4bpr(int n) : tpat_traj(n){
	initExtraParam();
	extraParam.at(0).reserve(n*6);
}//====================================================

/**
 *	@brief Copy the specified trajectory
 *	@param t a BCR4BPR trajectory object
 */
tpat_traj_bcr4bpr::tpat_traj_bcr4bpr(const tpat_traj_bcr4bpr &t) : tpat_traj(t){
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
tpat_traj_bcr4bpr& tpat_traj_bcr4bpr::operator= (const tpat_traj_bcr4bpr& t){
	tpat_traj::operator= (t);
	copyMe(t);
	return *this;
}//====================================================

/**
 *	@brief Copy the trajectory
 *	@param t a trajectory reference
 */
void tpat_traj_bcr4bpr::copyMe(const tpat_traj_bcr4bpr &t){
	sysData = t.sysData;
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
tpat_traj_bcr4bpr operator +(const tpat_traj_bcr4bpr &lhs, const tpat_traj_bcr4bpr &rhs){
	if(lhs.sysData != rhs.sysData){
		throw tpat_exception("Cannot sum two BCR4BPR trajectories from different systems!");
	}

	tpat_traj_bcr4bpr newTraj(lhs.sysData);
	basicConcat(&lhs, &rhs, &newTraj);

	return newTraj;
}//========================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@return the angle between the P1/P2 line and the inertial x-axis, radians
 */
double tpat_traj_bcr4bpr::getTheta0(){ return sysData.getTheta0(); }

/**
 *	@return the angle between the P2/P3 line (projected into the inertial XY plane)
 *	and the inertial x-axis, radians
 */
double tpat_traj_bcr4bpr::getPhi0(){ return sysData.getPhi0(); }

/**
 *	@return the inclination of the P2/P3 orbital plane relative to the P1/P2 orbital
 *	plane, radians
 */
double tpat_traj_bcr4bpr::getGamma(){ return sysData.getGamma(); }

/**
 *	@return the system data object describing this system
 */
tpat_sys_data_bcr4bpr tpat_traj_bcr4bpr::getSysData(){ return sysData; }

/**
 *	@brief Retrieve a pointer to the dqdT array for in-place editing.
 *	
 *	@return a pointer to the vector of dqdT values;
 */
std::vector<double>* tpat_traj_bcr4bpr::get_dqdT(){ return &(extraParam.at(0)); }

/**
 *	@param i the index of the dqdT vector to retrieve
 *	@return the i'th 6-element dqdT vector. If i is negative, the count
 *	will proceed from the end of the vector, i.e. -1 will return the final time, 
 *	-2 will give the second to last value, etc.
 */
std::vector<double> tpat_traj_bcr4bpr::get_dqdT(int i){
	if(i < 0)
		i += extraParam.at(0).size()/6;
	
	std::vector<double> temp(extraParam.at(0).begin()+i*6, extraParam.at(0).begin()+(i+1)*6);
	return temp;
}//===============================================

tpat_sys_data::system_t tpat_traj_bcr4bpr::getType() const{
	return sysData.getType();
}//===============================================

/**
 *	@brief Retrieve a pointer to the system data object
 *	@return a pointer to the system data object
 */
tpat_sys_data* tpat_traj_bcr4bpr::getSysDataPtr(){ return &sysData; }

/**
 *	@brief Set the system data object
 *	@param data a data object describing the BCR4BP
 */
void tpat_traj_bcr4bpr::setSysData(tpat_sys_data_bcr4bpr data){ sysData = data; }

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

void tpat_traj_bcr4bpr::initExtraParam(){
	numExtraParam = 1;	// one vector for dqdT values
	extraParamRowSize.push_back(6);	// each dqdT has six elements
	extraParam.push_back(std::vector<double>(0));
}//===========================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_traj_bcr4bpr::saveToMat(const char* filename){
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
		saveTime(matfp);
		saveSTMs(matfp);
		saveExtraParam(matfp, 0, 6, "dqdT");
		sysData.saveToMat(matfp);
	}

	Mat_Close(matfp);
}//========================================





//
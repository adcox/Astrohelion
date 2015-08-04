/**
 *	@file tpat_traj_cr3bp.cpp
 *	@brief Derivative of tpat_traj, specific to CR3BP
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

#include "tpat_traj_cr3bp.hpp"
#include "tpat_nodeset_cr3bp.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_utilities.hpp"

#include <cstring>

//-----------------------------------------------------
// 		Constructor Functions
//-----------------------------------------------------

/**
 *	@brief Default constructor; calls constructor for super-class tpat_traj
 *	and additionall initializes the jacobi matrix
 */
tpat_traj_cr3bp::tpat_traj_cr3bp() : tpat_traj(){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a CR3BP trajectory object for the specified system
 *	@param data a system data object describing the system
 */
tpat_traj_cr3bp::tpat_traj_cr3bp(tpat_sys_data_cr3bp data){
	initExtraParam();
	sysData = data;
}//====================================================

/**
 *	@brief Initialize all vectors to have size n; fill each vector with zeros.
 */
tpat_traj_cr3bp::tpat_traj_cr3bp(int n) : tpat_traj(n){
	initExtraParam();
	extraParam.at(0).reserve(n);
}//====================================================

/**
 *	@brief Copy the specified trajectory
 *	@param t trajectory
 */
tpat_traj_cr3bp::tpat_traj_cr3bp(const tpat_traj_cr3bp &t) : tpat_traj(t){
	copyMe(t);
}//====================================================

/**
 *	@brief Create a trajectory from a nodeset
 *
 *	This algorithm will concatenate trajectories integrated from each node in 
 *	the nodeset. It does not check to make sure the arcs are continuous; that
 *	is up to you
 *
 *	@param nodes a nodeset
 *	@return a trajectory formed from the integrated nodeset
 */
tpat_traj_cr3bp tpat_traj_cr3bp::fromNodeset(tpat_nodeset_cr3bp nodes){
	tpat_sys_data *sysData = nodes.getSysData();
	tpat_simulation_engine simEngine(sysData);
	tpat_traj_cr3bp totalTraj;

	for(int n = 0; n < nodes.getNumNodes()-1; n++){
		simEngine.runSim(nodes.getNode(n).getPosVelState(), nodes.getTOF(n));

		if(n == 0){
			totalTraj = simEngine.getCR3BP_Traj();
		}else{
			tpat_traj_cr3bp temp = totalTraj + simEngine.getCR3BP_Traj();
			totalTraj = temp;
		}
	}

	return totalTraj;
}//====================================================

//-----------------------------------------------------
// 		Operators
//-----------------------------------------------------

/**
 *	@brief Copy operator; copy a trajectory object into this one.
 *	@param t a trajectory object
 *	@return this trajectory object
 */
tpat_traj_cr3bp& tpat_traj_cr3bp::operator= (const tpat_traj_cr3bp& t){
	tpat_traj::operator= (t);
	copyMe(t);
	return *this;
}//=====================================================

/**
 *	@brief Copy all data objects specific to the CR3BP Trajectory
 *	@param t the source trajectory; copy its attributes into this one
 */
void tpat_traj_cr3bp::copyMe(const tpat_traj_cr3bp &t){
	sysData = t.sysData;
}//=====================================================

/**
 *	@brief Sum two CR3BP trajectories
 *
 *	Both trajectories must be propagated in the same system, or else an error will be thrown.
 *	The states from both trajectories are appended [lhs, rhs] without any modification, as is 
 *	the vector of Jacobi Constant values. The time vectors are appended such that time is 
 *	continuous throughout the trajectory; this doesn't affect the motion because the CR3BP is
 *	autonomous. STMs from <tt>rhs</tt> are multiplied by the final STM of <tt>lhs</tt> to 
 *	create new, accurate STMs. This assumes that there are no discontinuities in time,
 *	space, or velocity between the end of <tt>lhs</tt> and the beginning of <tt>rhs</tt>
 *
 *	@param lhs a trajectory
 *	@param rhs a trajectory
 *	@return a new trajectory
 */
tpat_traj_cr3bp operator +(const tpat_traj_cr3bp &lhs, const tpat_traj_cr3bp &rhs){
	if(lhs.sysData != rhs.sysData){
		throw tpat_exception("Cannot sum two CR3BP trajectories from different systems");
	}

	tpat_traj_cr3bp newTraj(lhs.sysData);
	basicConcat(&lhs, &rhs, &newTraj);

	return newTraj;
}//========================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Get the value of Jacobi constant at a particular point on the trajectory
 *	@param n the index of the point. If n is negative, the count
 *	will proceed from the end of the vector, i.e. -1 will return the final time, 
 *	-2 will give the second to last value, etc.
 *	@return the Jacobi constant (non-dimensional)
 */
double tpat_traj_cr3bp::getJacobi(int n) const{
	if(n < 0)
		n += extraParam.at(0).size();
	return extraParam.at(0)[n];
}

/**
 *	@brief Retrieve a pointer to the Jacobi array for in-place editing.
 *	
 *	@return a pointer to the vector of Jacobi constants;
 */
std::vector<double>* tpat_traj_cr3bp::getJacobi(){ return &(extraParam.at(0)); }

/**
 *	@brief Set the vector of Jacobi constant values for this trajectory
 *	@param j a vector of Jacobi constants
 */
void tpat_traj_cr3bp::setJacobi(std::vector<double> j){ extraParam.at(0) = j; }

/**
 *	@brief Retrieve data about the system this trajectory was propagated in
 *	@return the system data object
 */
tpat_sys_data_cr3bp tpat_traj_cr3bp::getSysData() const { return sysData; }


/**
 *	@brief Set the system data for this trajectory
 *	@param d a system data object
 */
void tpat_traj_cr3bp::setSysData(tpat_sys_data_cr3bp d){ sysData = d; }

tpat_sys_data::system_t tpat_traj_cr3bp::getType() const{
	return sysData.getType();
}//=====================================================

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

void tpat_traj_cr3bp::initExtraParam(){
	numExtraParam = 1;	// Jacobi
	extraParamRowSize.push_back(1);	// each jacobi value has one element
	extraParam.push_back(std::vector<double>(0));
}//=========================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_traj_cr3bp::saveToMat(const char* filename){
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
		saveExtraParam(matfp, 0, 1, "Jacobi");
		sysData.saveToMat(matfp);
	}

	Mat_Close(matfp);
}//========================================





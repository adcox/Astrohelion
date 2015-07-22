/**
 *	@file tpat_traj_cr3bp.cpp
 *
 *	CR3BP Trajectory
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
	numExtraParam = 1;
}//====================================================

/**
 *	@brief Create a CR3BP trajectory object for the specified system
 *	@param data a system data object describing the system
 */
tpat_traj_cr3bp::tpat_traj_cr3bp(tpat_sys_data_cr3bp data){
	numExtraParam = 1;
	sysData = data;
}//====================================================

/**
 *	@brief Initialize all vectors to have size n; fill each vector with zeros.
 */
tpat_traj_cr3bp::tpat_traj_cr3bp(int n) : tpat_traj(n){
	numExtraParam = 1;
	extraParam.reserve(n);
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
	tpat_simulation_engine simEngine(*sysData);
	tpat_traj_cr3bp totalTraj;

	for(int n = 0; n < nodes.getNumNodes()-1; n++){
		simEngine.runSim(nodes.getNode(n), nodes.getTOF(n));

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

	int skipShift = 1;
	double t1 = lhs.getTol();
	double t2 = rhs.getTol();
	double tol = t1 > t2 ? t1 : t2;
	if(tol == 0)
		tol = 1e-9;

	if(tpat_util::aboutEquals(lhs.getState(-1), rhs.getState(0), 100*tol)){
		skipShift = 1;
	}

	// create a new trajectory object with space for both sets of data to be combined
	tpat_traj_cr3bp newTraj(lhs.sysData);

	// Set the tolerance to the greater of the two
	newTraj.setTol(t1 > t2 ? t1 : t2);

	// Copy the states and times from the LHS into the new guy)
	newTraj.getState()->insert(newTraj.getState()->begin(), lhs.state.begin(), lhs.state.end());
	newTraj.getTime()->insert(newTraj.getTime()->begin(), lhs.times.begin(), lhs.times.end());
	newTraj.getAccel()->insert(newTraj.getAccel()->begin(), lhs.accel.begin(), lhs.accel.end());
	newTraj.getExtraParam()->insert(newTraj.getExtraParam()->begin(), lhs.extraParam.begin(), lhs.extraParam.end());

	// Append the rhs state to the end of the new guy's state vector
	newTraj.getState()->insert(newTraj.getState()->end(), rhs.state.begin()+skipShift*tpat_traj_cr3bp::STATE_SIZE, rhs.state.end());
	newTraj.getAccel()->insert(newTraj.getAccel()->end(), rhs.accel.begin()+skipShift*tpat_traj_cr3bp::ACCEL_SIZE, rhs.accel.end());
	newTraj.getExtraParam()->insert(newTraj.getExtraParam()->end(), rhs.extraParam.begin() +skipShift*lhs.numExtraParam, rhs.extraParam.end());

	// Append the rhs times, adjusted for continuity, to the new guy's time vector; adjustments
	// don't affect the result because system is autonomous
	for (int n = skipShift; n < rhs.numPoints; n++){
		newTraj.getTime()->push_back(lhs.times.back() + rhs.times.at(n) - rhs.times.at(0));
	}

	// Copy the lhs stm
	newTraj.getSTM()->insert(newTraj.getSTM()->begin(), lhs.allSTM.begin(), lhs.allSTM.end());
	// Assume the two are continuous, use matrix multiplication to shift the rhs STMs to be continuous
	for(size_t i = skipShift; i < rhs.allSTM.size(); i++){
		tpat_matrix shiftedSTM = rhs.getSTM(i)*lhs.getSTM(-1);
		newTraj.getSTM()->push_back(shiftedSTM);
	}

	newTraj.setLength();

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
double tpat_traj_cr3bp::getJC(int n) const{
	if(n < 0)
		n += extraParam.size();
	return extraParam[n];
}

/**
 *	@brief Retrieve a pointer to the Jacobi array for in-place editing.
 *	
 *	@return a pointer to the vector of Jacobi constants;
 */
std::vector<double>* tpat_traj_cr3bp::getJacobi(){ return &extraParam; }

void tpat_traj_cr3bp::setLength() { tpat_traj::setLength(); }

/**
 *	@brief Set the vector of Jacobi constant values for this trajectory
 *	@param j a vector of Jacobi constants
 */
void tpat_traj_cr3bp::setJacobi(std::vector<double> j){ extraParam = j; }

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





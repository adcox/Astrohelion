/**
 *  @file tpat_traj_bcr4bp.cpp
 *	@brief Derivative of tpat_traj, specific to BCR4BPR
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

#include "tpat.hpp"

#include "tpat_traj_bcr4bp.hpp"

#include "tpat_exceptions.hpp"
#include "tpat_nodeset_bcr4bp.hpp"
#include "tpat_simulation_engine.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_utilities.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Create a trajectory for a specific system
 *	@param sys a pointer to a system data object
 */
tpat_traj_bcr4bp::tpat_traj_bcr4bp(const tpat_sys_data_bcr4bpr *sys) : tpat_traj(sys){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from another trajectory
 *	@param t a trajectory reference
 */
tpat_traj_bcr4bp::tpat_traj_bcr4bp(const tpat_traj_bcr4bp &t) : tpat_traj(t){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from its base class
 *	@param a an arc data reference
 */
tpat_traj_bcr4bp::tpat_traj_bcr4bp(const tpat_arc_data &a) : tpat_traj(a){
	initExtraParam();
}//====================================================

/**
 *	@brief Create a trajectory from a nodeset
 *
 *	This algorithm will concatenate trajectories integrated from each node in 
 *	the nodeset. It does not check to make sure the arcs are continuous; that
 *	is up to you. The trajectory is constructed via a simulation engine that ignores
 *	crashes as we assume the initial nodeset has been propagated to either ignore
 *	or avoid the primaries; will not challenge that behavior. Each node is integrated
 *	for the associated time-of-flight and added (via operator +()) to a trajectory object.
 *
 *	@param nodes a nodeset
 *	@return a trajectory formed from the integrated nodeset
 *	
 *	@see tpat_traj::operator +()
 */
tpat_traj_bcr4bp tpat_traj_bcr4bp::fromNodeset(tpat_nodeset_bcr4bp nodes){
	const tpat_sys_data_bcr4bpr *sys = static_cast<const tpat_sys_data_bcr4bpr*>(nodes.getSysData());
	tpat_simulation_engine simEngine(sys);
	simEngine.clearEvents();	// don't trigger crashes; assume this has been taken care of already
	tpat_traj_bcr4bp totalTraj(sys);

	for(int n = 0; n < nodes.getNumNodes()-1; n++){
		simEngine.setRevTime(nodes.getTOF(n) < 0);
		simEngine.runSim(nodes.getState(n), nodes.getEpoch(n), nodes.getTOF(n));

		if(n == 0){
			totalTraj = simEngine.getBCR4BPR_Traj();
		}else{
			tpat_traj_bcr4bp temp = simEngine.getBCR4BPR_Traj();
			totalTraj += temp;
		}
	}

	return totalTraj;
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
double tpat_traj_bcr4bp::getTheta0(){
	const tpat_sys_data_bcr4bpr *bcSys = static_cast<const tpat_sys_data_bcr4bpr *>(sysData);
	return bcSys->getTheta0();
}//====================================================

/**
 *	@return the angle between the P2/P3 line (projected into the inertial XY plane)
 *	and the inertial x-axis, radians
 */
double tpat_traj_bcr4bp::getPhi0(){
	const tpat_sys_data_bcr4bpr *bcSys = static_cast<const tpat_sys_data_bcr4bpr *>(sysData);
	return bcSys->getPhi0();
}//====================================================

/**
 *	@return the inclination of the P2/P3 orbital plane relative to the P1/P2 orbital
 *	plane, radians
 */
double tpat_traj_bcr4bp::getGamma(){
	const tpat_sys_data_bcr4bpr *bcSys = static_cast<const tpat_sys_data_bcr4bpr *>(sysData);
	return bcSys->getGamma();
}//====================================================

/**
 *	@param ix the index of the dqdT vector to retrieve
 *	@return the i'th 6-element dqdT vector. If ix is negative, the count
 *	will proceed from the end of the vector, i.e. -1 will return the final time, 
 *	-2 will give the second to last value, etc.
 */
std::vector<double> tpat_traj_bcr4bp::get_dqdT(int ix){
	if(ix < 0)
		ix += steps.size();

	if(ix < 0 || ix > ((int)steps.size()))
		throw tpat_exception("tpat_traj_bcr4bp::getdqdT: invalid index");

	return getExtraParam(ix, 1);
}//====================================================

/**
 *	@brief Set the value of the dqdT vector for the specified step
 *	@param ix the index of the step; if < 0, it will count backwards from the end
 *	@param dqdT a pointer to the dqdT vector; this MUST have at least 6 elements,
 *	or the function will read unallocated memory.
 */
void tpat_traj_bcr4bp::set_dqdT(int ix, const double *dqdT){
	if(ix < 0)
		ix += steps.size();

	if(ix < 0 || ix > ((int)steps.size()))
		throw tpat_exception("tpat_traj_bcr4bp::setdqdT: invalid index");

	for(int i = 0; i < 6; i++)
		steps[ix].setExtraParam(1+i, dqdT[i]);
}//====================================================

/**
 *	@brief Set the value of the dqdT vector for the specified step
 *	@param ix the index of the step; if < 0, it will count backwards from the end
 *	@param dqdT a vector (6 elements) representing the dqdT vector
 */
void tpat_traj_bcr4bp::set_dqdT(int ix, std::vector<double> dqdT){
	if(dqdT.size() != 6)
		throw tpat_exception("tpat_traj_bcr4bp::set_dqdT: Cannot accept a dqdT with anything other than 6 elements");

	set_dqdT(ix, &(dqdT[0]));
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Initialize the extra param vector for info specific to this trajectory
 */
void tpat_traj_bcr4bp::initExtraParam(){
	// This function in tpat_traj was already called, so 
	// numExtraParam has been set to 1 and a row size has
	// been appended for the time variable

	// Add another variable for dqdT
	numExtraParam = 2;
	extraParamRowSize.push_back(6);
}//====================================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_traj_bcr4bp::saveToMat(const char* filename) const{
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
		saveExtraParam(matfp, 1, "dqdT");
		sysData->saveToMat(matfp);
	}

	Mat_Close(matfp);
}//========================================

/**
 *  @brief Populate data in this nodeset from a matlab file
 * 
 *  @param filepath the path to the matlab data file
 */
void tpat_traj_bcr4bp::readFromMat(const char *filepath){
	tpat_traj::readFromMat(filepath);
	
	// Load the matlab file
	mat_t *matfp = Mat_Open(filepath, MAT_ACC_RDONLY);
	if(NULL == matfp){
		throw tpat_exception("tpat_traj_bcr4bp: Could not load data from file");
	}

	readExtraParamFromMat(matfp, 1, "dqdT");
	
	Mat_Close(matfp);
}//====================================================
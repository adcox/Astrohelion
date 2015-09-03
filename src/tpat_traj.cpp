/**
 *  @file tpat_traj.cpp
 *	@brief Stores information about a trajectory
 *
 *	@author Andrew Cox
 *	@version August 30, 2015
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

#include "tpat_traj.hpp"

#include "tpat_traj_step.hpp"
#include "tpat_utilities.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

tpat_traj::tpat_traj(tpat_sys_data *data) : tpat_arc_data(data) {
	initExtraParam();
}//====================================================

tpat_traj::tpat_traj(const tpat_traj &t) : tpat_arc_data(t) {
	initExtraParam();
}//====================================================

tpat_traj::tpat_traj(const tpat_arc_data &a) : tpat_arc_data(a) {
	initExtraParam();
}

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

double tpat_traj::getTime(int ix) const {
	if(ix < 0)
		ix += steps.size();
	tpat_traj_step step(steps[ix]);
	return step.getTime();
}//====================================================

tpat_traj_step tpat_traj::getStep(int ix) const{
	if(ix < 0)
		ix += steps.size();
	return tpat_traj_step(steps[ix]);
}//====================================================

void tpat_traj::appendStep(tpat_traj_step s){ steps.push_back(s); }

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_traj::saveToMat(const char* filename){
	/*	Create a new Matlab MAT file with the given name and optional
	 *	header string. If no header string is given, the default string 
	 *	used containing the software, version, and date in it. If a header
	 *	string is specified, at most the first 116 characters are written to
	 *	the file. Arguments are:
	 *	const char *matname 	- 	the name of the file
	 *	const char *hdr_str 	- 	the 116 byte header string
	 *	enum mat_ft 			- 	matlab file version MAT_FT_MAT5 or MAT_FT_MAT4
	 */
	mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		printErr("tpat_traj::saveToMat: Error creating MAT file\n");
	}else{
		saveState(matfp);
		saveAccel(matfp);
		saveTime(matfp);
		saveSTMs(matfp);
	}

	Mat_Close(matfp);
}//====================================================

/**
 *	@brief Print a useful message describing this trajectory
 */
void tpat_traj::print() const {
	printf("This is a trajectory\n\tTODO - MAKE THIS MESSAGE USEFUL\n");
}//====================================================

/**
 *	@brief Save the times at each integration step
 *	@param file a pointer to an open Mat file
 */
void tpat_traj::saveTime(mat_t *file){
	saveExtraParam(file, 0, "Time");
}//====================================================

void tpat_traj::initExtraParam(){
	// ExtraParam = [time]
	numExtraParam = 1;
	extraParamRowSize.push_back(1);
}//====================================================
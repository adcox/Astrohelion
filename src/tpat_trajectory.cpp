/**
 *	@file tpat_trajectory.cpp
 *
 *	Trajectory Class: contains info about a trajectory
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
 
#include "tpat_trajectory.hpp"

#include "tpat_matrix.hpp"
#include "tpat_sys_data.hpp"
#include "tpat_utilities.hpp"
 
#include <cstdlib>

//-----------------------------------------------------
// 		Constructors
//-----------------------------------------------------

/**
 *	@brief Default constructor. 
 *
 *	Sets the number of points to 0
 */
tpat_trajectory::tpat_trajectory(){
	numPoints = 0;
	state.clear();
	times.clear();
	allSTM.clear();
}

/**
 *	@brief Create a trajectory to hold a specific number of points. 
 *
 *	The state, time, and STM
 *	vectors are all initialized with enough space to hold the number of points
 *	specified. State and Time are filled with zeros, the STMs are filled with identity
 *	matrices. If initializing the trajectory this way (recommended), use the in-place
 *	set() functions to manipulate the allocated state, time, and STM vectors using
 *	their pointers.
 *
 *	@param n the number of points
 */
tpat_trajectory::tpat_trajectory(int n){
	numPoints = n;
	state.reserve(n*STATE_WIDTH);
	times.reserve(n);
	allSTM.reserve(n);
}

/**
 *	@brief Create a copy of the input trajectory
 *	@param t a trajectory
 */
tpat_trajectory::tpat_trajectory(const tpat_trajectory& t){
	numPoints = t.numPoints;
	state = t.state;
	times = t.times;
	allSTM = t.allSTM;
}//=========================================

//-----------------------------------------------------
// 		Operators
//-----------------------------------------------------

/**
 *	@brief Assignment operator. Directly copy the vectors from one trajectory into this one
 */
tpat_trajectory& tpat_trajectory::operator= (const tpat_trajectory& t){
	numPoints = t.numPoints;
	state = t.state;
	times = t.times;
	allSTM = t.allSTM;

	return *this;
}//============================================


//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Retrieve the number of states, times, STMs, etc. on this trajectory
 *	@return the number of points along this trajectory
 */
int tpat_trajectory::getLength() const{ return numPoints; }

/**
 *	@return a pointer to the beginning of the state array.
 */
std::vector<double>* tpat_trajectory::getState(){ return &state;}

/**
 *	@brief Retrieve a single state along the trajectory
 *
 *	@param n the index of the state (0 is the first state, or IC). If n is negative,
 *	the countwill proceed from the end of the vector, i.e. -1 will return the 
 *	final time, -2 will give the second to last value, etc.
 *	@return a vector representing the full state (pos, vel, accel)
 */
std::vector<double> tpat_trajectory::getState(int n) const{
	if(n < 0)
		n += state.size()/STATE_WIDTH;

	std::vector<double>::const_iterator first = state.begin() + n*STATE_WIDTH;
	std::vector<double>::const_iterator last = state.begin() + (n+1)*STATE_WIDTH;
	std::vector<double> oneState(first, last);
	return oneState;
}

/**
 *	@brief Retrieve the time at a specific point along the trajectory
 *	@param n the index of the point (starts at 0). If n is negative, the count
 *	will proceed from the end of the vector, i.e. -1 will return the final time, 
 *	-2 will give the second to last value, etc.
 *	@return the non-dimensional time along the trajectory
 */
double tpat_trajectory::getTime(int n) const {
	if(n < 0)
		n += times.size();
	return times[n];
}

/**
 *	@return a pointer to the begining of the time vector
 */
std::vector<double>* tpat_trajectory::getTime(){ return &times; }

/**
 *	@brief Retrieve the STM at a specific point along the trajectory
 *	@param n the index of the point (starts at 0). If n is negative, the count
 *	will proceed from the end of the vector, i.e. -1 will return the final time, 
 *	-2 will give the second to last value, etc.
 *	@return the STM
 */
tpat_matrix tpat_trajectory::getSTM(int n) const {
	if(n < 0)
		n += allSTM.size();
	tpat_matrix temp(allSTM[n]);
	return temp;
}

/**
 *	@brief Useful for setting the STM (in place) after the trajectory has been initialized
 *	@return a pointer to the beginning of the vector of STMs
 */
std::vector<tpat_matrix>* tpat_trajectory::getSTM(){ return &allSTM;}

/**
 *	@return the system data type
 */
tpat_sys_data::system_t tpat_trajectory::getType() const { return tpat_sys_data::UNDEF_SYS; }

/** 
 *	@brief Set the state vector by copying a vector
 *	@param s a vector of non-dimensional states. The vector should be 1D in Row-Major
 *	format; all units are non-dimensional
 */
void tpat_trajectory::setState(std::vector<double> s){ state = s; }

/** 
 *	@brief Set the time vector by copying a vector
 *	@param t a vector of non-dimensional times along the trajectory. Length
 *	must match the length of the state vector
 */
void tpat_trajectory::setTime(std::vector<double> t){ times = t; }

/** 
 *	@brief Set the STM vector by copying a vector
 *	@param phi a vector of STMs, one for every point along the trajectory.
 */
void tpat_trajectory::setSTMs(std::vector<tpat_matrix> phi){ allSTM = phi; }

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	@brief Set the number of points by checking the number of data in 
 *	the state, time, and STM vectors.
 */
void tpat_trajectory::setLength(){
	int sL = state.size()/STATE_WIDTH;	// row-major format, 9 elements per row
	int tL = times.size();
	int pL = allSTM.size();

	if(sL == tL && tL == pL){
		numPoints = sL;
	}else{
		printWarn("Trajectory has vectors with different lengths\n");
		numPoints = sL;
	}
}//=======================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_trajectory::saveToMat(const char* filename){
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
	}

	Mat_Close(matfp);
}//========================================

/**
 *	@brief Save the state vector [pos, vel, accel] to a file
 *	@param matFile a pointer to the destination matlab file 
 */
void tpat_trajectory::saveState(mat_t *matFile){

	// We store data in row-major order, but the Matlab file-writing algorithm takes data
	// in column-major order, so we transpose our vector and split it into two smaller ones
	std::vector<double> posVel(numPoints*6);
	std::vector<double> accel(numPoints*3);

	for(int r = 0; r < numPoints; r++){
		for(int c = 0; c < STATE_WIDTH; c++){
			if(c < 6)
				posVel[c*numPoints + r] = state[r*STATE_WIDTH + c];
			else
				accel[(c-6)*numPoints + r] = state[r*STATE_WIDTH + c];
		}
	}

	// Next, create a matlab variable for the state and save it to the file
	/*	Create a matlab variable. Arguments are:
	 *	const char *name 	- "State"
	 *	enum matio_classes 	- MAT_C_DOUBLE, Matlab double-precision variable class
	 *	enum matio_types 	- MAT_T_DOUBLE, Matlab IEEE 754 double precision data type
	 * 	int rank 			- 2 - the variable rank. Must be 2 or more; not really sure what this does
	 *	size_t dims 		- dims - the dimensions of the variable (e.g. matrix size) {rows, cols}
	 *	void *data 			- data - the variable we're saving. The algorithm assumes data is in column-major 
	 *							format
	 *	int opt 			- 0, or bit-wise OR of the following options:
	 *							MAT_F_DONT_COPY_DATA: just use the pointer to the data, don't copy it. 
	 *								Note that the pointer should not be freed until you are done with 
	 *								the matvar. The Mat_VarFree function will NOT free data that was
	 *								created with MAT_F_DONT_COPY_DATA, so free it yourself.
	 *							MAT_F_COMPLEX: specify that the data is complex
	 *							MAT_F_GLOBAL: make the matlab variable global
	 *							MAT_F_LOGICAL: this variable is a logical variable
	 */
	size_t dims[2] = {static_cast<size_t>(numPoints), 6};
	matvar_t *matvar = Mat_VarCreate("State", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(posVel[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "State", MAT_COMPRESSION_NONE);

	// Repeat the procedure with the accelerations
	dims[1] = 3;
	matvar = Mat_VarCreate("Accel", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(accel[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "Accel", MAT_COMPRESSION_NONE);
}//==================================================

/**
 *	@brief Save the time vector to a file
 * 	@param matFile a pointer to the destination matlab file 
 */
void tpat_trajectory::saveTime(mat_t *matFile){
	size_t dims[2] = {static_cast<size_t>(numPoints), 1};
	matvar_t *matvar = Mat_VarCreate("Time", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(times[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "Time", MAT_COMPRESSION_NONE);
}//=================================================

/**
 *	@brief Save the STMs to a file; STMs are stored in a 6x6xn array for 
 *	compatibility with existing MATLAB scripts
 *	@param matFile a pointer to the destination matlab file 
 */
void tpat_trajectory::saveSTMs(mat_t *matFile){
	// Create one large vector to put all the STM elements in
	std::vector<double> allSTMEl(numPoints*36);
	for (int n = 0; n < numPoints; n++){
		// get the transpose of the STM matrix; we need to store it in column-major order
		// and it's currently in row-major order
		tpat_matrix P = trans(allSTM[n]);
		// Retrieve the data from the matrix
		double *matData = P.getDataPtr();
		// Store that data in our huge vector
		std::copy(matData, matData+36, &(allSTMEl[0]) + n*36);
	}

	size_t dims[3] = {6, 6, static_cast<size_t>(numPoints)};
	matvar_t *matvar = Mat_VarCreate("STM", MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims, &(allSTMEl[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "STM", MAT_COMPRESSION_NONE);
}//=========================================





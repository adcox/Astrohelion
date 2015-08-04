/**
 *	@file tpat_traj.cpp
 *	@brief stores information about a trajectory
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
 
#include "tpat_traj.hpp"

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
tpat_traj::tpat_traj(){
	numPoints = 0;
	state.clear();
	times.clear();
	allSTM.clear();
}

/**
 *	@brief Create a trajectory to hold a specific number of points. 
 *
 *	The state, time, and STM are initialized with enough space reserved to hold
 *	the number of points specified. No data is actually present; the data vectors
 *	 are empty, but pre-allocating the space may save computation time.
 *
 *	If initializing the trajectory this way (recommended), use the in-place
 *	set() functions to manipulate the allocated state, time, and STM vectors using
 *	their pointers.
 *
 *	@param n the number of points
 */
tpat_traj::tpat_traj(int n){
	numPoints = n;
	state.reserve(n*9);	// may need more, but this is a good starting spot
	times.reserve(n);
	allSTM.reserve(n);
}//====================================================

/**
 *	@brief Create a copy of the input trajectory
 *	@param t a trajectory reference
 */
tpat_traj::tpat_traj(const tpat_traj& t){
	copyMe(t);
}//====================================================

/**
 *	@brief Destructor
 */
tpat_traj::~tpat_traj(){
	state.clear();
	accel.clear();
	times.clear();
	extraParam.clear();
	allSTM.clear();
}//====================================================

//-----------------------------------------------------
// 		Operators
//-----------------------------------------------------

/**
 *	@brief Assignment operator. Directly copy the vectors from one trajectory into this one
 *	@param t a trajectory reference
 */
tpat_traj& tpat_traj::operator= (const tpat_traj& t){
	copyMe(t);
	return *this;
}//============================================

/**
 *	@brief Copy all the data from one trajectory to another
 *	@param t a trajectory reference
 */
void tpat_traj::copyMe(const tpat_traj &t){
	numPoints = t.numPoints;
	state = t.state;
	accel = t.accel;
	times = t.times;
	extraParam = t.extraParam;
	numExtraParam = t.numExtraParam;
	extraParamRowSize = t.extraParamRowSize;
	allSTM = t.allSTM;
	tol = t.tol;
}//====================================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Get one acceleration vector
 *	@param n the index of the step; if n < 0, it will count backwards from
 *	the end of the trajectory
 *	@return the nth acceleration vector [ax, ay, az]
 */
std::vector<double> tpat_traj::getAccel(int n) const {
	if(n < 0)
		n += accel.size()/ACCEL_SIZE;

	std::vector<double>::const_iterator first = accel.begin() + n*ACCEL_SIZE;
	std::vector<double>::const_iterator last = accel.begin() + (n+1)*ACCEL_SIZE;
	std::vector<double> oneAccel(first, last);
	return oneAccel;
}//======================================================

/**	
 *	@brief Retrieve a pointer to the acceleration vector
 *	@return a pointer to the acceleration vector
 */
std::vector<double>* tpat_traj::getAccel(){ return &accel; }

/**
 *	@brief Retrieve the number of states, times, STMs, etc. on this trajectory
 *	@return the number of points along this trajectory
 */
int tpat_traj::getLength() const{ return numPoints; }

/**
 *	@return a pointer to the beginning of the state array.
 */
std::vector<double>* tpat_traj::getState(){ return &state;}

/**
 *	@brief get a vector of one coordinate's evolution during the flight
 *	@param i the index of the coordinate, i.e. (0 = x, 1 = y, 2 = z, 3 = v_x, 4 = v_y, 5 = v_z)
 */
std::vector<double> tpat_traj::getCoord(int i) const{
	if(i >= STATE_SIZE)
		throw tpat_exception("Coordinate Index Out of Range");

	std::vector<double> coord;
	for(size_t n = 0; n < state.size()/STATE_SIZE; n++){
		coord.push_back(state[n*STATE_SIZE + i]);
	}

	return coord;
}//==========================================================

/**
 *	@brief Retrieve a single set of extra parameters
 *
 *	@param n the index of the set of parameters. If n is negative,
 *	the count will proceed from the end of the vector.
 *	@return a vector representing the values of all extra parameters
 *	at the nth integrated step
 */
std::vector<double> tpat_traj::getExtraParam(int n) const{
	if(n < 0)
		n += extraParam.size()/numExtraParam;

	return extraParam.at(n);
	// std::vector<double>::const_iterator first = extraParam.begin() + n*numExtraParam;
	// std::vector<double>::const_iterator last = extraParam.begin() + (n+1)*numExtraParam;
	// std::vector<double> oneSetParam(first, last);
	// return oneSetParam;
}//==========================================================

/**
 *	@brief Retrieve a pointer to the vector of extra parameters
 *	@param n the index of the vector of parameters
 *	@return a pointer to the vector of extra parameters
 */
std::vector<double>* tpat_traj::getExtraParamPtr(int n){ return &(extraParam.at(n)); }

/**
 *	@brief Retrieve a single state along the trajectory
 *
 *	@param n the index of the state (0 is the first state, or IC). If n is negative,
 *	the count will proceed from the end of the vector, i.e. -1 will return the 
 *	final time, -2 will give the second to last value, etc.
 *	@return a vector representing the full state (pos, vel)
 */
std::vector<double> tpat_traj::getState(int n) const{
	if(n < 0)
		n += state.size()/STATE_SIZE;

	std::vector<double>::const_iterator first = state.begin() + n*STATE_SIZE;
	std::vector<double>::const_iterator last = state.begin() + (n+1)*STATE_SIZE;
	std::vector<double> oneState(first, last);
	return oneState;
}//==========================================================

/**
 *	@brief Retrieve the time at a specific point along the trajectory
 *	@param n the index of the point (starts at 0). If n is negative, the count
 *	will proceed from the end of the vector, i.e. -1 will return the final time, 
 *	-2 will give the second to last value, etc.
 *	@return the non-dimensional time along the trajectory
 */
double tpat_traj::getTime(int n) const {
	if(n < 0)
		n += times.size();
	return times[n];
}//==========================================================

/**
 *	@return a pointer to the begining of the time vector
 */
std::vector<double>* tpat_traj::getTime(){ return &times; }

/**
 *	@brief Retrieve the tolerance used to generate this trajectory
 */
double tpat_traj::getTol() const { return tol; }

/**
 *	@brief Retrieve the STM at a specific point along the trajectory
 *	@param n the index of the point (starts at 0). If n is negative, the count
 *	will proceed from the end of the vector, i.e. -1 will return the final time, 
 *	-2 will give the second to last value, etc.
 *	@return the STM
 */
tpat_matrix tpat_traj::getSTM(int n) const {
	if(n < 0)
		n += allSTM.size();
	tpat_matrix temp(allSTM[n]);
	return temp;
}//==========================================================

/**
 *	@brief Useful for setting the STM (in place) after the trajectory has been initialized
 *	@return a pointer to the beginning of the vector of STMs
 */
std::vector<tpat_matrix>* tpat_traj::getSTM(){ return &allSTM;}

/**
 *	@return the system data type
 */
tpat_sys_data::system_t tpat_traj::getType() const { return tpat_sys_data::UNDEF_SYS; }

/**
 *	@brief Set the acceleration vector
 *	@param a a vector of non-dimensional accelerations. The vector should be 1D in
 *	Row-Major format.
 */
void tpat_traj::setAccel(std::vector<double> a){ accel = a; }

/** 
 *	@brief Set the state vector by copying a vector
 *	@param s a vector of non-dimensional states. The vector should be 1D in Row-Major
 *	format
 */
void tpat_traj::setState(std::vector<double> s){ state = s; }

/** 
 *	@brief Set the time vector by copying a vector
 *	@param t a vector of non-dimensional times along the trajectory. Length
 *	must match the length of the state vector
 */
void tpat_traj::setTime(std::vector<double> t){ times = t; }

/** 
 *	@brief Set the STM vector by copying a vector
 *	@param phi a vector of STMs, one for every point along the trajectory.
 */
void tpat_traj::setSTMs(std::vector<tpat_matrix> phi){ allSTM = phi; }

/**
 *	@brief Set the vector of extra parameters. 
 *	@param ix the index of the extra parameter
 *	@param e a vector of parameters
 */
void tpat_traj::setExtraParam(int ix, std::vector<double> e){ extraParam.at(ix) = e; }

/**
 *	@brief Set the tolerance for this trajectory
 *
 *	This information is useful later when comparing errors between states
 *	@param d the tolerance
 */
void tpat_traj::setTol(double d){ tol = d; }

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	@brief Concatenate two trajectories
 *	@param lhs a pointer to the left-hand side of the sum
 *	@param rhs a pointer to the right-hand side of the sum
 *	@param out a pointer to the output trajectory object
 */
void basicConcat(const tpat_traj *lhs, const tpat_traj *rhs, tpat_traj *out){
	int skipShift = 1;
	double t1 = lhs->getTol();
	double t2 = rhs->getTol();
	double tol = t1 > t2 ? t1 : t2;
	if(tol == 0)
		tol = 1e-9;

	if(tpat_util::aboutEquals(lhs->getState(-1), rhs->getState(0), 100*tol)){
		skipShift = 1;
	}
	
	// Copy all data form LHS into OUT
	out->state.insert(out->state.end(), lhs->state.begin(), lhs->state.end());
	out->accel.insert(out->accel.end(), lhs->accel.begin(), lhs->accel.end());
	out->times.insert(out->times.end(), lhs->times.begin(), lhs->times.end());

	for(size_t i = 0; i < lhs->extraParam.size(); i++){
		out->extraParam.at(i).insert(out->extraParam.at(i).end(), lhs->extraParam.at(i).begin(), lhs->extraParam.at(i).end());
	}

	// Append data from RHS to OUT
	out->state.insert(out->state.end(), rhs->state.begin()+skipShift*rhs->STATE_SIZE, rhs->state.end());
	out->accel.insert(out->accel.end(), rhs->accel.begin()+skipShift*rhs->ACCEL_SIZE, rhs->accel.end());
	out->times.insert(out->times.end(), rhs->times.begin()+skipShift, rhs->times.end());

	for(size_t i = 0; i < lhs->extraParam.size(); i++){
		int rowSize = rhs->extraParamRowSize.at(i);
		out->extraParam.at(i).insert(out->extraParam.at(i).end(), rhs->extraParam.at(i).begin()+skipShift*rowSize, 
			rhs->extraParam.at(i).end());
	}

	// Copy the lhs stm
	out->allSTM.insert(out->allSTM.end(), lhs->allSTM.begin(), lhs->allSTM.end());
	
	// Assume the two are continuous, use matrix multiplication to shift the rhs STMs to be continuous
	for(size_t i = skipShift; i < rhs->allSTM.size(); i++){
		tpat_matrix shiftedSTM = rhs->getSTM(i)*lhs->getSTM(-1);
		out->allSTM.push_back(shiftedSTM);
	}

	out->setLength();
}//========================================================

/**
 *	@brief Set the number of points by checking the number of data in 
 *	the state, time, and STM vectors.
 */
void tpat_traj::setLength(){
	size_t sL = state.size()/STATE_SIZE;
	size_t aL = accel.size()/ACCEL_SIZE;
	size_t tL = times.size();
	size_t pL = allSTM.size();
	bool stillTrue = true;

	if(sL == tL && tL == pL && aL == sL){
		for(size_t i = 0; i < extraParam.size(); i++){
			stillTrue = stillTrue && extraParam.at(i).size()/extraParamRowSize.at(i) == sL;
		}
	}

	numPoints = sL;

	if(!stillTrue)
		printWarn("Trajectory has vectors with different lengths\n");
}//=======================================

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void tpat_traj::saveToMat(const char* filename){
	// TODO: Check for propper file extension, add if necessary

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
		printErr("Error creating MAT file\n");
	}else{
		saveState(matfp);
		saveAccel(matfp);
		saveTime(matfp);
		saveSTMs(matfp);
	}

	Mat_Close(matfp);
}//========================================

/**
 *	@brief Save the state vector [pos, vel, accel] to a file
 *	@param matFile a pointer to the destination matlab file 
 */
void tpat_traj::saveState(mat_t *matFile){

	// We store data in row-major order, but the Matlab file-writing algorithm takes data
	// in column-major order, so we transpose our vector and split it into two smaller ones
	std::vector<double> posVel(state.size());

	for(size_t r = 0; r < state.size()/STATE_SIZE; r++){
		for(int c = 0; c < STATE_SIZE; c++){
			posVel[c*state.size()/STATE_SIZE + r] = state[r*STATE_SIZE + c];
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
	size_t dims[2] = {state.size()/STATE_SIZE, STATE_SIZE};
	matvar_t *matvar = Mat_VarCreate("State", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(posVel[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "State", MAT_COMPRESSION_NONE);
}//==================================================

/**
 *	@brief Save the acceleration vector to file
 *	@param matFile a pointer to the destination mat-file
 */
void tpat_traj::saveAccel(mat_t *matFile){
	// We store data in row-major order, but the Matlab file-writing algorithm takes data
	// in column-major order, so we transpose our vector and split it into two smaller ones
	std::vector<double> accel_colMaj(accel.size());

	for(size_t r = 0; r < accel.size()/ACCEL_SIZE; r++){
		for(int c = 0; c < ACCEL_SIZE; c++){
			accel_colMaj[c*accel.size()/ACCEL_SIZE + r] = state[r*ACCEL_SIZE + c];
		}
	}
	
	size_t dims[2] = {accel.size()/ACCEL_SIZE, ACCEL_SIZE};
	matvar_t *matvar = Mat_VarCreate("Accel", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(accel_colMaj[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "Accel", MAT_COMPRESSION_NONE);
}//=====================================================

/**
 *	@brief Save one of the extra parameters to file
 *	@param matFile a pointer to the destination mat-file
 *	@param varIx the index of the parameter
 *	@param width the number elements each row contains. For example, if I'm storing a 1x3 vector,
 *	varIx would be 0 and width would be three.
 *	@param name the name of the variable being saved
 */
void tpat_traj::saveExtraParam(mat_t *matFile, int varIx, int width, const char *name){
	if(varIx > numExtraParam || varIx < 0)
		throw tpat_exception("Could not save extra parameter; index out of bounds");

	// Get the specified coordinate
	std::vector<double> param;
	for(size_t i = 0; i < extraParam.at(varIx).size()/width; i++){
		param.insert(param.end(), extraParam.at(varIx).begin()+i*width, extraParam.at(varIx).begin()+(i+1)*width);
	}

	size_t dims[2] = {param.size()/width, static_cast<size_t>(width)};
	matvar_t *matvar = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(param[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, name, MAT_COMPRESSION_NONE);
}//======================================================


/**
 *	@brief Save the time vector to a file
 * 	@param matFile a pointer to the destination matlab file 
 */
void tpat_traj::saveTime(mat_t *matFile){
	size_t dims[2] = {static_cast<size_t>(numPoints), 1};
	matvar_t *matvar = Mat_VarCreate("Time", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(times[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "Time", MAT_COMPRESSION_NONE);
}//=================================================

/**
 *	@brief Save the STMs to a file; STMs are stored in a 6x6xn array for 
 *	compatibility with existing MATLAB scripts
 *	@param matFile a pointer to the destination matlab file 
 */
void tpat_traj::saveSTMs(mat_t *matFile){
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





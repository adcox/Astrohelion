/**
 *	Trajectory Class: contains info about a trajectory
 */

#include "adtk_trajectory.hpp"
#include <iostream>
#include <cstdlib>

using namespace std;

//-----------------------------------------------------
// 		Constructors
//-----------------------------------------------------

/**
 *	Default constructor. Sets the number of points to 1 and initializes all vectors
 *	to have one point of data with all values set to zero.
 */
adtk_trajectory::adtk_trajectory(){
	numPoints = 1;
	state.assign(STATE_WIDTH,0);	// assign 9 values of zero
	times.assign(1,0);
	allSTM.assign(1, adtk_matrix::Identity(6));
}

/**
 *	Create a trajectory to hold a specific number of points. The state, time, and STM
 *	vectors are all initialized with enough space to hold the number of points
 *	specified. State and Time are filled with zeros, the STMs are filled with identity
 *	matrices. If initializing the trajectory this way (recommended), use the in-place
 *	set() functions to manipulate the allocated state, time, and STM vectors using
 *	their pointers.
 *
 *	@param n the number of points
 */
adtk_trajectory::adtk_trajectory(int n){
	numPoints = n;
	state.assign(n*STATE_WIDTH, 0);
	times.assign(n, 0);
	allSTM.assign(n, adtk_matrix::Identity(6));
}

/**
 *	Destructor
 */
adtk_trajectory::~adtk_trajectory(){
	state.clear();
	times.clear();
	allSTM.clear();
}//==================================

//-----------------------------------------------------
// 		Operators
//-----------------------------------------------------

/**
 *	Copy operator. Directly copy the vectors from one trajectory into this one
 */
adtk_trajectory& adtk_trajectory::operator= (const adtk_trajectory& t){
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
 *	@return the number of points along this trajectory
 */
int adtk_trajectory::getLength(){ return numPoints; }

/**
 *	@return a pointer to the beginning of the state array.
 */
vector<double>* adtk_trajectory::getState(){ return &state;}

/**
 *	Retrieve a single state along the trajectory
 *
 *	@param n the index of the state (0 is the first state, or IC)
 *	@return a vector representing the full state (pos, vel, accel)
 */
vector<double> adtk_trajectory::getState(int n){
	vector<double>::const_iterator first = state.begin() + n*STATE_WIDTH;
	vector<double>::const_iterator last = state.begin() + (n+1)*STATE_WIDTH - 1;
	vector<double> oneState(first, last);
	return oneState;
}

/**
 *	Retrieve the time at a specific point along the trajectory
 *	@param n the index of the point (starts at 0)
 *	@return the non-dimensional time along the trajectory
 */
double adtk_trajectory::getTime(int n){ return times[n]; }

/**
 *	@return a pointer to the begining of the time vector
 */
vector<double>* adtk_trajectory::getTime(){ return &times; }

/**
 *	Retrieve the STM at a specific point along the trajectory
 *	@param n the index of the point (starts at 0)
 *	@return the STM
 */
adtk_matrix adtk_trajectory::getSTM(int n){ return allSTM[n]; }

/**
 *	Useful for setting the STM (in place) after the trajectory has been initialized
 *	@return a pointer to the beginning of the vector of STMs
 */
vector<adtk_matrix>* adtk_trajectory::getSTM(){ return &allSTM;}

/**
 *	Set the number of points by checking the number of data in 
 *	the state, time, and STM vectors.
 */
void adtk_trajectory::setLength(){
	int sL = state.size()/9;	// row-major format, 9 elements per row
	int tL = times.size();
	int pL = allSTM.size();

	if(sL == tL && tL == pL){
		numPoints = sL;
	}else{
		cout << "Warning: trajectory has vectors with different lengths:" << endl;
		cout << " state: " << sL << "\n time: " << tL << "\n STM: " << pL<< endl;
		numPoints = sL;
	}
}//=======================================

/** 
 *	Set the state vector by copying a vector
 *	@param s a vector of non-dimensional states. The vector should be 1D in Row-Major
 *	format; all units are non-dimensional
 */
void adtk_trajectory::setState(vector<double> s){ state = s; }

/** 
 *	Set the time vector by copying a vector
 *	@param t a vector of non-dimensional times along the trajectory. Length
 *	must match the length of the state vector
 */
void adtk_trajectory::setTime(vector<double> t){ times = t; }

/** 
 *	Set the STM vector by copying a vector
 *	@param phi a vector of STMs, one for every point along the trajectory.
 */
void adtk_trajectory::setSTMs(vector<adtk_matrix> phi){ allSTM = phi; }
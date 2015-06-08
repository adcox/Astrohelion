/**
 *	@file adtk_cr3bp_traj.cpp
 *
 *	CR3BP Trajectory
 */

#include "adtk_cr3bp_traj.hpp"

#include "adtk_utilities.hpp"

#include <cstring>
 
using namespace std;

//-----------------------------------------------------
// 		Constructor Functions
//-----------------------------------------------------

/**
 *	@brief Default constructor; calls constructor for super-class adtk_trajectory
 *	and additionall initializes the jacobi matrix
 */
adtk_cr3bp_traj::adtk_cr3bp_traj() : adtk_trajectory(){
	// jacobi.assign(1,0);
	jacobi.clear();
}

/**
 *	@brief Create a CR3BP trajectory object for the specified system
 *	@param data a system data object describing the system
 */
adtk_cr3bp_traj::adtk_cr3bp_traj(adtk_cr3bp_sys_data data){
	// jacobi.assign(1,0);
	jacobi.clear();
	sysData = data;
}

/**
 *	@brief Initialize all vectors to have size n; fill each vector with zeros.
 */
adtk_cr3bp_traj::adtk_cr3bp_traj(int n) : adtk_trajectory(n){
	// jacobi.assign(n, 0);
	jacobi.reserve(n);
}

/**
 *	@brief Copy the specified trajectory
 *	@param t trajectory
 */
adtk_cr3bp_traj::adtk_cr3bp_traj(const adtk_cr3bp_traj &t) : adtk_trajectory(t){
	sysData = t.sysData;
	jacobi = t.jacobi;
}

//-----------------------------------------------------
// 		Operators
//-----------------------------------------------------

/**
 *	@brief Copy operator; copy a trajectory object into this one.
 *	@param t a trajectory object
 *	@return this trajectory object
 */
adtk_cr3bp_traj& adtk_cr3bp_traj::operator= (const adtk_cr3bp_traj& t){
	adtk_trajectory::operator= (t);
	jacobi = t.jacobi;
	return *this;
}//=====================================

/**
 *	@brief Sum two CR3BP trajectories
 *
 *	Both trajectories must be propagated in the same system, or else an error will be thrown.
 *	The states from both trajectories are appended [lhs, rhs] without any modification, as is 
 *	the vector of Jacobi Constant values. The time vectors are appended such that time is 
 *	continuous throughout the trajectory; this doesn't affect the motion because the CR3BP is
 *	autonomous. Only the STMs from <tt>lhs</tt> are copied to the new trajectory since the 
 *	<tt>rhs</tt> STMs will have no meaning in the new, combined trajectory.
 *
 *	@param lhs a trajectory
 *	@param rhs a trajectory
 *	@return a new trajectory
 */
adtk_cr3bp_traj operator +(const adtk_cr3bp_traj &lhs, const adtk_cr3bp_traj &rhs){
	if(lhs.sysData.getPrimary(0).compare(rhs.sysData.getPrimary(0)) != 0 ||
			lhs.sysData.getPrimary(1).compare(rhs.sysData.getPrimary(1)) != 0){
		printErr("Cannot sum two CR3BP trajectories from different systems!\n");
		throw;
	}

	// create a new trajectory object with space for both sets of data to be combined
	adtk_cr3bp_traj newTraj(lhs.numPoints + rhs.numPoints);
	newTraj.setSysData(lhs.sysData);

	// Copy the states and times from the LHS into the new guy
	copy(lhs.state.begin(), lhs.state.end(), newTraj.getState()->begin());
	copy(lhs.times.begin(), lhs.times.end(), newTraj.getTime()->begin());
	
	// Append the rhs state to the end of the new guy's state vector
	copy(rhs.state.begin(), rhs.state.end(), newTraj.getState()->begin() + lhs.numPoints);

	// Append the rhs times, adjusted for continuity, to the new guy's time vector; adjustments
	// don't affect the result because system is autonomous
	double *newTimes = &(newTraj.getTime()->at(lhs.numPoints));
	for (int n = 0; n < rhs.numPoints; n++){
		*newTimes = lhs.times.back() + rhs.times.at(n) - rhs.times.at(0);
		newTimes++;
	}

	// Copy only the lhs STMs, because there is no gaurantee the two summed trajectories
	// will be continuous and smooth in time and state
	copy(lhs.allSTM.begin(), lhs.allSTM.end(), newTraj.getSTM()->begin());

	// Copy Jacobi constant
	copy(lhs.jacobi.begin(), lhs.jacobi.end(), newTraj.getJC()->begin());
	copy(rhs.jacobi.begin(), rhs.jacobi.end(), newTraj.getJC()->begin() + lhs.numPoints);

	return newTraj;
}//========================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Get the value of Jacobi constant at a particular point on the trajectory
 *	@param n the index of the point
 *	@return the Jacobi constant (non-dimensional)
 */
double adtk_cr3bp_traj::getJC(int n){ return jacobi[n]; }

/**
 *	@brief Retrieve a pointer to the Jacobi array for in-place editing.
 *	
 *	@return a pointer to the vector of Jacobi constants;
 */
vector<double>* adtk_cr3bp_traj::getJC(){ return &jacobi; }

/**
 *	@brief Set the vector of Jacobi constant values for this trajectory
 *	@param j a vector of Jacobi constants
 */
void adtk_cr3bp_traj::setJC(std::vector<double> j){ jacobi = j; }

/**
 *	@brief Retrieve data about the system this trajectory was propagated in
 *	@return the system data object
 */
adtk_cr3bp_sys_data adtk_cr3bp_traj::getSysData(){ return sysData; }


/**
 *	@brief Set the system data for this trajectory
 *	@param d a system data object
 */
void adtk_cr3bp_traj::setSysData(adtk_cr3bp_sys_data d){ sysData = d; }

adtk_sys_data::system_t adtk_cr3bp_traj::getType() const{
	return sysData.getType();
}

/**
 *	@brief Calls the basic trajectory setLength() method and implements extra catches
 *	specific to the CR3BP trajectory object
 */
void adtk_cr3bp_traj::setLength(){
	adtk_trajectory::setLength();
	
	if(jacobi.size() != times.size()){
		printErr("Warning: Jacobi vector has different length than time vector!\n");
	}
}
//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	@brief Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void adtk_cr3bp_traj::saveToMat(const char* filename){
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
		saveJacobi(matfp);
		sysData.saveToMat(matfp);
	}

	Mat_Close(matfp);
}//========================================

/**
 *	@brief Save the Jacobi vector to a file
 * 	@param matFile a pointer to the destination matlab file 
 */
void adtk_cr3bp_traj::saveJacobi(mat_t *matFile){
	size_t dims[2] = {static_cast<size_t>(numPoints), 1};
	matvar_t *matvar = Mat_VarCreate("Jacobi", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(jacobi[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "Jacobi", MAT_COMPRESSION_NONE);
}//=================================================





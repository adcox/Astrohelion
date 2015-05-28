/**
 *	@file adtk_bcr4bpr_traj.cpp
 *	
 */

#include "adtk_bcr4bpr_traj.hpp"
#include <iostream>

using namespace std;

//-----------------------------------------------------
// 		Constructor Functions
//-----------------------------------------------------

/**
 *	Construct a basic BCR4BPR trajectory object
 */
adtk_bcr4bpr_traj::adtk_bcr4bpr_traj() : adtk_trajectory() {
	dqdT.assign(6,0);
}

/**
 *	Construct a BCR4BPR trajectory object for the specified system
 *	@param data a system data object describing the BCR4BPR system
 */
adtk_bcr4bpr_traj::adtk_bcr4bpr_traj(adtk_bcr4bpr_sys_data data){
	dqdT.assign(6,0);
	sysData = data;
}

/**
 *	Construct a BCR4BPR trajectory object with room for a specified number of states
 *	@param n the number of states this trajectory will contain
 */
adtk_bcr4bpr_traj::adtk_bcr4bpr_traj(int n) : adtk_trajectory(n){
	dqdT.assign(n*6,0);
}

/**
 *	Copy the specified trajectory
 *	@param t a BCR4BPR trajectory object
 */
adtk_bcr4bpr_traj::adtk_bcr4bpr_traj(const adtk_bcr4bpr_traj &t) : adtk_trajectory(t){
	sysData = t.sysData;
	dqdT = t.dqdT;
}

//-----------------------------------------------------
// 		Operator Functions
//-----------------------------------------------------

/**
 *	Copy operator; copy a trajectory object into this one.
 *	@param t a trajectory object
 *	@return this trajectory object
 */
adtk_bcr4bpr_traj& adtk_bcr4bpr_traj::operator= (const adtk_bcr4bpr_traj& t){
	adtk_trajectory::operator= (t);
	sysData = t.sysData;
	dqdT = t.dqdT;
	return *this;
}//=====================================

/**
 *	Sum two BCR4BPR trajectories
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
adtk_bcr4bpr_traj operator +(const adtk_bcr4bpr_traj &lhs, const adtk_bcr4bpr_traj &rhs){
	if(lhs.sysData.getPrimary(0).compare(rhs.sysData.getPrimary(0)) != 0 ||
			lhs.sysData.getPrimary(1).compare(rhs.sysData.getPrimary(1)) != 0 || 
			lhs.sysData.getPrimary(2).compare(rhs.sysData.getPrimary(2)) != 0){
		fprintf(stderr, "Cannot sum two BCR4BPR trajectories from different systems!\n");
		throw;
	}

	// create a new trajectory object with space for both sets of data to be combined
	adtk_bcr4bpr_traj newTraj(lhs.numPoints + rhs.numPoints);
	
	// Copy the states and times from the LHS into the new guy
	copy(lhs.state.begin(), lhs.state.end(), newTraj.getState()->begin());
	copy(lhs.times.begin(), lhs.times.end(), newTraj.getTime()->begin());
	
	// Append the rhs state and timeto the end of the new guy's vectors; don't adjust the
	// time because the system is non-autonomous
	copy(rhs.state.begin(), rhs.state.end(), newTraj.getState()->begin() + lhs.numPoints);
	copy(rhs.times.begin(), rhs.times.end(), newTraj.getTime()->begin() + lhs.numPoints);

	// Copy only the lhs STMs, because there is no gaurantee the two summed trajectories
	// will be continuous and smooth in time and state; same for dqdT
	copy(lhs.allSTM.begin(), lhs.allSTM.end(), newTraj.getSTM()->begin());
	copy(lhs.dqdT.begin(), lhs.dqdT.end(), newTraj.get_dqdT()->begin());
	
	return newTraj;
}//========================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@return the angle between the P1/P2 line and the inertial x-axis, radians
 */
double adtk_bcr4bpr_traj::getTheta0(){ return sysData.getTheta0(); }

/**
 *	@return the angle between the P2/P3 line (projected into the inertial XY plane)
 *	and the inertial x-axis, radians
 */
double adtk_bcr4bpr_traj::getPhi0(){ return sysData.getPhi0(); }

/**
 *	@return the inclination of the P2/P3 orbital plane relative to the P1/P2 orbital
 *	plane, radians
 */
double adtk_bcr4bpr_traj::getGamma(){ return sysData.getGamma(); }

/**
 *	@return the system data object describing this system
 */
adtk_bcr4bpr_sys_data adtk_bcr4bpr_traj::getSysData(){ return sysData; }

/**
 *	Retrieve a pointer to the dqdT array for in-place editing.
 *	
 *	@return a pointer to the vector of dqdT values;
 */
vector<double>* adtk_bcr4bpr_traj::get_dqdT(){ return &dqdT; }

/**
 *	@param i the index of the dqdT vector to retrieve
 *	@return the i'th 6-element dqdT vector
 */
vector<double> adtk_bcr4bpr_traj::get_dqdT(int i){
	vector<double> temp(dqdT.begin()+i*6, dqdT.begin()+(i+1)*6);
	return temp;
}//===============================================

adtk_sys_data::system_t adtk_bcr4bpr_traj::getType() const{
	return sysData.getType();
}

/**
 *	Set the system data object
 *	@param data a data object describing the BCR4BP
 */
void adtk_bcr4bpr_traj::setSysData(adtk_bcr4bpr_sys_data data){ sysData = data; }

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	Calls the basic trajectory setLength() method and implements extra catches
 *	specific to the BCR4BPR trajectory object
 */
void adtk_bcr4bpr_traj::setLength(){
	adtk_trajectory::setLength();

	if(dqdT.size()/6 != times.size()){
		fprintf(stderr, "Warning: dqdT vector has different length than time vector!\n");
	}
}//================================

/**
 *	Save the trajectory to a file
 *	@param filename the name of the .mat file
 */
void adtk_bcr4bpr_traj::saveToMat(const char* filename){
	// TODO: Check for propper file extension, add if necessary

	/*	Create a new Matlab MAT file with the given name and optional
	 *	header string. If no header string is given, the default string 
	 *	used containing the software, version, and date in it. If a header
	 *	string is specified, at most the first 116 characters are written to
	 *	the file. Arguments are:
	 *	const char *matname 	- 	the name of the file
	 *	const char *hdr_str 	- 	the 116 byte header string
	 *	enum mat_ft 			- 	matlab file version: MAT_FT_MAT5 or MAT_FT_MAT4
	 */
	mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		fprintf(stderr, "Error creating MAT file\n");
	}else{
		saveState(matfp);
		saveTime(matfp);
		saveSTMs(matfp);
		save_dqdT(matfp);
		saveSysData(matfp);
	}

	Mat_Close(matfp);
}//========================================

/**
 *	Save the Jacobi vector to a file
 * 	@param matFile a pointer to the destination matlab file 
 */
void adtk_bcr4bpr_traj::save_dqdT(mat_t *matFile){
	size_t dims[2] = {static_cast<size_t>(numPoints), 6};
	matvar_t *matvar = Mat_VarCreate("dqdT", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(dqdT[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, matvar, "dqdT", MAT_COMPRESSION_NONE);
}//=================================================

/**
 *	Save system data, like the names of the primaries and the system mass ratio, to a .mat file
 *	@param matFile a pointer to the .mat file
 */
void adtk_bcr4bpr_traj::saveSysData(mat_t *matFile){
	size_t dims[2] = {1,1};

	// Initialize character array (larger than needed), copy in the name of the primary, then create a var.
	char p1_str[64];
	strcpy(p1_str, sysData.getPrimary(0).c_str());
	dims[1] = sysData.getPrimary(0).length();
	matvar_t *p1_var = Mat_VarCreate("P1", MAT_C_CHAR, MAT_T_UTF8, 2, dims, p1_str, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p1_var, "P1", MAT_COMPRESSION_NONE);

	char p2_str[64];
	strcpy(p2_str, sysData.getPrimary(1).c_str());
	dims[1] = sysData.getPrimary(1).length();
	matvar_t *p2_var = Mat_VarCreate("P2", MAT_C_CHAR, MAT_T_UTF8, 2, dims, &(p2_str[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p2_var, "P2", MAT_COMPRESSION_NONE);

	char p3_str[64];
	strcpy(p3_str, sysData.getPrimary(2).c_str());
	dims[1] = sysData.getPrimary(2).length();
	matvar_t *p3_var = Mat_VarCreate("P3", MAT_C_CHAR, MAT_T_UTF8, 2, dims, &(p3_str[0]), MAT_F_DONT_COPY_DATA);
	saveVar(matFile, p3_var, "P3", MAT_COMPRESSION_NONE);

	dims[1] = 1;
	double theta0_val = sysData.getTheta0();
	matvar_t *theta0_var = Mat_VarCreate("Theta0", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &theta0_val, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, theta0_var, "Theta0", MAT_COMPRESSION_NONE);

	double phi0_val = sysData.getPhi0();
	matvar_t *phi0_var = Mat_VarCreate("Phi0", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &phi0_val, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, phi0_var, "Phi0", MAT_COMPRESSION_NONE);

	double gamma_val = sysData.getGamma();
	matvar_t *gamma_var = Mat_VarCreate("Gamma", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &gamma_val, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, gamma_var, "Gamma", MAT_COMPRESSION_NONE);
}//===================================================
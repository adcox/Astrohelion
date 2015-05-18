/**
 *	@file adtk_cr3bp_traj.cpp
 *
 *	CR3BP Trajectory
 */

#include "adtk_cr3bp_traj.hpp"

using namespace std;

//-----------------------------------------------------
// 		Constructor Functions
//-----------------------------------------------------

/**
 *	Default constructor; calls constructor for super-class adtk_trajectory
 *	and additionall initializes the jacobi matrix
 */
adtk_cr3bp_traj::adtk_cr3bp_traj() : adtk_trajectory(){
	jacobi.assign(1,0);
}

/**
 *	Initialize all vectors to have size n; fill each vector with zeros.
 */
adtk_cr3bp_traj::adtk_cr3bp_traj(int n) : adtk_trajectory(n){
	jacobi.assign(n, 0);
}

//-----------------------------------------------------
// 		Operators
//-----------------------------------------------------

/**
 *	Copy operator; copy a trajectory object into this one.
 *	@param t a trajectory object
 *	@return this trajectory object
 */
adtk_cr3bp_traj& adtk_cr3bp_traj::operator= (const adtk_cr3bp_traj& t){
	adtk_trajectory::operator= (t);
	jacobi = t.jacobi;
	return *this;
}//=====================================

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	Get the value of Jacobi constant at a particular point on the trajectory
 *	@param n the index of the point
 *	@return the Jacobi constant (non-dimensional)
 */
double adtk_cr3bp_traj::getJC(int n){ return jacobi[n]; }

/**
 *	Retrieve a pointer to the Jacobi array for in-place editing.
 *	
 *	@return a pointer to the vector of Jacobi constants;
 */
vector<double>* adtk_cr3bp_traj::getJC(){ return &jacobi; }

/**
 *	Set the vector of Jacobi constant values for this trajectory
 *	@param j a vector of Jacobi constants
 */
void adtk_cr3bp_traj::setJC(std::vector<double> j){ jacobi = j; }

/**
 *	Retrieve data about the system this trajectory was propagated in
 *	@return 
 */
adtk_cr3bp_sys_data adtk_cr3bp_traj::getSysData(){ return sysData; }


/**
 *	Set the system data for this trajectory
 *	@param d
 */
void adtk_cr3bp_traj::setSysData(adtk_cr3bp_sys_data d){ sysData = d; }

//-----------------------------------------------------
// 		Utility Functions
//-----------------------------------------------------

/**
 *	Save the trajectory to a file
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
	 *	enum mat_ft 			- 	matlab file version: MAT_FT_MAT5 or MAT_FT_MAT4
	 */
	mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		fprintf(stderr, "Error creating MAT file\n");
	}else{
		saveState(matfp);
		saveTime(matfp);
		saveSTMs(matfp);
		saveJacobi(matfp);
	}

	Mat_Close(matfp);
}//========================================

/**
 *	Save the Jacobi vector to a file
 * 	@param matFile a pointer to the destination matlab file 
 */
void adtk_cr3bp_traj::saveJacobi(mat_t *matFile){
	size_t dims[2] = {static_cast<size_t>(numPoints), 1};
	matvar_t *matvar = Mat_VarCreate("Jacobi", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(jacobi[0]), MAT_F_DONT_COPY_DATA);
	if(NULL == matvar){
		fprintf(stderr, "Error creating variable for 'Jacobi'\n");
	}else{
		Mat_VarWrite(matFile, matvar, MAT_COMPRESSION_NONE);
		Mat_VarFree(matvar);
	}
}//=================================================
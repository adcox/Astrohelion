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

adtk_cr3bp_traj::adtk_cr3bp_traj(adtk_cr3bp_sys_data data){
	jacobi.assign(1,0);
	sysData = data;
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
		saveSysData(matfp);
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
	saveVar(matFile, matvar, "Jacobi", MAT_COMPRESSION_NONE);
}//=================================================

/**
 *	Save system data, like the names of the primaries and the system mass ratio, to a .mat file
 *	@param matFile a pointer to the .mat file
 */
void adtk_cr3bp_traj::saveSysData(mat_t *matFile){
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

	double mu_val = sysData.getMu();
	dims[1] = 1;	
	matvar_t *mu_var = Mat_VarCreate("Mu", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &mu_val, MAT_F_DONT_COPY_DATA);
	saveVar(matFile, mu_var, "Mu", MAT_COMPRESSION_NONE);
}//===================================================





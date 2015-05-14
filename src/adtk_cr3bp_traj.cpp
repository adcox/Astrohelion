/**
 *	CR3BP Trajectory
 */

#include "adtk_cr3bp_traj.hpp"

using namespace std;

//-----------------------------------------------------
// 		Constructor Functions
//-----------------------------------------------------

adtk_cr3bp_traj::adtk_cr3bp_traj() : adtk_trajectory(){
	jacobi.assign(1,0);
}


adtk_cr3bp_traj::adtk_cr3bp_traj(int n) : adtk_trajectory(n){
	jacobi.assign(n, 0);
}

//-----------------------------------------------------
// 		Operators
//-----------------------------------------------------

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
void adtk_cr3bp_traj::setJC(vector<double> j){ jacobi = j; }

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
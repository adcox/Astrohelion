/**
 *	@file adtk_bcr4bpr_traj.cpp
 *	
 */

#include "adtk_bcr4bpr_traj.hpp"

using namespace std;

//-----------------------------------------------------
// 		Constructor Functions
//-----------------------------------------------------

adtk_bcr4bpr_traj::adtk_bcr4bpr_traj() : adtk_trajectory() {
	theta0 = 0;
	phi0 = 0;
	gamma = 0;
}

adtk_bcr4bpr_traj::adtk_bcr4bpr_traj(int n) : adtk_trajectory(n){
	theta0 = 0;
	phi0 = 0;
	gamma = 0;
}

//-----------------------------------------------------
// 		Set and Get Functions
//-----------------------------------------------------

/**
 *	@return the angle between the P1/P2 line and the inertial x-axis, radians
 */
double adtk_bcr4bpr_traj::getTheta0(){ return theta0; }

/**
 *	@return the angle between the P2/P3 line (projected into the inertial XY plane)
 *	and the inertial x-axis, radians
 */
double adtk_bcr4bpr_traj::getPhi0(){ return phi0; }

/**
 *	@return the inclination of the P2/P3 orbital plane relative to the P1/P2 orbital
 *	plane, radians
 */
double adtk_bcr4bpr_traj::getGamma(){ return gamma; }

/**
 *	@return the system data object describing this system
 */
adtk_bcr4bpr_sys_data adtk_bcr4bpr_traj::getSysData(){ return sysData; }

/**
 *	Set <tt>theta0</tt> to <tt>t</tt>
 *	@param t the angle between the P1/P2 line and the inertial x-axis, radians
 */
void adtk_bcr4bpr_traj::setTheta0(double t){ theta0 = t; }

/**
 *	Set <tt>phi0</tt> to <tt>p</tt>
 *	@param p the angle between the P2/P3 line (projected into the inertial XY plane)
 *	and the inertial x-axis, radians
 */
void adtk_bcr4bpr_traj::setPhi0(double p){ phi0 = p; }

/**
 *	Set <tt>gamma</tt> to <tt>g</tt>
 *	@param t the inclination of the P2/P3 orbital plane relative to the P1/P2 orbital
 *	plane, radians
 */
void adtk_bcr4bpr_traj::setGamma(double g){ gamma = g; }

/**
 *	Set the system data object
 *	@param data a data object describing the BCR4BP
 */
void adtk_bcr4bpr_traj::setSysData(adtk_bcr4bpr_sys_data data){ sysData = data; }
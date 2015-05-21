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

adtk_bcr4bpr_traj::adtk_bcr4bpr_traj() : adtk_trajectory() {
	dqdT.assign(1,0);
}

adtk_bcr4bpr_traj::adtk_bcr4bpr_traj(adtk_bcr4bpr_sys_data data){
	dqdT.assign(1,0);
	sysData = data;
}

adtk_bcr4bpr_traj::adtk_bcr4bpr_traj(int n) : adtk_trajectory(n){
	dqdT.assign(n,0);
}

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
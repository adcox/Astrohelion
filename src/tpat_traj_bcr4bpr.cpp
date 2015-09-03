/**
 *  @file tpat_traj_bcr4bpr.cpp
 *	@brief Derivative of tpat_traj, specific to BCR4BPR
 *
 *	@author Andrew Cox
 *	@version September 2, 2015
 *	@copyright GNU GPL v3.0
 */
 
/*
 *  Trajectory Propagation and Analysis Toolkit 
 *  Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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

#include "tpat_traj_bcr4bpr.hpp"

#include "tpat_exceptions.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_utilities.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

tpat_traj_bcr4bpr::tpat_traj_bcr4bpr(tpat_sys_data_bcr4bpr *sys) : tpat_traj(sys){
	initExtraParam();
}//====================================================

tpat_traj_bcr4bpr::tpat_traj_bcr4bpr(const tpat_traj_bcr4bpr &t) : tpat_traj(t){
	initExtraParam();
}//====================================================

tpat_traj_bcr4bpr::tpat_traj_bcr4bpr(const tpat_arc_data &a) : tpat_traj(a){
	initExtraParam();
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@return the angle between the P1/P2 line and the inertial x-axis, radians
 */
double tpat_traj_bcr4bpr::getTheta0(){
	tpat_sys_data_bcr4bpr *bcSys = static_cast<tpat_sys_data_bcr4bpr *>(sysData);
	return bcSys->getTheta0();
}//====================================================

/**
 *	@return the angle between the P2/P3 line (projected into the inertial XY plane)
 *	and the inertial x-axis, radians
 */
double tpat_traj_bcr4bpr::getPhi0(){
	tpat_sys_data_bcr4bpr *bcSys = static_cast<tpat_sys_data_bcr4bpr *>(sysData);
	return bcSys->getPhi0();
}//====================================================

/**
 *	@return the inclination of the P2/P3 orbital plane relative to the P1/P2 orbital
 *	plane, radians
 */
double tpat_traj_bcr4bpr::getGamma(){
	tpat_sys_data_bcr4bpr *bcSys = static_cast<tpat_sys_data_bcr4bpr *>(sysData);
	return bcSys->getGamma();
}//====================================================

/**
 *	@param ix the index of the dqdT vector to retrieve
 *	@return the i'th 6-element dqdT vector. If ix is negative, the count
 *	will proceed from the end of the vector, i.e. -1 will return the final time, 
 *	-2 will give the second to last value, etc.
 */
std::vector<double> tpat_traj_bcr4bpr::get_dqdT(int ix){
	if(ix < 0)
		ix += steps.size();

	return getExtraParam(ix, 2);
}//====================================================

/**
 *	@brief Set the value of the dqdT vector for the specified step
 *	@param ix the index of the step; if < 0, it will count backwards from the end
 *	@param dqdT a pointer to the dqdT vector; this MUST have at least 6 elements,
 *	or the function will read unallocated memory.
 */
void tpat_traj_bcr4bpr::set_dqdT(int ix, double *dqdT){
	if(ix < 0)
		ix += steps.size();

	for(int i = 0; i < 6; i++)
		steps[ix].setExtraParam(1+i, dqdT[i]);
}//====================================================

/**
 *	@brief Set the value of the dqdT vector for the specified step
 *	@param ix the index of the step; if < 0, it will count backwards from the end
 *	@param dqdT a vector (6 elements) representing the dqdT vector
 */
void tpat_traj_bcr4bpr::set_dqdT(int ix, std::vector<double> dqdT){
	if(dqdT.size() != 6)
		throw tpat_exception("tpat_traj_bcr4bpr::set_dqdT: Cannot accept a dqdT with anything other than 6 elements");

	set_dqdT(ix, &(dqdT[0]));
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

void tpat_traj_bcr4bpr::initExtraParam(){
	// This function in tpat_traj was already called, so 
	// numExtraParam has been set to 1 and a row size has
	// been appended for the time variable

	// Add another variable for dqdT
	numExtraParam = 2;
	extraParamRowSize.push_back(6);
}//====================================================
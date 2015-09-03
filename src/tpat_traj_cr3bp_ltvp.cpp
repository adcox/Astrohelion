/**
 *  @file tpat_traj_cr3bp_ltvp.cpp
 *	@brief 
 *
 *	@author Andrew Cox
 *	@version 
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

#include "tpat_traj_cr3bp_ltvp.hpp"

#include "tpat_sys_data_cr3bp_ltvp.hpp"
//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

tpat_traj_cr3bp_ltvp::tpat_traj_cr3bp_ltvp(tpat_sys_data_cr3bp_ltvp* sys) : tpat_traj(sys){
	initExtraParam();
}//====================================================

tpat_traj_cr3bp_ltvp::tpat_traj_cr3bp_ltvp(const tpat_traj_cr3bp_ltvp &t) : tpat_traj(t){
	initExtraParam();
}//====================================================

tpat_traj_cr3bp_ltvp::tpat_traj_cr3bp_ltvp(const tpat_arc_data &a) : tpat_traj(a){
	initExtraParam();
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

double tpat_traj_cr3bp_ltvp::getJacobi(int ix) const{
	if(ix < 0)
		ix += steps.size();
	tpat_arc_step step = steps[ix];
	return step.getExtraParam(1);
}//====================================================

double tpat_traj_cr3bp_ltvp::getMass(int ix) const{
	if(ix < 0)
		ix += steps.size();
	tpat_arc_step step = steps[ix];
	return step.getExtraParam(2);
}//====================================================

void tpat_traj_cr3bp_ltvp::setJacobi(int ix, double val){
	if(ix < 0)
		ix += steps.size();

	steps[ix].setExtraParam(1, val);
}//====================================================

void tpat_traj_cr3bp_ltvp::setMass(int ix, double val){
	if(ix < 0)
		ix += steps.size();

	steps[ix].setExtraParam(2, val);
}//====================================================

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

void tpat_traj_cr3bp_ltvp::initExtraParam(){
	// This function in tpat_traj was already called, so 
	// numExtraParam has been set to 1 and a row size has
	// been appended for the time variable

	// Add another variable for Jacobi Constant, and one for mass
	numExtraParam = 3;
	extraParamRowSize.push_back(1);	// add var for Jacobi
	extraParamRowSize.push_back(1); // add var for Mass
}//====================================================
/**
 *  @file tpat_traj_step.cpp
 *	@brief Data object that stores information about a node
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

#include "tpat_traj_step.hpp"

//-----------------------------------------------------
//      *structors
//-----------------------------------------------------

/**
 *	@brief Construct a step from the state and time
 *	@param state a 6-element vector containing the non-dimensional
 *	state at one step along an integrated arc
 *	@param t the non-dimensional time associated with the state
 */
tpat_traj_step::tpat_traj_step(double *state, double t){
	initArrays();
	std::copy(state, state+6, posVelState);
	extraParam[0] = t;
}//====================================================

/**
 *	@brief Construct a step with all parameters specified
 *	@param state a 6-element vector containing the non-dimensional
 *	state at one step along an integrated arc
 *	@param t the non-dimensional time associated with the state
 *	@param accel a 3-element vector containing the non-dimensional 
 *	acceleration at this step on the arc
 *	@param stm a 36-element vector containing the STM elements for
 *	this step on the arc
 */
tpat_traj_step::tpat_traj_step(double *state, double t,
	double *accel, double *stm){
	
	initArrays();
	std::copy(state, state+6, posVelState);
	extraParam[0] = t;
	std::copy(accel, accel+3, this->accel);
	std::copy(stm, stm+36, this->stm);
}//====================================================

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

double tpat_traj_step::getTime() const { return extraParam[0]; }

void tpat_traj_step::setTime(double t){ extraParam[0] = t; }

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

void tpat_traj_step::initArrays(){
	tpat_arc_step::initArrays();
	extraParam.assign(1,0);
}
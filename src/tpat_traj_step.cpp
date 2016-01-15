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
tpat_traj_step::tpat_traj_step(const double *state, double t){
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
tpat_traj_step::tpat_traj_step(const double *state, double t, const double *accel, const double *stm){
	initArrays();
	std::copy(state, state+6, posVelState);
	extraParam[0] = t;
	std::copy(accel, accel+3, this->accel);
	std::copy(stm, stm+36, this->stm);
}//====================================================

/**
 *	@brief Copy constructor
 *	@param s a tpat_traj_step reference
 */
tpat_traj_step::tpat_traj_step(const tpat_traj_step &s) : tpat_arc_step(s) {}

/**
 *	@brief Create a trajectory step object from a generic arc
 *	step object.
 *
 *	This is permissible because traj_step only defines
 *	new access methods, not new data objects
 *	@param s a tpat_arc_step reference
 */
tpat_traj_step::tpat_traj_step(const tpat_arc_step &s) : tpat_arc_step(s) {}

//-----------------------------------------------------
//      Operators
//-----------------------------------------------------

//-----------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------

/**
 *	@brief Get the time at this step
 *	@return the non-dimensional integration time at this step
 */
double tpat_traj_step::getTime() const { return extraParam[0]; }

/**
 *	@brief Set the time for this step
 *	@param t non-dimensional integration time at this step
 */
void tpat_traj_step::setTime(double t){ extraParam[0] = t; }

//-----------------------------------------------------
//      Utility Functions
//-----------------------------------------------------

/**
 *	@brief initialize data arrays for this object
 */
void tpat_traj_step::initArrays(){
	tpat_arc_step::initArrays();
	extraParam.assign(1,0);
}
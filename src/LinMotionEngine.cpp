/**
 *	\file LinMotionEngine.cpp
 *	\brief Uses linear EOMS near libration points to generate trajectories
 *	
 *	\author Andrew Cox
 *	\version September 28, 2017
 *	\copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of Astrohelion
 *
 *  Astrohelion is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Astrohelion is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Astrohelion.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "LinMotionEngine.hpp"

namespace astrohelion{
//-----------------------------------------------------------------------------
//      *structors
//-----------------------------------------------------------------------------

	// Included in header file

//-----------------------------------------------------------------------------
//      Set and Get Functions
//-----------------------------------------------------------------------------

/**
 *	\brief Retrieve the number of revolutions to simulate for
 *
 *	Note that this rev count applies to the in-plane oscillations only; the out-of-plane
 *	motion will most likely have a different period and therefore perform a different
 *	number of revs during the same time period
 *
 *	\return the number of rotatins to simulate for
 */
int LinMotionEngine::getNumRevs() const { return revs; }

/**
 *	\brief Retrieve the step size (in non-dimensional units) for the time vector
 *	\return he step size (in non-dimensional units) for the time vector
 */
double LinMotionEngine::getTimeStep() const { return t_step; }

/**
 *	\brief Retrieve the acceptable numerical tolerance for computations
 *
 *	This tolerance is used to target Lagrange point locations and to determine
 *	where or not numbers are "equal" to zero (or any other value)
 *	\return the acceptable numerical tolerance for computations
 */
double LinMotionEngine::getTol() const { return tol; }

/**
 *	\brief Set the number of revolutions to simulate for.
 *
 *	The length of one revolutions is determined by the in-plane period, not the
 *	out-of-plane period
 *	\param numRevs number of revolutions
 */
void LinMotionEngine::setNumRevs(int numRevs) { revs = numRevs; }

/**
 *	\brief Set the step size for the time vector
 *	\param dt non-dimensional time step
 */
void LinMotionEngine::setTimeStep(double dt) { t_step = dt; }

/**
 *	\brief set the tolerance to use
 *	\param t tolerance, non-dimensional units
 */
void LinMotionEngine::setTol(double t){ tol = t; }

/**
 *	\brief get a human-readable string for a motion type
 *	\param type the motion type
 *	\return a human-redable string
 */
const char* LinMotionEngine::getTypeStr(unsigned int type) const{
	switch(type){
		case LinMotion_tp::NONE: return "NONE";
		case LinMotion_tp::HYP: return "HYPERBOLIC";
		case LinMotion_tp::OSC: return "OSCILLATORY";
		case LinMotion_tp::STAB_OSC: return "STABLE OSCILLATORY";
		case LinMotion_tp::UNSTAB_OSC: return "UNSTABLE OSCILLATORY";
		default: return "Unrecognized type...";
	}
}//====================================================

//-----------------------------------------------------------------------------
//      Utility Functions
//-----------------------------------------------------------------------------

void LinMotionEngine::cleanEngine(){
	bIsClean = true;
	// Nothing to do here
}//====================================================

void LinMotionEngine::reset(){
	if(!bIsClean)
		cleanEngine();

	t_step = 0.001;
	revs = 1;
	tol = 1e-14;
}//====================================================

}// END of Astrohelion namespace
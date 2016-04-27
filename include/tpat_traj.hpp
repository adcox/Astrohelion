/*
 *	Trajectory Propagation and Analysis Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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

#ifndef H_TRAJECTORY
#define H_TRAJECTORY

#include "tpat_arc_data.hpp"
#include "matio.h"

// Forward Declarations
class tpat_traj_step;
class tpat_nodeset;

/**
 *	@brief Contains information about a series of continuous states along a trajectory.
 *	
 *	The trajecotory object holds information about a trajectory in one package
 *	to make passing the data between engines and analysis tools simple. This class acts
 *	as a base class for model-specific trajectory objects. This general object stores the
 *	following information about a trajectory:
 *	* A vector of trajectory steps, each of which contains:
 *		* State  (x, y, z, vx, vy, vz [non-dimensional])
 *		* Accleration (ax, ay, az [non-dimensional])
 *		* State Transition Matrix (6x6 [non-dimensional])
 *		* Time since first step [non-dimensional]
 * 	* A system data object that describes the system the trajectory was integrated in
 * 	* A tolerance value that describes the minimum numerical accuracy used to create the trajectory
 * 	
 *	@author Andrew Cox
 *	@version August 29, 2015
 *	@copyright GNU GPL v3.0
 *	
 *	@see tpat_arc_data
 */
class tpat_traj : public tpat_arc_data{

public:
	// *structors
	tpat_traj(const tpat_sys_data*);
	tpat_traj(const tpat_traj&);
	tpat_traj(const tpat_arc_data&);

	// Operators

	// Set and Get Functions
	double getTime(int) const;
	double getTotalTOF() const;
	tpat_traj_step getStep(int) const;

	void appendStep(tpat_traj_step);
	void setTime(int, double);

	void shiftAllTimes(double);
	
	// Utility Functions
	tpat_nodeset discretize(int) const;
	void initExtraParam();
	void print() const;
	virtual void readFromMat(const char*);
	virtual void saveToMat(const char*) const;

protected:
	void saveTime(mat_t*) const;
};

#endif
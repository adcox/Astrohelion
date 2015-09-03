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

/**
 *	@brief Contains information about a series of continuous states along a trajectory.
 *
 *	This object holds information about a trajectory in one package to
 *	make passing the data between engines and analysis tools symple.
 *	This class acts as a template for derivative classes that apply
 *	to specific systems. For example, the CR3BP has a specific
 *	trajectory class which includes additional information specific
 *	to the CR3BP, like Jacobi constant.
 *
 *	@author Andrew Cox
 *	@version August 29, 2015
 *	@copyright GNU GPL v3.0
 */
class tpat_traj : public tpat_arc_data{

public:
	// *structors
	tpat_traj(tpat_sys_data*);
	tpat_traj(const tpat_traj&);
	tpat_traj(const tpat_arc_data&);
	
	// Set and Get Functions
	double getTime(int) const;
	tpat_traj_step getStep(int) const;

	void appendStep(tpat_traj_step);
	
	// Utility Functions
	void saveToMat(const char*);
	void print() const;
	void initExtraParam();

protected:

	void saveTime(mat_t*);
};

#endif
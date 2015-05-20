/**
 *	This object holds information about a trajectory in one package to
 *	make passing the data between engines and analysis tools symple.
 *	This class acts as a template for derivative classes that apply
 *	to specific systems. For example, the CR3BP has a specific
 *	trajectory class which includes additional information specific
 *	to the CR3BP, like Jacobi constant.
 *
 *	Author: Andrew Cox
 *
 *	Version: May 15, 2015
 *
 *	TODO:
 *		- Overload +, += operators for concatenating trajectories
 */

/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (ADTK).
 *
 *  ADTK is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ADTK is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ATDK.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __H_TRAJECTORY_
#define __H_TRAJECTORY_

#include "matio.h"
 
#include <vector>

// Forward Declarations
class adtk_matrix;

class adtk_trajectory{
	static const int STATE_WIDTH = 9;

	public:
		// *structors
		adtk_trajectory();
		adtk_trajectory(int);
		virtual ~adtk_trajectory();

		// Operators
		adtk_trajectory& operator= (const adtk_trajectory&);

		// Set and Get functions
		int getLength();
		std::vector<double> getState(int);
		std::vector<double>* getState();
		double getTime(int);
		std::vector<double>* getTime();
		adtk_matrix getSTM(int);
		std::vector<adtk_matrix>* getSTM();

		void setLength();
		void setState(std::vector<double>);
		void setTime(std::vector<double>);
		void setSTMs(std::vector<adtk_matrix>);

		// Utility functions
		void saveToMat(const char*);

	protected:
		/** Number of points along integrated path */
		int numPoints;

		/** Holds state info: [pos, vel, accel] in 1D form, so every 9 elements
		 * 	constitutes a new "row" 
		 */
		std::vector<double> state;

		/** Holds time info */
		std::vector<double> times;

		/** An array containing the STM at each step of the integration */
		std::vector<adtk_matrix> allSTM;

		void saveState(mat_t*);
		void saveTime(mat_t*);
		void saveSTMs(mat_t*);
		void saveVar(mat_t*, matvar_t*, const char*, matio_compression);
};

#endif

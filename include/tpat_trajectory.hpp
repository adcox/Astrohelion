/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (TPAT).
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
#ifndef __H_TRAJECTORY_
#define __H_TRAJECTORY_

#include "tpat_matrix.hpp"
#include "tpat_sys_data.hpp"

#include "matio.h"
#include <vector>


/**
 *	@brief Contains a vector of states and times that fully describe a trajectory
 *
 *	This object holds information about a trajectory in one package to
 *	make passing the data between engines and analysis tools symple.
 *	This class acts as a template for derivative classes that apply
 *	to specific systems. For example, the CR3BP has a specific
 *	trajectory class which includes additional information specific
 *	to the CR3BP, like Jacobi constant.
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 *	@copyright GNU GPL v3.0
 *
 *	TODO:
 *		- Overload +, += operators for concatenating trajectories
 */
class tpat_trajectory{
	public:
		/** The width of one row of the state vector (state vector stored in row-major order) */
		static const int STATE_WIDTH = 9;
		
		// *structors
		tpat_trajectory();
		tpat_trajectory(int);
		tpat_trajectory(const tpat_trajectory&);
		virtual ~tpat_trajectory(){}

		// Operators
		tpat_trajectory& operator= (const tpat_trajectory&);

		// Set and Get functions
		int getLength() const;
		virtual tpat_sys_data::system_t getType() const;

		std::vector<double> getCoord(int) const;
		std::vector<double> getState(int) const;
		std::vector<double> getState_6(int) const;
		std::vector<double>* getState();
		tpat_matrix getSTM(int) const;
		std::vector<tpat_matrix>* getSTM();
		double getTime(int) const;
		std::vector<double>* getTime();

		void setState(std::vector<double>);
		void setSTMs(std::vector<tpat_matrix>);
		void setTime(std::vector<double>);

		// Utility functions
		void saveToMat(const char*);
		virtual void setLength();
	protected:
		int numPoints = 0;	//!< Number of points along integrated path

		/** Holds state info: [pos, vel, accel] in 1D form, so every 9 elements
		 * 	constitutes a new "row" 
		 */
		std::vector<double> state;
		std::vector<double> times;			//!< Holds time info
		std::vector<tpat_matrix> allSTM;	//!< An array containing the STM at each step of the integration

		void saveState(mat_t*);
		void saveTime(mat_t*);
		void saveSTMs(mat_t*);
};

#endif

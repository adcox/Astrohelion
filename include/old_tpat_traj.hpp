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
 *		- It probably makes more sense to make extraParams a vector of 
 *		  vectors to avoid confusion when reading and writing values that
 *		  are not 1-dimensional
 */
class tpat_traj{
	public:
		static const int STATE_SIZE = 6;	//!< Number of elements stored in a state vector
		static const int ACCEL_SIZE = 3;	//!< Number of elements stored in an acceleration vector

		// *structors
		tpat_traj();
		tpat_traj(int);
		tpat_traj(const tpat_traj&);
		virtual ~tpat_traj();

		// Operators
		tpat_traj& operator= (const tpat_traj&);

		// Set and Get functions
		std::vector<double> getAccel(int) const;
		std::vector<double>* getAccel();
		std::vector<double> getCoord(int) const;
		int getLength() const;
		std::vector<double> getExtraParam(int) const;
		std::vector<double>* getExtraParamPtr(int);
		std::vector<double> getState(int) const;
		std::vector<double>* getState();
		tpat_matrix getSTM(int) const;
		std::vector<tpat_matrix>* getSTM();
		virtual tpat_sys_data* getSysDataPtr();
		double getTime(int) const;
		std::vector<double>* getTime();
		double getTol() const;
		virtual tpat_sys_data::system_t getType() const;

		void setAccel(std::vector<double>);
		void setExtraParam(int, std::vector<double>);
		void setState(std::vector<double>);
		void setSTMs(std::vector<tpat_matrix>);
		void setTime(std::vector<double>);
		void setTol(double);
		
		// Utility functions
		friend void basicConcat(const tpat_traj*, const tpat_traj*, tpat_traj*);
		void saveToMat(const char*);
		void setLength();
	protected:
		int numPoints = 0;		//!< Number of points along integrated path
		double tol = 0;			//!< tolerance used to compute this trajectory
		int numExtraParam = 0;	//!< Number of variables stored in the extraParam vector

		std::vector<double> state;			//!< Holds [pos, vel] (6d) for every step
		std::vector<double> accel;			//!< Holds accelerations in three principle directions at every step
		std::vector<double> times;			//!< Holds time info
		std::vector< std::vector<double> > extraParam;	//!< Holds data for any extra parameters
		std::vector<int> extraParamRowSize;	//!< Number of elements in one row of each of the extra parameter vectors
		std::vector<tpat_matrix> allSTM;	//!< An array containing the STM at each step of the integration

		void copyMe(const tpat_traj&);
		void saveAccel(mat_t*);
		void saveExtraParam(mat_t*, int, int, const char*);
		void saveState(mat_t*);
		void saveTime(mat_t*);
		void saveSTMs(mat_t*);
};

#endif

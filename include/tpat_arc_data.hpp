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

#ifndef H_ARC_DATA
#define H_ARC_DATA

#include "tpat_arc_step.hpp"
#include "tpat_eigen_defs.hpp"
#include "tpat_sys_data.hpp"

#include "matio.h"
#include <vector>

/**
 *	@brief Abstract class that provides the framework for trajectories and nodesets
 *	
 *	The arc_data object specifies default and mandatory behaviors for all derivative
 *	classes (i.e. tpat_traj and tpat_nodeset). All variables and data for an arc or 
 *	one of its derivative classes are declared and stored here; in other words, no 
 *	derivative classes declare class-specific data objects. This architecture has been
 *	chosen to facilitate easy casting between model-specific derivative classes with
 *	the added bonus of being able to cast easily between, say, a trajectory and a
 *	nodeset.
 *	
 *	This class contains all information about any trajectory or nodeset in a few objects:
 *	* steps - a vector of tpat_arc_step objects, each of which contains information about
 *		the state, acceleration, STM, any any other parameter values at one single step
 *	* sysData - a pointer to a system data object that describes the system this arc
 *		has been generated in
 *	* numExtraParam and extraParamRowSize describe the number of extra parameters and 
 *		the number of elements in each parameter; this various for different dynamical
 *		model-specific derivative classes
 *	* tol - The maximum numerical tolerance with which the data in this object has been 
 *		computed
 *	
 *	The following behavior is mandatory for all derivative classes:
 *	* Ability to add two arcs together via operator +()
 *	* Ability to save to a matlab file via saveToMat()
 *	* Ability to display a textual representation of the object via print()
 *	
 *	Additionally, the following behavior is defined for all derivative classes, though
 *	they may override the default:
 *	* Assignment operator
 *	* Access to position and velocity values at each step via getState()
 *	* Access to acceleration values at each step via getAccel()
 *	* Access to any extra parameters that evolve each step via getExtraParam()
 *	* Access to the STM at each step via getSTM()
 *	* Access to individual step objects via getStep()
 *	* Access to the system data object pointer that describes the system this arc was integrated in
 *	
 */
class tpat_arc_data{

public:
	// *structors
	tpat_arc_data(const tpat_sys_data*);
	tpat_arc_data(const tpat_arc_data&);
	virtual ~tpat_arc_data();

	// Operators
	tpat_arc_data& operator =(const tpat_arc_data&);
	virtual tpat_arc_data& operator +(const tpat_arc_data&);

	// Set and Get functions
	std::vector<double> getAccel(int) const;
	std::vector<double> getCoord(int) const;
	int getLength() const;
	std::vector<double> getExtraParam(int, int) const;
	std::vector<double> getState(int) const;
	tpat_arc_step getStep(int) const;
	MatrixXRd getSTM(int) const;
	const tpat_sys_data* getSysData() const;
	double getTol() const;

	void appendStep(tpat_arc_step);
	void setAccel(int, std::vector<double>);
	void setState(int, std::vector<double>);
	void setSTM(int, MatrixXRd);

	void setTol(double);

	// Utility Functions

	/**
	 *	@brief Saves the object to a Matlab binary file
	 */
	virtual void saveToMat(const char*) const = 0;

	/**
	 *	@brief Displays a useful messages about the object
	 */
	virtual void print() const = 0;

	void updateCons();
	
protected:
	/** Contains all integration steps */
	std::vector<tpat_arc_step> steps;

	/** A pointer to the system data object that describes the system this arc exists in */
	const tpat_sys_data *sysData;

	/** 
	 *	Number of variables stored in the extraParam vector. This
	 *	parameter should be set by the constructor of all derived
	 *	classes
	 */
	int numExtraParam = 0;

	/** 
	 *	Number of elements in each extra parameter. This parameter
	 *	should be set by the constructor of all derived classes
	 */
	std::vector<int> extraParamRowSize;

	/** Tolerance used to compute this data */
	double tol = 0;

	void copyMe(const tpat_arc_data&);
	void saveAccel(mat_t*) const;
	void saveExtraParam(mat_t*, int, const char*) const;
	void saveState(mat_t*) const;
	void saveState(mat_t*, const char*) const;
	void saveSTMs(mat_t*) const;
};


#endif
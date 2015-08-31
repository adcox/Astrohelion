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
#include "tpat_sys_data.hpp"

#include "matio.h"
#include <vector>

/**
 *	@brief Abstract class that provides the framework for trajectories and nodesets
 */
class tpat_arc_data{

public:
	// *structors
	tpat_arc_data(tpat_sys_data*);
	tpat_arc_data(const tpat_arc_data&);
	virtual ~tpat_arc_data();

	// Operators
	tpat_arc_data& operator =(const tpat_arc_data&);
	friend tpat_arc_data& operator +(const tpat_arc_data&, const tpat_arc_data&);

	// Set and Get functions
	std::vector<double> getAccel(int) const;
	std::vector<double> getCoord(int) const;
	int getLength() const;
	std::vector<double> getExtraParam(int, int) const;
	std::vector<double> getState(int) const;
	tpat_arc_step getStep(int) const;
	tpat_matrix getSTM(int) const;
	tpat_sys_data* getSysData();
	double getTol() const;

	void setTol(double);

	// Utility Functions

	/**
	 *	@brief Saves the object to a Matlab binary file
	 */
	virtual void saveToMat(const char*) = 0;

	/**
	 *	@brief Displays a useful messages about the object
	 */
	virtual void print() const = 0;

protected:
	/** Contains all integration steps */
	std::vector<tpat_arc_step> steps;

	/** A pointer to the system data object that describes the system this arc exists in */
	tpat_sys_data *sysData;

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
	void saveAccel(mat_t*);
	void saveExtraParam(mat_t*, int, const char*);
	void saveState(mat_t*);
	void saveSTMs(mat_t*);
};


#endif
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

#ifndef H_ARC_STEP
#define H_ARC_STEP

#include "tpat_constraint.hpp"

// Forward Declarations
class tpat_matrix;

/**
 *	@brief Base class that represents a single integration step or node
 *	
 *	This class defines the default behavior for either of the more specific
 *	step types, tpat_traj_step and tpat_node. At their most basic level, both
 *	objects store the same information, hence this class defines all data objects
 *	and derivative classes provide only access functions to the data for methods
 *	more specific to their applications.
 *	
 *	An arc step stores the following data:
 *	* State (x, y, z, vx, vy, vz [non-dimensional]), accessed via getPosVelState()
 *	* Acceleration (ax, ay, az [non-dimensional]), accessed via getAccel()
 *	* Constraints (only applies to tpat_node), accessed via getConstraints()
 *	* State Transition Matrix (6x6), accessed via getSTM() or getSTMElements()
 *	* Extra Parameters, which may be useful for model-specific derivative classes.
 *	  	Parameters like Jacobi Constant, epoch, or mass are stored here
 *	* Flags - a vector of booleans that may be used as status variables. The 
 *		tpat_node class uses this data object to store information about velocity
 *		continuity
 */
class tpat_arc_step{

public:
	// *structors
	tpat_arc_step();
	tpat_arc_step(const tpat_arc_step&);
	virtual ~tpat_arc_step();

	// Operators
	tpat_arc_step& operator =(const tpat_arc_step&);
	friend bool operator ==(const tpat_arc_step&, const tpat_arc_step&);
	friend bool operator !=(const tpat_arc_step&, const tpat_arc_step&);

	// Set and Get functions
	std::vector<double> getAccel() const;
	std::vector<tpat_constraint> getConstraints() const;
	double getExtraParam(int) const;
	std::vector<double> getExtraParams() const;
	std::vector<double> getPosVelState() const;
	tpat_matrix getSTM() const;
	std::vector<double> getSTMElements() const;

	void addConstraint(tpat_constraint);
	void clearConstraints();
	void removeConstraint(int);
	void setConstraints(std::vector<tpat_constraint>);
	
	void setAccel(double*);
	void setAccel(std::vector<double>);
	void setExtraParam(int, double);
	void setExtraParams(std::vector<double>);
	void setPosVelState(double*);
	void setPosVelState(std::vector<double>);
	void setSTM(tpat_matrix);
	void setSTM(double*);
	void setSTM(std::vector<double>);

protected:
	virtual void copyMe(const tpat_arc_step&);
	virtual void initArrays();

	double posVelState[6]; 	//!< Stores 6 position and velocity states
	double accel[3];		//!< Stores 3 acceleration states
	double stm[36];			//!< Stores 36 STM elements (row-major order)

	/** Stores extra parameters like integration time, time-of-flight, epoch, etc. */
	std::vector<double> extraParam;

	/** Stores flags, which may be interpreted by derived classes */
	std::vector<bool> flags;

	/** Stores constraints on this arc step (especially usefull in nodesets) */
	std::vector<tpat_constraint> constraints;
};


#endif
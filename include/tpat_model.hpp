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
#ifndef H_MODEL_DEF
#define H_MODEL_DEF

// Forward Declarations
class tpat_event;
class tpat_traj;
class tpat_sys_data;

/**
 *	@brief A base class that defines the behavior of a dynamical model and provides
 *	functions used in simulation and correction algorithms
 */
class tpat_model{

public:
	enum dynamicModel_t{
		MODEL_NULL,			//!< Undefined type
		MODEL_CR3BP,		//!< Circular Restricted, 3-Body Problem
		MODEL_CR3BP_LTVP,	//!< Circular Restricted, 3-Body Problem with Velocity-Pointing Low Thrust (constant)
		MODEL_BCR4BPR		//!< Bi-Circular Restricted, 4-Body Problem in a rotating reference frame
	};

	typedef int (*eom_fcn)(double, const double[], double[], void*);

	// *structors
	tpat_model(dynamicModel_t);
	tpat_model(const tpat_model&);
	virtual ~tpat_model();

	// Operators
	tpat_model& operator =(const tpat_model&);

	// All derived classes MUST implement these methods
	/**
	 *	@brief Takes an input state and time and saves the data to the trajectory
	 *	@param y an array containing the core state and any extra states integrated
	 *	by the EOM function, including STM elements.
	 *	@param traj a pointer to the trajectory we should store the data in
	 */
	virtual void saveIntegratedData(double *y, double t, tpat_traj* traj) = 0;
	
	/**
	 *  @brief Use a correction algorithm to check and see if an event has occurred
	 *
	 *  The simulation engine calls this function if and when it determines that an event 
	 *  has been crossed. To accurately locate the event, we employ differential corrections
	 *  and find the exact event occurence in space and time.
	 *
	 *  @param event the event we're looking for
	 *  @param traj a pointer to the trajectory the event should occur on
	 *  @param model the dynamical model we're working in
	 *  @param ic the core state vector for this system
	 *  @param t0 non-dimensional time at the beginning of the search arc
	 *  @param tof the time-of-flight for the arc to search over
	 *  @param verbose whether or not we should be verbose with output messages
	 *
	 *  @param return wether or not the event has been located. If it has, a new point
	 *  has been appended to the trajectory's data vectors.
	 */
	virtual bool locateEvent(tpat_event event, tpat_traj *traj, tpat_model* model,
    	double *ic, double t0, double tof, bool verbose) = 0;

	/**
	 *	@brief Retrieve a pointer to the EOM function that computes derivatives
	 *	for only the core states (i.e. simple)
	 *
	 *	The EOM function must be a non-member function; we store them in the 
	 *	tpat_calculations.cpp file
	 */
	virtual eom_fcn getSimpleEOM_fcn() = 0;

	/**
	 *	@brief Retrieve a pointer to the EOM function that computes derivatives
	 *	for all states (i.e. full)
	 *
	 *	The EOM function must be a non-member function; we store them in the 
	 *	tpat_calculations.cpp file
	 */
	virtual eom_fcn getFullEOM_fcn() = 0;

	// Set and Get Functions
	int getCoreStateSize() const;
	int getSTMStateSize() const;
	int getExtraStateSize() const;

protected:
	dynamicModel_t modelType = MODEL_NULL;	//!< Describes the model type
	int coreStates = 6;		//!< The number of "core" states; these are computed in the simple EOM function; default is 6
	int stmStates = 36;		//!< The number of states used to store the STM; will always be 36
	int extraStates = 0;		//!< The number of extra states stored after the core states and STM states; default is zero.

	void copyMe(const tpat_model&);
};

#endif
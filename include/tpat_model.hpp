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

#include "tpat_constraint.hpp"

#include <vector>

// Forward Declarations
class tpat_constraint;
class tpat_event;
class tpat_nodeset;
class tpat_traj;
class tpat_sys_data;
struct iterationData;

/**
 *	@brief A base class that defines the behavior of a dynamical model and provides
 *	functions used in simulation and correction algorithms
 */
class tpat_model{

public:
	/**
	 *	@brief Describes the type of dynamic model; used for easy identification
	 */
	enum dynamicModel_t{
		MODEL_NULL,			//!< Undefined type
		MODEL_CR3BP,		//!< Circular Restricted, 3-Body Problem
		MODEL_CR3BP_LTVP,	//!< Circular Restricted, 3-Body Problem with Velocity-Pointing Low Thrust (constant)
		MODEL_BCR4BPR		//!< Bi-Circular Restricted, 4-Body Problem in a rotating reference frame
	};

	/**
	 *	@brief A function pointer to an EOM function
	 *
	 *	All EOM functions must take this form to work with GSL's 
	 *	integrators.
	 */
	typedef int (*eom_fcn)(double, const double[], double[], void*);

	// *structors
	tpat_model(dynamicModel_t);
	tpat_model(const tpat_model&);
	virtual ~tpat_model();

	// Operators
	tpat_model& operator =(const tpat_model&);

	/**
	 *	@brief Takes an input state and time and saves the data to the trajectory
	 *	@param y an array containing the core state and any extra states integrated
	 *	by the EOM function, including STM elements.
	 *	@param t the time at the current integration state
	 *	@param traj a pointer to the trajectory we should store the data in
	 */
	virtual void saveIntegratedData(double *y, double t, tpat_traj* traj) = 0;
	
	/**
	 *  @brief Use a correction algorithm to accurately locate an event crossing
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
	 *  @return wether or not the event has been located. If it has, a new point
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

	/**
	 *	@brief Compute the positions of all primaries
	 *
	 *	@param t the time at which the computations occur (only important for non-autonomous systems)
	 *	@param sysData object describing the specific system
	 *	@return an n x 3 vector (row-major order) containing the positions of
	 *	n primaries; each row is one position vector in non-dimensional units
	 */
	virtual std::vector<double> getPrimPos(double t, tpat_sys_data *sysData) = 0;

	/**
	 *	@brief Compute the velocities of all primaries
	 *
	 *	@param t the time at which the computations occur (only important for non-autonomous systems)
	 *	@param sysData object describing the specific system
	 *	@return an n x 3 vector (row-major order) containing the velocities of
	 *	n primaries; each row is one velocity vector in non-dimensional units
	 */
	virtual std::vector<double> getPrimVel(double t, tpat_sys_data *sysData) = 0;

	virtual void corrector_initDesignVec(iterationData*, tpat_nodeset*);
	virtual void corrector_createContCons(iterationData*, tpat_nodeset*);
	virtual void corrector_getSimICs(iterationData*, tpat_nodeset*, int, double*, double*, double*);
	virtual void corrector_applyConstraint(iterationData*, tpat_constraint, int);
	virtual void corrector_targetPosVelCons(iterationData*, tpat_constraint, int);
	virtual void corrector_targetExContCons(iterationData*, tpat_constraint, int);
	virtual void corrector_targetState(iterationData*, tpat_constraint, int);
	virtual void corrector_targetMatchAll(iterationData*, tpat_constraint, int);
	virtual void corrector_targetMatchCust(iterationData*, tpat_constraint, int);
	virtual void corrector_targetDist(iterationData*, tpat_constraint, int);
	virtual void corrector_targetDeltaV(iterationData*t, tpat_constraint, int);
	virtual void corrector_targetTOF(iterationData*, tpat_constraint, int);

	/**
	 *  @brief Take the final, corrected free variable vector <tt>X</tt> and create an output 
	 *  nodeset
	 *
	 *  If <tt>findEvent</tt> is set to true, the
	 *  output nodeset will contain extra information for the simulation engine to use. Rather than
	 *  returning only the position and velocity states, the output nodeset will contain the STM 
	 *  and other values for the final node; this information will be appended to the extraParameter
	 *  vector in the final node.
	 *
	 *  @param it an iteration data object containing all info from the corrections process
	 *	@param nodes_in a pointer to the original, uncorrected nodeset
	 *	@param findEvent whether or not this correction process is locating an event
	 *
	 *  @return a pointer to a nodeset containing the corrected nodes
	 */
	virtual tpat_nodeset* corrector_createOutput(iterationData* it, tpat_nodeset *nodes_in, bool findEvent) = 0;

	// Set and Get Functions
	int getCoreStateSize() const;
	int getSTMStateSize() const;
	int getExtraStateSize() const;
	bool supportsCon(tpat_constraint::constraint_t) const;
protected:
	dynamicModel_t modelType = MODEL_NULL;	//!< Describes the model type
	int coreStates = 6;		//!< The number of "core" states; these are computed in the simple EOM function; default is 6
	int stmStates = 36;		//!< The number of states used to store the STM; will always be 36
	int extraStates = 0;	//!< The number of extra states stored after the core states and STM states; default is zero.

	/** A vector containing the all the types of constraints this model supports */
	std::vector<tpat_constraint::constraint_t> allowedCons {tpat_constraint::NONE, tpat_constraint::STATE,
		tpat_constraint::MATCH_ALL, tpat_constraint::MATCH_CUST,
		tpat_constraint::DIST, tpat_constraint::MIN_DIST, tpat_constraint::MAX_DIST,
		tpat_constraint::MAX_DELTA_V, tpat_constraint::DELTA_V,
		tpat_constraint::TOF,
		tpat_constraint::CONT_PV, tpat_constraint::CONT_EX};

	void copyMe(const tpat_model&);
};

#endif
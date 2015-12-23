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

#include "tpat_constants.hpp"
#include "tpat_constraint.hpp"
#include "tpat_event.hpp"

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
 *
 *	This class provides the flexibility that allows the engine objects to operate on
 *	any dynamical model. This class is abstract and cannot be instantiated as an object,
 *	and, as such, only provides a framework and some common methods for other derived
 *	methods.
 *
 *	The tpat_simulation_engine and tpat_correction_engine make heavy use of dynamic 
 *	models to generalize their code. This allows the developer to easily implement
 *	new dynamic models for use in the simulator and corrector without making major
 *	modifications to their code.
 *
 *	The simulation engine calls the getSimpleEOM_fcn() and getFullEOM_fcn() to 
 *	obtain function pointers to the equations of motion for a dynamic model. It
 *	also calls the sim_locateEvent() and sim_saveIntegratedData() functions.
 *
 *	The correction engine calls functions defined in the dynamic model to 
 *	populate the design vector, constraint vector, and Jacobi matrix. Again,
 *	placing these functions in a model class allows the developer to add new
 *	models to the system without making huge changes to the complex correction
 *	algorithm.
 *
 *	In addition to useful functions, the dynamic model stores information about
 *	the number of states the equations of motion compute, and also contains a 
 *	list of all constraints the model supports in a corrections environment. 
 *	Derived classes can modify these variables to their needs.
 *
 *	@author Andrew Cox
 *	@version August 3, 2015
 *	@copyright GNU GPL v3.0
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

	virtual void multShoot_initDesignVec(iterationData*, tpat_nodeset*);
	virtual void multShoot_createContCons(iterationData*, tpat_nodeset*);
	virtual void multShoot_getSimICs(iterationData*, tpat_nodeset*, int, double*, double*, double*);
	virtual void multShoot_applyConstraint(iterationData*, tpat_constraint, int);
	virtual double multShoot_getSlackVarVal(iterationData*, tpat_constraint);

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
	virtual tpat_nodeset* multShoot_createOutput(iterationData* it, tpat_nodeset *nodes_in, bool findEvent) = 0;

	/**
	 *  @brief Use a correction algorithm to accurately locate an event crossing
	 *
	 *  The simulation engine calls this function if and when it determines that an event 
	 *  has been crossed. To accurately locate the event, we employ differential corrections
	 *  and find the exact event occurence in space and time.
	 *
	 *  @param event the event we're looking for
	 *  @param traj a pointer to the trajectory the event should occur on
	 *  @param ic the core state vector for this system
	 *  @param t0 non-dimensional time at the beginning of the search arc
	 *  @param tof the time-of-flight for the arc to search over
	 *  @param verbose whether or not we should be verbose with output messages
	 *
	 *  @return wether or not the event has been located. If it has, a new point
	 *  has been appended to the trajectory's data vectors.
	 */
	virtual bool sim_locateEvent(tpat_event event, tpat_traj *traj,
    	double *ic, double t0, double tof, verbosity_t verbose) = 0;

	/**
	 *	@brief Takes an input state and time and saves the data to the trajectory
	 *	@param y an array containing the core state and any extra states integrated
	 *	by the EOM function, including STM elements.
	 *	@param t the time at the current integration state
	 *	@param traj a pointer to the trajectory we should store the data in
	 */
	virtual void sim_saveIntegratedData(double *y, double t, tpat_traj* traj) = 0;

	// Set and Get Functions
	int getCoreStateSize() const;
	int getSTMStateSize() const;
	int getExtraStateSize() const;
	bool supportsCon(tpat_constraint::constraint_t) const;
	bool supportsEvent(tpat_event::event_t) const;
	
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
		tpat_constraint::TOF, tpat_constraint::APSE,
		tpat_constraint::CONT_PV, tpat_constraint::CONT_EX};

	/** A vector containing all the types of events this model supports */
	std::vector<tpat_event::event_t> allowedEvents {tpat_event::NONE, tpat_event::XY_PLANE, tpat_event::XZ_PLANE,
		tpat_event::YZ_PLANE, tpat_event::CRASH, tpat_event::APSE, tpat_event::DIST};

	void copyMe(const tpat_model&);
	virtual void multShoot_targetApse(iterationData*, tpat_constraint, int);
	virtual void multShoot_targetDeltaV(iterationData*, tpat_constraint, int);
	virtual double multShoot_targetDeltaV_compSlackVar(iterationData*, tpat_constraint);
	virtual void multShoot_targetDist(iterationData*, tpat_constraint, int);
	virtual double multShoot_targetDist_compSlackVar(iterationData*, tpat_constraint);
	virtual void multShoot_targetExContCons(iterationData*, tpat_constraint, int);
	virtual void multShoot_targetMatchAll(iterationData*, tpat_constraint, int);
	virtual void multShoot_targetMatchCust(iterationData*, tpat_constraint, int);
	virtual void multShoot_targetPosVelCons(iterationData*, tpat_constraint, int);
	virtual void multShoot_targetState(iterationData*, tpat_constraint, int);
	virtual void multShoot_targetTOF(iterationData*, tpat_constraint, int);
};

#endif
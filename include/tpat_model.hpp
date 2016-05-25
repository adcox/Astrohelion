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

#include "tpat.hpp"
 
#include "tpat_constants.hpp"
#include "tpat_constraint.hpp"
#include "tpat_event.hpp"
#include "tpat_sys_data.hpp"
 
#include <vector>

// Forward Declarations
class TPAT_Constraint;
class TPAT_Event;
class TPAT_Nodeset;
class TPAT_Traj;
class TPAT_Sys_Data;
class TPAT_MultShoot_Data;

/**
 *  @brief Container for EOM parameters
 *  @details At the current time, this object stores only the system data object pointer.
 *  Since the GSL functions demand a null pointer and the system data pointers owned by
 *  most objects are const-modified, this object serves as a non-const wrapper for 
 *  the system data pointers.
 *  
 *  @param sys A system data object
 *  @return a reference to this struct
 */
struct eomParamStruct{
	/**
	 *  @brief Construct an EOM Parameter structure
	 * 
	 *  @param sys a pointer to a system data object
	 */
	eomParamStruct(const TPAT_Sys_Data *sys) : sysData(sys) {}
	const TPAT_Sys_Data *sysData;	//!< Pointer to a system data object that will be passed into the EOMs
};

/**
 *	@brief A base class that defines the behavior of a dynamical model and provides
 *	functions used in simulation and correction algorithms
 *
 *	This class provides the flexibility that allows the engine objects to operate on
 *	any dynamical model. This class is abstract and cannot be instantiated as an object,
 *	and, as such, only provides a framework and some common methods for other derived
 *	methods.
 *
 *	The TPAT_Sim_Engine and TPAT_Correction_Engine make heavy use of dynamic 
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
class TPAT_Model : public TPAT{

public:
	/**
	 *	@brief Describes the type of dynamic model; used for easy identification
	 */
	enum tpat_dynamicsModel_tp{
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
	TPAT_Model(tpat_dynamicsModel_tp);
	TPAT_Model(const TPAT_Model&);
	virtual ~TPAT_Model();

	// Operators
	TPAT_Model& operator =(const TPAT_Model&);

	/**
	 *	@brief Retrieve a pointer to the EOM function that computes derivatives
	 *	for only the core states (i.e. simple)
	 *
	 *	The EOM function must be a non-member function; we store them in the 
	 *	tpat_calculations.cpp file
	 */
	virtual eom_fcn getSimpleEOM_fcn() const = 0;

	/**
	 *	@brief Retrieve a pointer to the EOM function that computes derivatives
	 *	for all states (i.e. full)
	 *
	 *	The EOM function must be a non-member function; we store them in the 
	 *	tpat_calculations.cpp file
	 */
	virtual eom_fcn getFullEOM_fcn() const = 0;

	/**
	 *	@brief Compute the positions of all primaries
	 *
	 *	@param t the time at which the computations occur (only important for non-autonomous systems)
	 *	@param sysData object describing the specific system
	 *	@return an n x 3 vector (row-major order) containing the positions of
	 *	n primaries; each row is one position vector in non-dimensional units
	 */
	virtual std::vector<double> getPrimPos(double t, const TPAT_Sys_Data *sysData) const = 0;

	/**
	 *	@brief Compute the velocities of all primaries
	 *
	 *	@param t the time at which the computations occur (only important for non-autonomous systems)
	 *	@param sysData object describing the specific system
	 *	@return an n x 3 vector (row-major order) containing the velocities of
	 *	n primaries; each row is one velocity vector in non-dimensional units
	 */
	virtual std::vector<double> getPrimVel(double t, const TPAT_Sys_Data *sysData) const = 0;

	virtual double getRDot(int, double, const double*, const TPAT_Sys_Data*) const;

	/**
	 *  @brief Do any model-specific initializations for the TPAT_MultShoot_Data object
	 *  @param it a pointer to the TPAT_MultShoot_Data object for the multiple shooting process
	 */
	virtual void multShoot_initIterData(TPAT_MultShoot_Data *it) const = 0;

	virtual void multShoot_initDesignVec(TPAT_MultShoot_Data*, const TPAT_Nodeset*) const;
	virtual void multShoot_scaleDesignVec(TPAT_MultShoot_Data*, const TPAT_Nodeset*) const;
	virtual void multShoot_createContCons(TPAT_MultShoot_Data*, const TPAT_Nodeset*) const;
	virtual void multShoot_getSimICs(const TPAT_MultShoot_Data*, const TPAT_Nodeset*, int, double*, double*, double*) const;
	virtual void multShoot_applyConstraint(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	virtual double multShoot_getSlackVarVal(const TPAT_MultShoot_Data*, TPAT_Constraint)const ;

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
	 *	@param nodesOut pointer to a nodeset object that will contain the results of
	 *	the shooting process
	 *  @return a pointer to a nodeset containing the corrected nodes
	 */
	virtual void multShoot_createOutput(const TPAT_MultShoot_Data* it, const TPAT_Nodeset *nodes_in, bool findEvent, TPAT_Nodeset *nodesOut) const = 0;

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
	virtual bool sim_locateEvent(TPAT_Event event, TPAT_Traj *traj,
    	const double *ic, double t0, double tof, TPAT_Verbosity_Tp verbose) const = 0;

	/**
	 *	@brief Takes an input state and time and saves the data to the trajectory
	 *	@param y an array containing the core state and any extra states integrated
	 *	by the EOM function, including STM elements.
	 *	@param t the time at the current integration state
	 *	@param traj a pointer to the trajectory we should store the data in
	 */
	virtual void sim_saveIntegratedData(const double *y, double t, TPAT_Traj* traj) const = 0;

	// Set and Get Functions
	int getCoreStateSize() const;
	int getSTMStateSize() const;
	int getExtraStateSize() const;
	bool supportsCon(TPAT_Constraint_Tp) const;
	bool supportsEvent(TPAT_Event_Tp) const;
	
protected:
	tpat_dynamicsModel_tp modelType = MODEL_NULL;	//!< Describes the model type
	int coreStates = 6;		//!< The number of "core" states; these are computed in the simple EOM function; default is 6
	int stmStates = 36;		//!< The number of states used to store the STM; will always be 36
	int extraStates = 0;	//!< The number of extra states stored after the core states and STM states; default is zero.

	/** A vector containing the all the types of constraints this model supports */
	std::vector<TPAT_Constraint_Tp> allowedCons {TPAT_Constraint_Tp::NONE, TPAT_Constraint_Tp::STATE,
		TPAT_Constraint_Tp::MATCH_ALL, TPAT_Constraint_Tp::MATCH_CUST,
		TPAT_Constraint_Tp::DIST, TPAT_Constraint_Tp::MIN_DIST, TPAT_Constraint_Tp::MAX_DIST,
		TPAT_Constraint_Tp::MAX_DELTA_V, TPAT_Constraint_Tp::DELTA_V,
		TPAT_Constraint_Tp::TOF, TPAT_Constraint_Tp::APSE,
		TPAT_Constraint_Tp::CONT_PV, TPAT_Constraint_Tp::CONT_EX,
		TPAT_Constraint_Tp::SEG_CONT_PV, TPAT_Constraint_Tp::SEG_CONT_EX};

	/** A vector containing all the types of events this model supports */
	std::vector<TPAT_Event_Tp> allowedEvents {TPAT_Event_Tp::NONE, TPAT_Event_Tp::XY_PLANE, TPAT_Event_Tp::XZ_PLANE,
		TPAT_Event_Tp::YZ_PLANE, TPAT_Event_Tp::CRASH, TPAT_Event_Tp::APSE, TPAT_Event_Tp::DIST};

	void copyMe(const TPAT_Model&);
	virtual void multShoot_targetApse(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	virtual void multShoot_targetDeltaV(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	virtual double multShoot_targetDeltaV_compSlackVar(const TPAT_MultShoot_Data*, TPAT_Constraint) const;
	virtual void multShoot_targetDist(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	virtual double multShoot_targetDist_compSlackVar(const TPAT_MultShoot_Data*, TPAT_Constraint) const;
	virtual void multShoot_targetCont_Ex(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	virtual void multShoot_targetCont_Ex_Seg(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	virtual void multShoot_targetMatchAll(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	virtual void multShoot_targetMatchCust(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	virtual void multShoot_targetCont_PosVel(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	virtual void multShoot_targetCont_PosVel_Seg(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	virtual void multShoot_targetState(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
	virtual void multShoot_targetTOF(TPAT_MultShoot_Data*, TPAT_Constraint, int) const;
};

#endif
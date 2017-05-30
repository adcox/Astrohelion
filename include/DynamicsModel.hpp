/**
 *  \file DynamicsModel.hpp
 *	\brief 
 *	
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of Astrohelion
 *
 *  Astrohelion is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Astrohelion is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Astrohelion.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "Core.hpp"
 
#include "Common.hpp"
#include "Constraint.hpp"
#include "Event.hpp"
#include "SysData.hpp"
 
#include <vector>


namespace astrohelion{

// Forward Declarations
class Constraint;
class Event;
class Node;
class Arcset;
class Arcset;
class Segment;
class SysData;
class MultShootData;

/**
 *  \ingroup model
 *  \brief Container for EOM parameters
 *  \details At the current time, this object stores only the system data object pointer.
 *  Since the GSL functions demand a null pointer and the system data pointers owned by
 *  most objects are const-modified, this object serves as a non-const wrapper for 
 *  the system data pointers.
 *  
 *  \param sys A system data object
 *  \return a reference to this struct
 */
struct EOM_ParamStruct{

	/**
	 *  \brief Construct an EOM Parameter structure
	 * 
	 *  \param sys a pointer to a system data object
	 *  \param lawID an ID for the control law to use during the integration
	 */
	EOM_ParamStruct(const SysData *sys, ControlLaw *law) : pSysData(sys), pCtrlLaw(law) {}
	
	const SysData *pSysData;			//!< Pointer to a system data object that will be passed into the EOMs
	ControlLaw *pCtrlLaw;				//!< Pointer to a control law object that will be passed into the EOMs
};

/**
 *	\brief Describes the type of dynamic model; used for easy identification
 */
enum class DynamicsModel_tp{
	MODEL_NULL,			//!< Undefined type
	MODEL_2BP,			//!< Relative 2-Body Problem
	MODEL_CR3BP,		//!< Circular Restricted, 3-Body Problem
	MODEL_CR3BP_LT,		//!< Circular Restricted, 3-Body Problem with Velocity-Pointing Low Thrust (constant)
	MODEL_BCR4BPR		//!< Bi-Circular Restricted, 4-Body Problem in a rotating reference frame
};

/**
 *	\brief A base class that defines the behavior of a dynamical model and provides
 *	functions used in simulation and correction algorithms
 *
 *	This class provides the flexibility that allows the engine objects to operate on
 *	any dynamical model. This class is abstract and cannot be instantiated as an object,
 *	and, as such, only provides a framework and some common methods for other derived
 *	methods.
 *
 *	The SimEngine and MultShootEngine make heavy use of dynamic 
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
 *	\author Andrew Cox
 *	\version August 3, 2015
 *	\copyright GNU GPL v3.0
 */
class DynamicsModel : public Core{

public:
	/**
	 *	\brief A function pointer to an EOM function
	 *	\details All EOM functions must take this form to work with GSL's 
	 *	integrators.
	 *	
	 *	\param t time associated with the current integration step; this value is passed
	 *	in by the integrator
	 *	\param q full state vector; the dimension of this vector is set by SimEngine
	 *	and is equal to <code>coreDim + stmStates + extraDim</code>; this vector is passed
	 *	in by the integrator
	 *	\param qdot full state derivative vector, same dimension as <code>q</code>. This vector
	 *	must be initialized by the function (i.e., set all values to zero or another value), or
	 *	invalid memory access will occur within GSL's integration methods
	 *	\param params a pointer to a set of extra parameters required for integration. All
	 *	of the EOM functions in this software take a pointer to a SysData object
	 *	consistent with the DynamicsModel as this parameter. The parameter set can be modified 
	 *	between integration steps (i.e., change model parameters), but the ode functions must be reset
	 *	via <code>gsl_odeiv2_driver_reset</code>, <code>gsl_odeiv2_evolve_reset</code>, or
	 *	<code>gsl_odeiv2_step_reset</code> before continuing with an updated parameter set
	 */
	typedef int (*eom_fcn)(double t, const double q[], double qdot[], void *params);

	/**
	 *  \name *structors
	 *  \{
	 */
	DynamicsModel(DynamicsModel_tp);
	DynamicsModel(const DynamicsModel&);
	virtual ~DynamicsModel();
	//\}

	/**
	 *  \name Operators
	 *  \{
	 */
	DynamicsModel& operator =(const DynamicsModel&);
	//\}

	/**
	 *  \name Core Functions
	 *  \{
	 */

	/**
	 *	\brief Retrieve a pointer to the EOM function that computes derivatives
	 *	for only the core states (i.e. simple)
	 *
	 *	The EOM function must be a non-member function; we store them in the 
	 *	Calculations.cpp file
	 */
	virtual eom_fcn getSimpleEOM_fcn() const = 0;

	/**
	 *	\brief Retrieve a pointer to the EOM function that computes derivatives
	 *	for all states (i.e. full)
	 *
	 *	The EOM function must be a non-member function; we store them in the 
	 *	Calculations.cpp file
	 */
	virtual eom_fcn getFullEOM_fcn() const = 0;

	/**
	 *	\brief Compute the positions of all primaries
	 *
	 *	\param t the time at which the computations occur (only important for non-autonomous systems)
	 *	\param sysData object describing the specific system
	 *	\return an n x 3 vector (row-major order) containing the positions of
	 *	n primaries; each row is one position vector in non-dimensional units
	 */
	virtual std::vector<double> getPrimPos(double t, const SysData *pSysData) const = 0;

	/**
	 *  \brief Compute the position of a specified primary
	 *  \details This is the faster alternative to getPrimPos(t, pSysData).
	 * 
	 *  \param t Nondimensional time
	 *  \param pSysData pointer to system data object
	 *  \param pIx Index of the primary; a value of -1 will return the positions of all primaries,
	 *  in order of largest to smallest mass
	 *  \param pos An array to store the primary position(s) in with all elements initialized to zero.
	 *  For a single primary position, the array must have at least three elements allocated. For all 
	 *  primaries (i.e., pIx = -1), the array must have n*3 elements allocated where n is the number 
	 *  of primaries.
	 */
	virtual void getPrimPos(double t, const SysData *pSysData, int pIx, double *pos) const = 0;

	/**
	 *	\brief Compute the velocities of all primaries
	 *
	 *	\param t the time at which the computations occur (only important for non-autonomous systems)
	 *	\param sysData object describing the specific system
	 *	\return an n x 3 vector (row-major order) containing the velocities of
	 *	n primaries; each row is one velocity vector in non-dimensional units
	 */
	virtual std::vector<double> getPrimVel(double t, const SysData *pSysData) const = 0;

	/**
	 *  \brief Compute the velocity of a specified primary
	 *  \details This is the faster alternative to getPrimVel(t, pSysData).
	 * 
	 *  \param t Nondimensional time
	 *  \param pSysData pointer to system data object
	 *  \param pIx Index of the primary; a value of -1 will return the velocities of all primaries,
	 *  in order of largest to smallest mass
	 *  \param vel An array to store the primary velocity(s) in with all elements initialized to zero. 
	 *  For a single primary velocity, the array must have at least three elements allocated. For all 
	 *  primaries (i.e., pIx = -1), the array must have n*3 elements allocated where n is the number 
	 *  of primaries.
	 */
	virtual void getPrimVel(double t, const SysData *pSysData, int pIx, double *vel) const = 0;

	virtual double getRDot(int, double, const double*, const SysData*) const;

	virtual std::vector<double> getStateDeriv(double t, std::vector<double> state, EOM_ParamStruct *params) const = 0;

	//\}

	/**
	 *  \name Multiple Shooting Support Functions
	 *  \{
	 */

	/**
	 *  \brief Do any model-specific initializations for the MultShootData object
	 *  \param it a pointer to the MultShootData object for the multiple shooting process
	 */
	virtual void multShoot_initIterData(MultShootData *it) const = 0;

	virtual void multShoot_initDesignVec(MultShootData*) const;
	virtual void multShoot_createContCons(MultShootData*) const;
	virtual void multShoot_getSimICs(const MultShootData*, int, double*, double*, double*, double*) const;
	virtual void multShoot_applyConstraint(MultShootData*, const Constraint&, int) const;
	virtual double multShoot_getSlackVarVal(const MultShootData*, const Constraint&)const ;

	/**
	 *  \brief Take the final, corrected free variable vector <tt>X</tt> and create an output 
	 *  nodeset
	 *
	 *  If <tt>findEvent</tt> is set to true, the
	 *  output nodeset will contain extra information for the simulation engine to use. Rather than
	 *  returning only the position and velocity states, the output nodeset will contain the STM 
	 *  and other values for the final node; this information will be appended to the extraParameter
	 *  vector in the final node.
	 *
	 *  \param it an iteration data object containing all info from the corrections process
	 *  
	 *  \return a pointer to a nodeset containing the corrected nodes
	 */
	virtual void multShoot_createOutput(const MultShootData* it) const = 0;
	//\}

	/**
	 *  \name Simulation Support Functions
	 *  \{
	 */

	virtual std::vector<Event> sim_makeDefaultEvents(const SysData*) const;
	virtual int sim_addNode(Node &node, const double *y, double t, Arcset* traj, EOM_ParamStruct *params, Event_tp tp) const;
	virtual int sim_addSeg(Segment &seg, const double *y, double t, Arcset* traj, EOM_ParamStruct *params) const;

	//\}

	/**
	 *  \name Set and Get Functions
	 *  \{
	 */
	unsigned int getCoreStateSize() const;
	unsigned int getExtraStateSize() const;
	bool supportsCon(Constraint_tp) const;
	virtual bool supportsControl(const ControlLaw*) const;
	bool supportsEvent(Event_tp) const;
	//\}

	/**
	 *  \name Utility Functions
	 *  \{
	 */
	virtual ControlLaw* createControlLaw() const;
	//\}
protected:
	DynamicsModel_tp modelType = DynamicsModel_tp::MODEL_NULL;	//!< Describes the model type
	unsigned int coreDim = 6;		//!< The number of "core" states; these are computed in the simple EOM function; default is 6; STM is an nxn matrix with n = coreDim
	unsigned int extraDim = 0;	//!< The number of extra states stored after the core states and STM states; default is zero.

	/** A vector containing the all the types of constraints this model supports */
	std::vector<Constraint_tp> allowedCons {Constraint_tp::NONE,
		Constraint_tp::STATE,
		Constraint_tp::MATCH_ALL,
		Constraint_tp::MATCH_CUST,
		Constraint_tp::DIST,
		Constraint_tp::MIN_DIST,
		Constraint_tp::MAX_DIST,
		Constraint_tp::MAX_DELTA_V,
		Constraint_tp::DELTA_V,
		Constraint_tp::TOF,
		Constraint_tp::APSE,
		Constraint_tp::CONT_CTRL,
		Constraint_tp::CONT_PV,
		Constraint_tp::CONT_EX,
		Constraint_tp::SEG_CONT_PV,
		Constraint_tp::SEG_CONT_EX,
		Constraint_tp::RM_STATE,
		Constraint_tp::ENDSEG_STATE};

	/** A vector containing all the types of events this model supports */
	std::vector<Event_tp> allowedEvents {Event_tp::NONE,
		Event_tp::SIM_TOF,
		Event_tp::SIM_COMPTIME,
		Event_tp::SIM_ERR,
		Event_tp::XY_PLANE,
		Event_tp::XZ_PLANE,
		Event_tp::YZ_PLANE,
		Event_tp::CRASH,
		Event_tp::APSE,
		Event_tp::DIST};

	void copyMe(const DynamicsModel&);

	/**
	 *  \name Multiple Shooting Support Functions
	 *  \{
	 */
	virtual void multShoot_targetApse(MultShootData*, const Constraint&, int) const;
	virtual void multShoot_targetDeltaV(MultShootData*, const Constraint&, int) const;
	virtual double multShoot_targetDeltaV_compSlackVar(const MultShootData*, const Constraint&) const;
	virtual void multShoot_targetDist(MultShootData*, const Constraint&, int) const;
	virtual double multShoot_targetDist_compSlackVar(const MultShootData*, const Constraint&) const;
	virtual void multShoot_targetCont_Ctrl(MultShootData*, const Constraint&, int) const;
	virtual void multShoot_targetCont_Ex(MultShootData*, const Constraint&, int) const;
	virtual void multShoot_targetCont_Ex_Seg(MultShootData*, const Constraint&, int) const;
	virtual void multShoot_targetMatchAll(MultShootData*, const Constraint&, int) const;
	virtual void multShoot_targetMatchCust(MultShootData*, const Constraint&, int) const;
	virtual void multShoot_targetCont_State(MultShootData*, const Constraint&, int) const;
	virtual void multShoot_targetCont_State_Seg(MultShootData*, const Constraint&, int) const;
	virtual void multShoot_targetState(MultShootData*, const Constraint&, int) const;
	virtual void multShoot_targetState_endSeg(MultShootData*, const Constraint&, int) const;
	virtual void multShoot_targetTOF(MultShootData*, const Constraint&, int) const;
	//\}
};


}// END of Astrohelion namespace
/**
 *  @file SimEngine.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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

#include <gsl/gsl_odeiv2.h>
#include <vector>

#include "Core.hpp"
#include "Engine.hpp"

#include "ControlLaw.hpp"
#include "DynamicsModel.hpp"
#include "Event.hpp"
#include "Arcset.hpp"


namespace astrohelion{

// Forward declarations
class SysData;

/**
 *  @brief DEPRECATED Wrapper object to use GSL EOM functions with Boost integrators
 */
class boost_eom_wrapper{
	/**
	 *  @brief Definition of a function that the boost integrators use
	 * 
	 *  @param t integration time
	 *  @param q state vector
	 *  @param qdot state time-derivative vector
	 *  @param params void pointer that stores parameters relevant to the integration
	 */
	typedef int (*eom_fcn)(double t, const double q[], double qdot[], void *params);
	
	eom_fcn gslFcn;		//!< The GSL-ready function that evaluates equations of motion
	EOM_ParamStruct *eomParams = nullptr;	//!< Parameters to be passed to the GSL EOM function

public:
	/**
	 *  @brief DEPRECATED Construct an EOM wrapper for 
	 *  @details 
	 * 
	 *  @param f 
	 *  @param params 
	 */
	boost_eom_wrapper(eom_fcn f, EOM_ParamStruct *params) : gslFcn(f){
		eomParams = params;
	}//================================================

	/**
	 *  @brief Evaluate the equations of motion.
	 *  @details This function is called by the Boost integrator. 
	 * 
	 *  @param q Reference to the current state
	 *  @param qdot Reference to a vector in which the state time-derivative is stored
	 *  @param t time
	 */
	void operator() (const std::vector<double> &q, std::vector<double> &qdot, const double t){
		qdot.assign(q.size(), 0);
		gslFcn(t, &(q.front()), &(qdot.front()), eomParams);	// Evaluate the EOMs from the GSL function
	}//================================================
};// END OF BOOST_EOM_WRAPPER--------------------------------------------------


/**
 *  @brief DEPRECATED Wrapper object to use already-written functions to save trajectory data from 
 *  Boost integrators.
 *  @details This wrapper does not handle any event detection or other, more advanced, 
 *  observer behaviors.
 *  
 */
class boost_observer{
	Arcset *traj;		//!< Trajectory being integrated
	const DynamicsModel *model;	//!< DynamicsModel associated with the problem
	EOM_ParamStruct *eomParams = nullptr;	//!< Data that is passed to the EOMs

public:
	/**
	 *  @brief DEPRECATED
	 * 
	 *  @param m 
	 *  @param t 
	 *  @param params 
	 */
	boost_observer(const DynamicsModel *m, Arcset *t, EOM_ParamStruct *params) : traj(t), model(m){
		eomParams = params;
	}//================================================

	/**
	 *  @brief Observe the propagation.
	 *  @details This function is called by the Boost integration library after each
	 *  step. Capabilities such as saving data and event detection should be handled here.
	 * 
	 *  @param q reference to the state vector
	 *  @param t time at the state
	 */
	void operator() (const std::vector<double> &q, double t){
		(void) q;
		(void) t;
		// model->sim_saveIntegratedData(&(q.front()), t, traj, eomParams);
	}//================================================
};// END OF BOOST_OBSERVER------------------------------------------------------

/**
 *  @brief Classify integrator types
 */
enum class Integ_tp{
	RKCK,			//!< Explicit embedded Runge-Kutta Cash-Karp (4,5); variable step propagations
	RK8PD,			//!< Explicit embedded Runge-Kutta Dormance-Prince (8,9); variable step propagations
	MSADAMS,		//!< Variable coefficient linear multistep Adams method in Nordisieck form; uses explicit Adams-Bashforth (predictor) and implicit Adams-Moulton (corrector); fixed step propagations
	BOOST_RKCK,		//!< Explicit Runge-Kutta Cash-Karp (4,5); variable step propagations. NOTE: limited capabilities (no event detection); "beta" testing
	BOOST_RKF,		//!< Explicit Runge-Kutta-Felburg (7,8); variable step propagations. NOTE: limited capabilities (no event detection); "beta" testing
	BOOST_BS 		//!< Explicit Bulirsch-Stoer; variable step propagations. NOTE: limited capabilities (no event detection); "beta" testing
};

/**
 *	\ingroup engine
 *	@brief Performs numerical integration on any system type and produces an
 *	Arcset object
 *
 *	The simulation engine is the workhorse object for the Core. It
 *	holds functions to integrate equations of motion and is called by the 
 *	`MultShootEngine` to compute arcs between nodes.
 *
 *	Creating a simulation engine is simple; it can either be instantiated with no
 *	arguments, or by specifying a system data object. Further settings can be applied
 *	via the "set and get" methods, and the `runSim()` method can be called to 
 *	integrate a trajectory for a specific amount of time. The integration will most likely
 *	run for the specified interval, but crash-detecting event functions are included by default
 *	and will end the integration if triggered. To run a simulation without these events,
 *	use the `setMakesCrashEvents()` to tell the engine to skip creating these events. Alternatively, 
 *	more event functions can be added to the simulation to end (or simply flag) the simulation
 *	at different event occurrences.
 *
 * 	Once the integration has completed, a trajectory object can be obtained by calling one of the
 *	`get_*Arcset()` functions. To reuse the simulation engine, call the `reset()`
 *  function and re-run the simulation with different initial conditions and time-of-flight.
 *
 *	<b>Events</b>
 *
 *	The user can tell the simulation to watch for certain types of events during the simulation; 
 *	if such an event occurs, the simulation can be made to end immediately, or record the 
 *	occurence and continue. The `addEvent()` function adds events to the simulation and
 *	the `clearEvents()` function removes all events from the simulation. In order for 
 *	the engine to locate events, the simulation must be run with a sufficient number of points.
 *	The engine compares points along the integrated path to the location of the event and flags
 *	an event occurence when the direction to the event changes. If too few points are used,
 *	the engine may never detect a change in direction.
 *	
 *	<b>Algorithms</b>
 *	
 *	Two integration algorithms are used in this engine to numerically integrate a model's
 *	equations of motion. One variable determines which integrator will be used: step size
 *	variability. If steps are allowed to vary, a Runge-Kutta Cash-Karp 4-5 method or a
 *	Runge-Kutta Dormand-Prince 8-9 method is used.
 *	If step size is fixed, then an Adams-Bashforth Adams-Moulton method is used. These 
 *	integrators performed best in an error analysis of Jacobi constant in an unscientific
 *	test comparison of multiple methods. Perhaps someday we'll actually do a propper comparison.
 *	
 *	Because the Adams-Bashforth method requires a driver object, it is not possible (to our knowledge)
 *	to retrieve intermediate points between the integration bounds. The driver will simply integrate from
 *	t0 to tf and return the final state. To obtain a more complete trajectory history, the simulation
 *	engine allows the user to specify the number of steps, which will be split evenly between the initial
 *	and final integration time limits. The Adams-Bashforth method is run between each step, effectively
 *	producing a series of states along a trajectory.
 *	
 *	The Runge-Kutta method, on the other hand, will return intermediate states because it is run using
 *	more detailed commands (we do not use a driver object). The user can allow step size to vary, and
 *	the integrator will return steps between t0 and tf at intervals that the integrator "likes". Typically
 *	this means smaller steps are taken near dynamically volatile areas (e.g. near a primary body) and larger
 *	steps are taken where the dynamics are less volatile.
 *	
 *	
 *	@author Andrew Cox
 *	@version June 1, 2015
 *	@copyright GNU GPL v3.0
 */
class SimEngine : public Core, public Engine{
	public:
		/**
		 *  \name *structors
		 *  \{
		 */
		SimEngine();
		SimEngine(const SimEngine&);	//copy constructor
		~SimEngine();
		//\}

		/**
		 *  \name Operators
		 *  \{
		 */
		SimEngine& operator =(const SimEngine&);
		//\}

		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		double getAbsTol() const;
		std::vector<Event> getEvents() const;
		int getNumSteps() const;
		double getRelTol() const;
		bool makesDefaultEvents() const;
		bool usesSimpleInt() const;
		bool usesRevTime() const;
		bool usesVarStepSize() const;
		
		void setAbsTol(double);
		void setFixStepInteg(Integ_tp);
		void setMakeDefaultEvents(bool);
		void setMaxCompTime(double);
		void setNumSteps(int);
		void setSimpleInt(bool);
		void setRelTol(double);
		void setRevTime(bool);
		void setVarStepInteg(Integ_tp);
		void setVarStepSize(bool);
		//\}

		/**
		 *  \name Analysis Functions
		 *  \{
		 */
		int addEvent(Event);

		// Assume t0, ctrl0, stm0
		void runSim(std::vector<double> ic, double tof, Arcset *arcset, ControlLaw *pLaw = nullptr);
		void runSim(const double* ic, double tof, Arcset *arcset, ControlLaw *pLaw = nullptr);
		
		// Assume ctrl0, stm0
		void runSim(std::vector<double> ic, double t0, double tof, Arcset *arcset, ControlLaw *pLaw = nullptr);
		void runSim(const double* ic, double t0, double tof, Arcset *arcset, ControlLaw *pLaw = nullptr);
		
		// Assume stm0
		void runSim(std::vector<double> ic, std::vector<double> ctrl0, double t0, double tof, Arcset *arcset, ControlLaw *pLaw);

		// No assumptions
		void runSim(std::vector<double> ic, std::vector<double> ctrl0, const MatrixXRd &stm, double t0, double tof, Arcset *arcset, ControlLaw *pLaw);
		void runSim(const double* ic, const double* ctrl0, const double* stm, double t0, double tof, Arcset *arcset, ControlLaw *pLaw);
		
		// No assumptions, final fcn before calling integrate()
		void runSim(const double *ic, const double *ctrl0, const double *stm, std::vector<double> tspan, Arcset*, ControlLaw *pLaw);

		// Assume t0, ctrl0, stm0
		void runSim_manyNodes(const double* ic, double tof, int numNodes, Arcset *arcset, ControlLaw *pLaw = nullptr);
		void runSim_manyNodes(std::vector<double> ic, double tof, int numNodes, Arcset *arcset, ControlLaw *pLaw = nullptr);

		// Assume ctrl0, stm0
		void runSim_manyNodes(const double* ic, double t0, double tof, int numNodes, Arcset *arcset, ControlLaw *pLaw = nullptr);
		void runSim_manyNodes(std::vector<double> ic, double t0, double tof, int numNodes, Arcset *arcset, ControlLaw *pLaw = nullptr);

		// Assume stm0
		void runSim_manyNodes(std::vector<double> ic, std::vector<double> ctrl0, double t0, double tof, int numNodes, Arcset *arcset, ControlLaw *pLaw);
		//\}
		
		/**
		 *  \name Utility Functions
		 *  \{
		 */
		void clearEvents();
		void reset();
		//\}

	private:
		/** Vector of events to consider during integration */
		std::vector<Event> events {};
		
		/** a pointer to a data structure with parameters for the EOM function */
		EOM_ParamStruct *eomParams = nullptr;

		/** Whether or not to run the simulation in reverse time */
		bool bRevTime = false;

		/** Whether or not to use variable step size when integrating */
		bool bVarStepSize = true;

		/** Whether or not to use "simple" integration; simple means no STM or other quantites are integrated
		 with the EOMs */
		bool bSimpleIntegration = false;

		/** Whether or not default events should be created for the simulation */
		bool bMakeDefaultEvents = true;

		/** Whether or not the default events have been created */
		bool bMadeDefaultEvents = false;

		/** Absolute tolerance for integrated data, units are same as integrated data */
		double absTol = 1e-12;

		/** Relative tolerance for integrated data, units are same as integrated data */
		double relTol = 1e-14;

		/** Initial guess for time step size */
		double dtGuess = 1e-6;

		/** Number of steps to take (only applies when varStepSize = false) */
		double numSteps = 1000;

		/** Number of seconds allowed for integration; set to -1 to allow unlimitted computation time */
		double maxCompTime = -1;

		/** 
		 * Maximum number of steps the fixed-step integration driver can take. This cap only limits
		 * integrations performed when bVarStepSize = false
		 */
		unsigned long int maxDriverSteps = 250000;

		/** Timestamp at integration start */
		time_t startTimestamp = 0;

		/** Integrator to use for variable step-size propagations; Default is RK8PD */
		Integ_tp varStep_integ = Integ_tp::RK8PD;
		
		/** Integrator to use for fixed step-size propagations; Default is MSADAMS */
		Integ_tp fixStep_integ = Integ_tp::MSADAMS;

		/**
		 *  \name Utility Functions
		 *  \{
		 */		
		void cleanEngine();
		void copyMe(const SimEngine&);
		void free_odeiv2(gsl_odeiv2_step*, gsl_odeiv2_control*, gsl_odeiv2_evolve*, gsl_odeiv2_driver*);
		void reportPropErrs(int, double);
		//\}

		/**
		 *  \name Analysis Functions
		 *  \{
		 */
		void createDefaultEvents(const SysData*);
		void integrate(const double*, const double*, const double*, const double*, unsigned int, Arcset*);
		bool locateEvents(const double*, double, Arcset*, int);
		bool locateEvent_multShoot(const double*, double, int, Arcset*);
		//\}
		
};


}// END of Astrohelion namespace
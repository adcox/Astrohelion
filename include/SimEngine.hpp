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
 *	Copyright 2016, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "DynamicsModel.hpp"
#include "Event.hpp"
#include "Traj.hpp"


namespace astrohelion{

// Forward declarations
class Traj_bc4bp;
class Traj_cr3bp;
class Traj_cr3bp_ltvp;
class SysData;

// struct gsl_odeiv2_step;
// struct gsl_odeiv2_control;
// struct gsl_odeiv2_evolve;
// struct gsl_odeiv2_driver;

/**
 *	@brief A small structure to store event occurrence records
 */
struct SimEventRecord{
public:
	/**
	 *	@brief Construct an event record
	 *	@param e the event index within the simulation engine event vector
	 *	@param s the step index; which step did this event occur at?
	 */
	SimEventRecord(int e, int s) : eventIx(e), stepIx(s) {}
	int eventIx;	//!< The index of the event (index from simulation engine vector of events)
	int stepIx;		//!< The index of the integration step the event occured at
};

enum class Integ_tp{
	RKCK,			//!< Explicit embedded Runge-Kutta Cash-Karp (4,5); variable step propagations
	RK8PD,			//!< Explicit embedded Runge-Kutta Dormance-Prince (8,9); variable step propagations
	MSADAMS			//!< Variable coefficient linear multistep Adams method in Nordisieck form; uses explicit Adams-Bashforth (predictor) and implicit Adams-Moulton (corrector); fixed step propagations
};

/**
 *	@ingroup engine
 *	@brief Performs numerical integration on any system type and produces an
 *	Traj object
 *
 *	The simulation engine is the workhorse object for the Core. It
 *	holds functions to integrate equations of motion and is called by the 
 *	<tt>CorrectionEngine</tt> to compute arcs between nodes.
 *
 *	Creating a simulation engine is simple; it can either be instantiated with no
 *	arguments, or by specifying a system data object. Further settings can be applied
 *	via the "set and get" methods, and the <tt>runSim()</tt> method can be called to 
 *	integrate a trajectory for a specific amount of time. The integration will most likely
 *	run for the specified interval, but crash-detecting event functions are included by default
 *	and will end the integration if triggered. To run a simulation without these events,
 *	call the <tt>clearEvents()</tt> function before running the simulation. Alternatively, 
 *	more event functions can be added to the simulation to end (or simply flag) the simulation
 *	at different event occurrences.
 *
 * 	Once the integration has completed, a trajectory object can be obtained by calling one of the
 *	<tt>get_*Traj()</tt> functions. To reuse the simulation engine, call the <tt>reset()</tt>
 *  function and re-run the simulation with different initial conditions and time-of-flight.
 *
 *	<b>Events</b>
 *
 *	The user can tell the simulation to watch for certain types of events during the simulation; 
 *	if such an event occurs, the simulation can be made to end immediately, or record the 
 *	occurence and continue. The <tt>addEvent()</tt> function adds events to the simulation and
 *	the <tt>clearEvents()</tt> function removes all events from the simulation. In order for 
 *	the engine to locate events, the simulation must be run with a sufficient number of points.
 *	The engine compares points along the integrated path to the location of the event and flags
 *	an event occurence when the direction to the event changes. If too few points are used,
 *	the engine may never detect a change in direction.
 *	
 *	<b>Algorithms</b>
 *	
 *	Two integration algorithms are used in this engine to numerically integrate a model's
 *	equations of motion. One variable determines which integrator will be used: step size
 *	variability. If steps are allowed to vary, a Runge-Kutta Cash-Karp 4-5 method is used.
 *	If step size is fixed, then an Adams-Bashforth Adams-Moulton method is used. These 
 *	integrators performed best in an error analysis of Jacobi constant in an unscientific
 *	test comparison of multiple methods. Perhaps someday we'll actually do a propper comparison.
 *	
 *	The Runge-Kutta method is similar to Matlab's ode45, which uses a Dormand-Prince 8-9 variation
 *	on the traditional Runge-Kutta method. The Cash-Karp variation used here should be nearly
 *	identical, but it performed better in the error analysis.
 *	
 *	The Adams-Bashforth method is similar to Matlab's ode113.
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
		// Constructors
		SimEngine();
		SimEngine(const SimEngine&);	//copy constructor
		
		//Destructor
		~SimEngine();

		//Operators
		SimEngine& operator =(const SimEngine&);

		/**
		 *  @name Set and Get Functions
		 *  @{
		 */
		void addEvent(Event);
		double getAbsTol() const;
		std::vector<Event> getEndEvents(Traj*) const;
		std::vector<Event> getEvents() const;
		std::vector<SimEventRecord> getEventRecords() const;
		int getNumSteps() const;
		double getRelTol() const;
		bool makesCrashEvents() const;
		bool usesSimpleInt() const;
		bool usesRevTime() const;
		bool usesVarStepSize() const;
		
		void setAbsTol(double);
		void setFixStepInteg(Integ_tp);
		void setMakeCrashEvents(bool);
		void setMaxCompTime(int);
		void setNumSteps(int);
		void setSimpleInt(bool);
		void setRelTol(double);
		void setRevTime(bool);
		void setVarStepInteg(Integ_tp);
		void setVarStepSize(bool);
		//@}

		/**
		 *  @name Simulation Methods
		 *  @{
		 */
		void runSim(const double*, double, Traj*);
		void runSim(std::vector<double>, double, Traj*);
		void runSim(const double*, double, double, Traj*);
		void runSim(std::vector<double>, double, double, Traj*);
		void runSim(const double*, MatrixXRd, double, double, Traj*);
		void runSim(std::vector<double>, MatrixXRd, double, double, Traj*);
		void runSim(const double *ic, MatrixXRd, std::vector<double>, Traj*);
		//@}
		
		// Utility Functions
		void clearEvents();
		void reset();
		

	private:
		/** Vector of events to consider during integration */
		std::vector<Event> events {};

		/**
		 *	Contains data recroding which events happened and at which step in the integration
		 */
		std::vector<SimEventRecord> eventOccurs {};
		
		/** a void pointer to some data object that contains data for the EOM function */
		void *eomParams = 0;

		/** Whether or not to run the simulation in reverse time */
		bool bRevTime = false;

		/** Whether or not to use variable step size when integrating */
		bool bVarStepSize = true;

		/** Whether or not to use "simple" integration; simple means no STM or other quantites are integrated
		 with the EOMs */
		bool bSimpleIntegration = false;

		/** Whether or not crash events should be created for the simulation */
		bool bMakeCrashEvents = true;

		/** Whether or not the default crash events have been created */
		bool bMadeCrashEvents = false;

		/** Absolute tolerance for integrated data, units are same as integrated data */
		double absTol = 1e-12;

		/** Relative tolerance for integrated data, units are same as integrated data */
		double relTol = 1e-14;

		/** Initial guess for time step size */
		double dtGuess = 1e-6;

		/** Number of steps to take (only applies when varStepSize = false) */
		double numSteps = 1000;

		/** Number of seconds allowed for integration; set to -1 to allow unlimitted computation time */
		int maxCompTime = -1;

		/** Timestamp at integration start */
		int startTimestamp = 0;

		/** Integrator to use for variable step-size propagations; Default is RK8PD */
		Integ_tp varStep_integ = Integ_tp::RK8PD;
		
		/** Integrator to use for fixed step-size propagations; Default is MSADAMS */
		Integ_tp fixStep_integ = Integ_tp::MSADAMS;

		void cleanEngine();
		void copyMe(const SimEngine&);
		void createCrashEvents(const SysData*);
		void free_odeiv2(gsl_odeiv2_step*, gsl_odeiv2_control*, gsl_odeiv2_evolve*, gsl_odeiv2_driver*);
		void integrate(const double*, MatrixXRd, const double*, int, Traj*);
		bool locateEvents(const double*, double, Traj*);
};


}// END of Astrohelion namespace
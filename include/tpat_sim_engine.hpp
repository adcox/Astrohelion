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
#ifndef H_SIMENGINE
#define H_SIMENGINE

#include "tpat.hpp"
 
#include "tpat_model.hpp"
#include "tpat_event.hpp"
#include "tpat_traj.hpp"
 
#include <vector>

// Forward declarations
class TPAT_Traj_BC4BP;
class TPAT_Traj_CR3BP;
class TPAT_Traj_CR3BP_LTVP;
class TPAT_Sys_Data;

/**
 *	@brief A small structure to store event occurrence records
 */
struct TPAT_Sim_EventRecord{
public:
	/**
	 *	@brief Construct an event record
	 *	@param e the event index within the simulation engine event vector
	 *	@param s the step index; which step did this event occur at?
	 */
	TPAT_Sim_EventRecord(int e, int s) : eventIx(e), stepIx(s) {}
	int eventIx;	//!< The index of the event (index from simulation engine vector of events)
	int stepIx;		//!< The index of the integration step the event occured at
};

/**
 *	@brief Performs numerical integration on any system type and produces an
 *	TPAT_Traj object
 *
 *	The simulation engine is the workhorse object for the TPAT. It
 *	holds functions to integrate equations of motion and is called by the 
 *	<tt>TPAT_Correction_Engine</tt> to compute arcs between nodes.
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
class TPAT_Sim_Engine : public TPAT{
	public:
		// Constructors
		TPAT_Sim_Engine();
		TPAT_Sim_Engine(const TPAT_Sim_Engine&);	//copy constructor
		
		//Destructor
		~TPAT_Sim_Engine();

		//Operators
		TPAT_Sim_Engine& operator =(const TPAT_Sim_Engine&);

		// Set and get functions
		void addEvent(TPAT_Event);
		double getAbsTol() const;
		std::vector<TPAT_Event> getEndEvents(TPAT_Traj*) const;
		std::vector<TPAT_Event> getEvents() const;
		std::vector<TPAT_Sim_EventRecord> getEventRecords() const;
		int getNumSteps() const;
		double getRelTol() const;
		bool makesCrashEvents() const;
		bool usesRevTime() const;
		TPAT_Verbosity_Tp getVerbosity() const;
		bool usesVarStepSize() const;
		
		void setAbsTol(double);
		void setMakeCrashEvents(bool);
		void setNumSteps(int);
		void setRelTol(double);
		void setRevTime(bool);
		void setVerbose(TPAT_Verbosity_Tp);
		void setVarStepSize(bool);

		// Simulation Methods
		void runSim(const double*, double, TPAT_Traj*);
		void runSim(std::vector<double>, double, TPAT_Traj*);
		void runSim(const double*, double, double, TPAT_Traj*);
		void runSim(std::vector<double>, double, double, TPAT_Traj*);
		
		// Utility Functions
		void clearEvents();
		void reset();
		

	private:
		/** Vector of events to consider during integration */
		std::vector<TPAT_Event> events {};

		/**
		 *	Contains data recroding which events happened and at which step in the integration
		 */
		std::vector<TPAT_Sim_EventRecord> eventOccurs {};
		
		/** a void pointer to some data object that contains data for the EOM function */
		void *eomParams = 0;

		/** Whether or not to run the simulation in reverse time */
		bool revTime = false;

		/** Describes the verbosity of this engine */
		TPAT_Verbosity_Tp verbose = TPAT_Verbosity_Tp::NO_MSG;

		/** Whether or not to use variable step size when integrating */
		bool varStepSize = true;

		/** Whether or not to use "simple" integration; simple means no STM or other quantites are integrated
		 with the EOMs */
		bool simpleIntegration = false;

		/** "Clean" means there is no data in the trajectory object. Once a simulation has been
		 run, the engine is no longer clean and will need to be cleaned before running another sim */
		bool isClean = true;

		/** Whether or not crash events should be created for the simulation */
		bool makeCrashEvents = true;

		/** Whether or not the default crash events have been created */
		bool madeCrashEvents = false;

		/** Absolute tolerance for integrated data, units are same as integrated data */
		double absTol = 1e-12;

		/** Relative tolerance for integrated data, units are same as integrated data */
		double relTol = 1e-14;

		/** Initial guess for time step size */
		double dtGuess = 1e-6;

		/** Number of steps to take (only applies when varStepSize = false) */
		double numSteps = 1000;

		void cleanEngine();
		void copyMe(const TPAT_Sim_Engine&);
		void createCrashEvents(const TPAT_Sys_Data*);
		void integrate(const double*, const double*, int, TPAT_Traj*);
		bool locateEvents(const double*, double, TPAT_Traj*);
};

#endif
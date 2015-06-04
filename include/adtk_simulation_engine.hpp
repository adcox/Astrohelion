/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (ADTK).
 *
 *  ADTK is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ADTK is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ADTK.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __H_SIMENGINE_
#define __H_SIMENGINE_

#include "adtk_event.hpp"
#include "adtk_sys_data.hpp"
#include "adtk_trajectory.hpp"
 
#include <vector>

// Forward declarations
class adtk_bcr4bpr_traj;
class adtk_cr3bp_traj;

/**
 *	@brief Performs numerical integration on any system type and produces an
 *	adtk_trajectory object
 *
 *	The simulation engine is the workhorse object for the ADTK. It
 *	holds functions to integrate equations of motion and is called by the 
 *	<tt>adtk_correction_engine</tt> to compute arcs between nodes.
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
 *	@author Andrew Cox
 *	@version June 1, 2015
 *	@copyright GNU GPL v3.0
 */
class adtk_simulation_engine{
	public:
		// Constructors
		adtk_simulation_engine();
		adtk_simulation_engine(adtk_sys_data*);
		adtk_simulation_engine(const adtk_simulation_engine&);	//copy constructor
		void createCrashEvents();
		
		//Destructor
		~adtk_simulation_engine();

		//Operators
		adtk_simulation_engine& operator =(const adtk_simulation_engine&);

		// Set and get functions
		bool usesRevTime() const;
		bool isVerbose() const;
		bool usesVarStepSize() const;
		double getAbsTol() const;
		double getRelTol() const;
		int getNumSteps() const;
		adtk_trajectory getTraj() const;
		adtk_cr3bp_traj getCR3BPTraj() const;
		adtk_bcr4bpr_traj getBCR4BPRTraj() const;

		void addEvent(adtk_event::event_t, int, bool);
		void addEvent(adtk_event);
		void setSysData(adtk_sys_data*);
		void setRevTime(bool);
		void setVerbose(bool);
		void setVarStepSize(bool);
		void setAbsTol(double);
		void setRelTol(double);
		void setNumSteps(int);

		// Simulation Methods
		void runSim(double*, double);
		void runSim(double*, double, double);

		// Utility Functions
		void reset();
		void clearEvents();

	private:
		/** Pointer to a system data object; contains characteristic quantities, among other things */
		adtk_sys_data *sysData = 0;	// set null pointers for now

		/** Pointer to a trajectory object; is set to non-null value when integration occurs */
		adtk_trajectory *traj = 0;

		/** Vector of events to consider during integration */
		std::vector<adtk_event> events;

		/** Contains data recording which events happened and at which step in the 
		 integration. Data is stored in sets of two, with the first value representing the 
		 step # along the integration and the second representing the index of the event
		 in the events vector */
		std::vector<int> eventOccurs;

		/** a void pointer to some data object that contains data for the EOM function */
		void *eomParams = 0;

		/** Whether or not to run the simulation in reverse time */
		bool revTime = false;

		/** Whether or not to print out lots of messages */
		bool verbose = false;

		/** Whether or not to use variable step size when integrating */
		bool varStepSize = true;

		/** Whether or not to use "simple" integration; simple means no STM or other quantites are integrated
		 with the EOMs */
		bool simpleIntegration = false;

		/** "Clean" means there is no data in the trajectory object. Once a simulation has been
		 run, the engine is no longer clean and will need to be cleaned before running another sim */
		bool isClean = true;

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

		void integrate(double ic[], double t[], int t_dim);
		void saveIntegratedData(double *y, double t);
		void setEOMParams();
		bool locateEvents(double*, double);
		void cleanEngine();
		void copyEngine(const adtk_simulation_engine&);
};

#endif

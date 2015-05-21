/**
 *	The simulation engine is the workhorse object for the ATDK. It
 *	holds functions to integrate equations of motion, compute
 *	zero-velocity surfaces, manifold arcs, and other useful
 *	structures
 *
 *	Author: Andrew Cox
 *
 *	Version: May 15, 2015
 */

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
 *  along with ATDK.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __H_SIMENGINE_
#define __H_SIMENGINE_

#include "adtk_sys_data.hpp"
#include "adtk_trajectory.hpp"
 
// Forward declarations
class adtk_bcr4bpr_traj;
class adtk_cr3bp_traj;

class adtk_simulation_engine{
	public:
		// Constructors
		adtk_simulation_engine();
		adtk_simulation_engine(adtk_sys_data*);
		adtk_simulation_engine(const adtk_simulation_engine&);	//copy constructor

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
		adtk_cr3bp_traj getCR3BPTraj();
		adtk_bcr4bpr_traj getBCR4BPRTraj();

		void setSysData(adtk_sys_data *d);
		void setRevTime(bool b);
		void setVerbose(bool b);
		void setVarStepSize(bool b);
		void setAbsTol(double t);
		void setRelTol(double t);
		void setNumSteps(int n);

		// Simulation Methods
		void runSim(double *ic, double tf);
		void runSim(double *ic, double t0, double tf);

		// Utility Functions
		void reset();
		
	private:
		/** Pointer to a system data object; contains characteristic quantities, among other things */
		adtk_sys_data *sysData = 0;	// set null pointers for now

		/** Pointer to a trajectory object; is set to non-null value when integration occurs */
		adtk_trajectory *traj = 0;

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

		/** Absolute tolerance for integrated data, units are same as integrated data */
		double absTol = 1e-12;

		/** Relative tolerance for integrated data, units are same as integrated data */
		double relTol = 1e-14;

		/** Initial guess for time step size */
		double dtGuess = 1e-6;

		/** Number of steps to take (only applies when varStepSize = false) */
		double numSteps = 1000;

		void integrate(double ic[], double t[], int t_dim);
		void saveIntegratedData(double *y, double t, bool);
		void setEOMParams();
		void cleanEngine();
		void printMessage(const char*,...);
};

#endif

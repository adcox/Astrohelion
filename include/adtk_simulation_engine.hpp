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

		//Destructor
		~adtk_simulation_engine();

		// Set and get functions
		bool usesRevTime();
		bool isVerbose();
		double getAbsTol();
		double getRelTol();
		adtk_cr3bp_traj getCR3BPTraj();
		adtk_bcr4bpr_traj getBCR4BPRTraj();

		void setSysData(adtk_sys_data *d);
		void setRevTime(bool b);
		void setVerbose(bool b);
		void setAbsTol(double t);
		void setRelTol(double t);

		// Simulation Methods
		void runSim(double *ic, double tf);
		void runSim(double *ic, double t0, double tf);

		// Utility Functions

	private:
		adtk_sys_data *sysData;
		adtk_trajectory *traj;

		bool revTime;
		bool verbose;
		bool simpleIntegration;
		double absTol, relTol, dtGuess;

		// void cr3bp_integrate(double ic[], double t[], double mu, int t_dim);
		void integrate(double ic[], double t[], int t_dim);
		void saveIntegratedData(double *y, double t, bool);
};

#endif

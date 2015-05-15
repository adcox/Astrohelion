/**
 *	Simulation Engine Header
 */

#ifndef __H_SIMENGINE_
#define __H_SIMENGINE_

#include "adtk_sys_data.hpp"
#include "adtk_trajectory.hpp"
#include "adtk_cr3bp_traj.hpp"

#include <vector>

class adtk_simulation_engine{
	public:
		// Constructors
		adtk_simulation_engine();

		//Destructor
		~adtk_simulation_engine();

		// Set and get functions
		adtk_sys_data* getSysData();
		bool usesRevTime();
		bool isVerbose();
		double getAbsTol();
		double getRelTol();
		adtk_cr3bp_traj getCR3BPTraj();

		void setSysData(adtk_sys_data *d);
		void setRevTime(bool b);
		void setVerbose(bool b);
		void setAbsTol(double t);
		void setRelTol(double t);

		// Simulation Methods
		void runSim(double *ic, double tf);
		void runSim(double *ic, double t0, double tf);

	private:
		adtk_sys_data *sysData;
		adtk_trajectory *traj;

		bool revTime;
		bool verbose;
		double absTol, relTol, dtGuess;

		void cr3bp_integrate(double ic[], double t[], double mu, int t_dim);
		void saveIntegratedData(double *y, double t);
};

#endif

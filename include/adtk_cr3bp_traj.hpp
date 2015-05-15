/**
 *	Header for CR3BP Trajectory data object
 */

#ifndef __H_CR3BP_TRAJ_
#define __H_CR3BP_TRAJ_

#include "adtk_trajectory.hpp"
#include "adtk_cr3bp_sys_data.hpp"

class adtk_cr3bp_traj : public adtk_trajectory{
	private:
		/** Vector to hold jacobi constants along the path */
		std::vector<double> jacobi;

		/** A system data object specific to the CR3BP */
		adtk_cr3bp_sys_data sysData;

	public:
		adtk_cr3bp_traj();
		adtk_cr3bp_traj(int);
		
		adtk_cr3bp_traj& operator= (const adtk_cr3bp_traj&);

		double getJC(int);
		std::vector<double>* getJC();
		
		void setJC(std::vector<double>);

		adtk_cr3bp_sys_data getSysData();
		void setSysData(adtk_cr3bp_sys_data);
};

#endif
//END
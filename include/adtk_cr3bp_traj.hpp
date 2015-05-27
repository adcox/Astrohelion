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
#ifndef __H_CR3BP_TRAJ_
#define __H_CR3BP_TRAJ_

#include "adtk_cr3bp_sys_data.hpp"
#include "adtk_trajectory.hpp"

/**
 *	A derivative class of the adtk_trajectory super-class. This object
 *	contains trajectory information specific to the CR3BP
 *
 *	Author: Andrew Cox
 *
 *	Version: May 15, 2015
 */
class adtk_cr3bp_traj : public adtk_trajectory{
	public:
		// *structors
		adtk_cr3bp_traj();
		adtk_cr3bp_traj(int);
		adtk_cr3bp_traj(adtk_cr3bp_sys_data);
		adtk_cr3bp_traj(const adtk_cr3bp_traj&);

		// Operators
		adtk_cr3bp_traj& operator= (const adtk_cr3bp_traj&);
		friend adtk_cr3bp_traj operator +(const adtk_cr3bp_traj &lhs, const adtk_cr3bp_traj &rhs);

		// Set and Get Functions
		double getJC(int);
		std::vector<double>* getJC();
		void setJC(std::vector<double>);
		adtk_cr3bp_sys_data getSysData();
		void setSysData(adtk_cr3bp_sys_data);
		adtk_sys_data::system_t getType() const;

		// Utility functions
		void saveToMat(const char*);
		void setLength();
	private:
		/** Vector to hold jacobi constants along the path */
		std::vector<double> jacobi;

		/** A system data object specific to the CR3BP */
		adtk_cr3bp_sys_data sysData;
		
		void saveJacobi(mat_t*);
		void saveSysData(mat_t*);
};

#endif
//END
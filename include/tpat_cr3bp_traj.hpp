/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (TPAT).
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
#ifndef H_CR3BP_TRAJ
#define H_CR3BP_TRAJ

#include "tpat_cr3bp_sys_data.hpp"
#include "tpat_trajectory.hpp"

#include "matio.h"

// forward declarations
class tpat_cr3bp_nodeset;

/**
 *	@brief A derivative class of the tpat_trajectory object, which
 *	contains trajectory information specific to the CR3BP
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 * 	@copyright GNU GPL v3.0
 */
class tpat_cr3bp_traj : public tpat_trajectory{
	public:
		// *structors
		tpat_cr3bp_traj();
		tpat_cr3bp_traj(int);
		tpat_cr3bp_traj(tpat_cr3bp_sys_data);
		tpat_cr3bp_traj(const tpat_cr3bp_traj&);
		static tpat_cr3bp_traj fromNodeset(tpat_cr3bp_nodeset); 

		// Operators
		tpat_cr3bp_traj& operator= (const tpat_cr3bp_traj&);
		friend tpat_cr3bp_traj operator +(const tpat_cr3bp_traj &lhs, const tpat_cr3bp_traj &rhs);
		tpat_cr3bp_traj& operator +=(const tpat_cr3bp_traj&);

		// Set and Get Functions
		double getJC(int) const;
		std::vector<double>* getJC();
		tpat_cr3bp_sys_data getSysData() const;
		tpat_sys_data::system_t getType() const;

		void setJC(std::vector<double>);
		void setSysData(tpat_cr3bp_sys_data);
		// Utility functions
		void saveToMat(const char*);
		void setLength();
	private:
		/** Vector to hold jacobi constants along the path */
		std::vector<double> jacobi;

		/** A system data object specific to the CR3BP */
		tpat_cr3bp_sys_data sysData;
		
		void saveJacobi(mat_t*);
		void copyMe(const tpat_cr3bp_traj&);
};

#endif
//END
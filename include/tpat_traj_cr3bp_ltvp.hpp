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
#ifndef H_CR3BP_LTVP_TRAJ
#define H_CR3BP_LTVP_TRAJ

#include "tpat_sys_data_cr3bp_ltvp.hpp"
#include "tpat_traj.hpp"

#include "matio.h"

// forward declarations
class tpat_nodeset_cr3bp;

/**
 *	@brief A derivative class of the tpat_traj object, which
 *	contains trajectory information specific to the Low Thrust, 
 *	Velocity-Pointing CR3BP
 *
 *	@author Andrew Cox
 *	@version July 20, 2015
 * 	@copyright GNU GPL v3.0
 */
class tpat_traj_cr3bp_ltvp : public tpat_traj{
	public:
		// *structors
		tpat_traj_cr3bp_ltvp();
		tpat_traj_cr3bp_ltvp(int);
		tpat_traj_cr3bp_ltvp(tpat_sys_data_cr3bp_ltvp);
		tpat_traj_cr3bp_ltvp(const tpat_traj_cr3bp_ltvp&);

		// Operators
		tpat_traj_cr3bp_ltvp& operator= (const tpat_traj_cr3bp_ltvp&);

		// Set and Get Functions
		std::vector<double>* getJacobi();
		double getJacobi(int) const;
		tpat_sys_data_cr3bp_ltvp getSysData() const;
		tpat_sys_data::system_t getType() const;

		void setSysData(tpat_sys_data_cr3bp_ltvp);
		
		// Utility functions
		void saveToMat(const char*);

	private:
		/** A system data object specific to the Low Thrust, Velocity-Pointing CR3BP */
		tpat_sys_data_cr3bp_ltvp sysData;
		
		void saveJacobi(mat_t*);
		void copyMe(const tpat_traj_cr3bp_ltvp&);
};

#endif
//END
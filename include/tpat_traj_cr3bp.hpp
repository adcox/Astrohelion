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
#ifndef H_CR3BP_TRAJ
#define H_CR3BP_TRAJ

#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_traj.hpp"

#include "matio.h"

// forward declarations
class tpat_nodeset_cr3bp;

/**
 *	@brief A derivative class of the tpat_traj object, which
 *	contains trajectory information specific to the CR3BP
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 * 	@copyright GNU GPL v3.0
 */
class tpat_traj_cr3bp : public tpat_traj{
	public:
		// *structors
		tpat_traj_cr3bp();
		tpat_traj_cr3bp(int);
		tpat_traj_cr3bp(tpat_sys_data_cr3bp);
		tpat_traj_cr3bp(const tpat_traj_cr3bp&);
		static tpat_traj_cr3bp fromNodeset(tpat_nodeset_cr3bp); 

		// Operators
		tpat_traj_cr3bp& operator= (const tpat_traj_cr3bp&);
		friend tpat_traj_cr3bp operator +(const tpat_traj_cr3bp &lhs, const tpat_traj_cr3bp &rhs);

		// Set and Get Functions
		double getJacobi(int) const;
		std::vector<double>* getJacobi();
		tpat_sys_data_cr3bp getSysData() const;
		tpat_sys_data::system_t getType() const;

		void setJacobi(std::vector<double>);
		void setSysData(tpat_sys_data_cr3bp);
		
		// Utility functions
		void saveToMat(const char*);
	private:
		/** A system data object specific to the CR3BP */
		tpat_sys_data_cr3bp sysData;
		
		void copyMe(const tpat_traj_cr3bp&);
};

#endif
//END
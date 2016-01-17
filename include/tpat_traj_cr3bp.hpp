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

#include "tpat_traj.hpp"

// Forward Declarations
class tpat_nodeset_cr3bp;
class tpat_sys_data_cr3bp;

/**
 *	@brief A derivative class of the tpat_traj object, contains 
 *	trajectory information specific to the CR3BP
 *
 *	@author Andrew Cox
 *	@version August 30, 2015
 *	@copyright GNU GPL v3.0
 */
class tpat_traj_cr3bp : public tpat_traj{

public:
	// *structors
	tpat_traj_cr3bp(const tpat_sys_data_cr3bp*);
	tpat_traj_cr3bp(const tpat_traj_cr3bp&);
	tpat_traj_cr3bp(const tpat_arc_data&);
	
	static tpat_traj_cr3bp fromNodeset(tpat_nodeset_cr3bp);
	
	// Set and Get Functions
	double getJacobi(int) const;

	void setJacobi(int, double);
	void saveToMat(const char*) const;
private:
	void initExtraParam();
	
};

#endif
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
class TPAT_Nodeset_CR3BP;
class TPAT_Sys_Data_CR3BP;

/**
 *	@brief A derivative class of the TPAT_Traj object, contains 
 *	trajectory information specific to the CR3BP
 *
 *	@author Andrew Cox
 *	@version August 30, 2015
 *	@copyright GNU GPL v3.0
 */
class TPAT_Traj_CR3BP : public TPAT_Traj{

public:
	// *structors
	TPAT_Traj_CR3BP(const TPAT_Sys_Data_CR3BP*);
	TPAT_Traj_CR3BP(const TPAT_Traj_CR3BP&);
	TPAT_Traj_CR3BP(const TPAT_Base_Arcset&);
	baseArcsetPtr create(const TPAT_Sys_Data*) const;
	baseArcsetPtr clone() const;

	// Operators
	TPAT_Traj& operator +=(const TPAT_Traj&);
	
	// Set and Get Functions
	double getJacobiByIx(int) const;
	void setJacobiByIx(int, double);

	// Utility
	void readFromMat(const char*);
	void saveToMat(const char*) const;
private:

	void initExtraParam();
};

#endif
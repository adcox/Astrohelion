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

#include "tpat_traj.hpp"

// Forward Declarations
class TPAT_Sys_Data_CR3BP_LTVP;

/**
 *	@brief A derivative class of the TPAT_Traj object, which
 *	contains trajectory information specific to the Low Thrust, 
 *	Velocity-Pointing CR3BP
 *
 *	@author Andrew Cox
 *	@version September 2, 2015
 *	@copyright GNU GPL v3.0
 */
class TPAT_Traj_CR3BP_LTVP : public TPAT_Traj{

public:
	// *structors
	TPAT_Traj_CR3BP_LTVP(const TPAT_Sys_Data_CR3BP_LTVP*);
	TPAT_Traj_CR3BP_LTVP(const TPAT_Traj_CR3BP_LTVP&);
	TPAT_Traj_CR3BP_LTVP(const TPAT_Base_Arcset&);
	baseArcsetPtr create(const TPAT_Sys_Data*) const;
	baseArcsetPtr clone() const;

	// Set and Get Functions
	double getJacobiByIx(int) const;
	double getMassByIx(int) const;

	void setJacobiByIx(int, double);
	void setMassByIx(int, double);
	
	void saveToMat(const char*) const;
private:
	void initExtraParam();
};

#endif
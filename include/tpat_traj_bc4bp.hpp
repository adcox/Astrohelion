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

#ifndef H_BCR4BPR_TRAJ
#define H_BCR4BPR_TRAJ

#include "tpat_traj.hpp"

// Forward Declarations
class TPAT_Sys_Data_BC4BP;
class TPAT_Nodeset_BC4BP;

/**
 *	@brief A derivative class of the TPAT_Traj object that
 *	contains trajectory information specific to the CR3BP
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 *	@copyright GNU GPL v3.0
 */
class TPAT_Traj_BC4BP : public TPAT_Traj{

public:
	// *structors
	TPAT_Traj_BC4BP(const TPAT_Sys_Data_BC4BP*);
	TPAT_Traj_BC4BP(const TPAT_Traj_BC4BP&);
	TPAT_Traj_BC4BP(const TPAT_Base_Arcset&);
	TPAT_Traj_BC4BP(const char*);
	baseArcsetPtr create(const TPAT_Sys_Data*) const;
	baseArcsetPtr clone() const;
	//static TPAT_Traj_BC4BP fromNodeset(TPAT_Nodeset_BC4BP);

	// Set and Get Functions
	double getTheta0();
	double getPhi0();
	double getGamma();
	std::vector<double> get_dqdTByIx(int);

	void set_dqdTByIx(int, const double*);
	void set_dqdTByIx(int, std::vector<double>);

	void readFromMat(const char*);
	void saveToMat(const char*) const;
private:

	void initExtraParam();
};

#endif
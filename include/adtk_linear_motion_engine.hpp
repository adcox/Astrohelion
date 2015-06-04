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
 *  along with ADTK.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __H_LIN_MOTION_
#define __H_LIN_MOTION_

// Forward delcarations
class adtk_cr3bp_traj;

/**
 *	@brief An engine that will generate a trajectory from the linearized EOMs
 */
class adtk_linear_motion_engine{
	public:
		enum motion_t {NONE, HYP, ELLIP, SPO, LPO, MPO, CONVERGE, DIVERGE};

		// *structors
		adtk_linear_motion_engine();

		// Operators

		// Set and Get
		void setTol(double);
		const char* getTypeStr(motion_t) const;

		// Misc
		adtk_cr3bp_traj getCR3BPLinear(int, double[3], const char*, const char*);
		adtk_cr3bp_traj getCR3BPLinear(int, double[3], motion_t, const char*, const char*);
	private:
		double t_step = 0.001;

		double rots = 1;

		double tol = 1e-14;

		double nu = 1;
};

#endif
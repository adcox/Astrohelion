/**
 *	A derivative class of the adtk_trajectory super-class. This object
 *	contains trajectory information specific to the CR3BP
 *
 *	Author: Andrew Cox
 *
 *	Version: May 15, 2015
 */

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
#ifndef __H_BCR4BPR_TRAJ_
#define __H_BCR4BPR_TRAJ_

#include "adtk_bcr4bpr_sys_data.hpp"
#include "adtk_trajectory.hpp"

class adtk_bcr4bpr_traj : public adtk_trajectory{
	private:

		/** Angle (radians) between the P1/P2 line and the inertial x-axis at time t = 0 */
		double theta0;

		/** Angle (radians) between the P2/P3 line (projected onto the inertial XY plane) and
		 the inertial x-axis at time t = 0 */
		double phi0;

		/** Inclination (radians) of P2/P3 orbital plane relative to the P1/P2 orbital plane;
		 Gamma is constant during the integration */
		double gamma;

		/** A system data object specific to the BCR4BPR */
		adtk_bcr4bpr_sys_data sysData;

	public:
		adtk_bcr4bpr_traj();
		adtk_bcr4bpr_traj(int);

		// Set and Get Functions
		double getTheta0();
		double getPhi0();
		double getGamma();
		adtk_bcr4bpr_sys_data getSysData();

		void setTheta0(double t);
		void setPhi0(double p);
		void setGamma(double g);
		void setSysData(adtk_bcr4bpr_sys_data);
};

#endif
//END
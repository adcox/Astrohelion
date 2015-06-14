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
#ifndef __H_BCR4BPR_SYS_DATA_
#define __H_BCR4BPR_SYS_DATA_

#include "tpat_sys_data.hpp"
#include "matio.h"
#include "tpat_constants.hpp"

/**
 *	@brief This derivative class of the tpat_sys_data super-class
 *	contains information specific to the BCR4BPR
 *
 *	@author Andrew Cox
 *	@version May 18, 2015
 *	@copyright GNU GPL v3.0
 */
class tpat_bcr4bpr_sys_data : public tpat_sys_data{
	public:
		tpat_bcr4bpr_sys_data();
		tpat_bcr4bpr_sys_data(std::string, std::string, std::string);
		tpat_bcr4bpr_sys_data(const tpat_bcr4bpr_sys_data&);

		tpat_bcr4bpr_sys_data& operator=(const tpat_bcr4bpr_sys_data&);
		
		double getMu() const;
		double getNu() const;
		double getK() const;
		double getCharLRatio() const;

		double getTheta0() const;
		double getPhi0() const;
		double getGamma() const;

		void setTheta0(double t);
		void setPhi0(double p);
		void setGamma(double g);

		void saveToMat(mat_t*);
		
		/** Time when geometry is at reference orientation (theta = phi = 0), seconds, J2000, UTC */
		static double REF_EPOCH;	// 2005/06/21 18:21:35

	private:
		double mu = 0;			//!< Mass ratio between P2 + P3 and total
		double nu = 0;			//!< Mass ratio between P3 and total
		double k = 1/100;		//!< Scaling constant, non-dim 

		/** Ratio between P3's orbital radius and P2's orbital radius */
		double charLRatio = 0;

		/** Angle between the P1/P2 line and the inertial x-axis, radians */
		double theta0 = 0;

		/** Angle between the P2/P3 line (projected into inertial XY plane) and the
		x-axis, radians */
		double phi0 = 0;

		/** Inclination of the P2/P3 orbital plane relative to the P1/P2 orbital plane;
		This angle does not change over time; radians */
		double gamma = 5.14*PI/180;

		void copyData(const tpat_bcr4bpr_sys_data&);
};

#endif
//END
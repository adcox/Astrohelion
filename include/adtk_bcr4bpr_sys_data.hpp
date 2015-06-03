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
#ifndef __H_BCR4BPR_SYS_DATA_
#define __H_BCR4BPR_SYS_DATA_

#include "adtk_sys_data.hpp"

#include "adtk_constants.hpp"

/**
 *	@brief This derivative class of the adtk_sys_data super-class
 *	contains information specific to the BCR4BPR
 *
 *	@author Andrew Cox
 *	@version May 18, 2015
 *	@copyright GNU GPL v3.0
 */
class adtk_bcr4bpr_sys_data : public adtk_sys_data{
	public:
		adtk_bcr4bpr_sys_data();
		adtk_bcr4bpr_sys_data(std::string, std::string, std::string);

		adtk_bcr4bpr_sys_data& operator=(const adtk_bcr4bpr_sys_data&);
		
		double getMu() const;
		double getNu() const;
		double getK() const;
		double getCharLRatio() const;
		std::string getPrimary(int n) const;	//We override this function, so re-declare it

		double getTheta0() const;
		double getPhi0() const;
		double getGamma() const;

		void setTheta0(double t);
		void setPhi0(double p);
		void setGamma(double g);
	private:
		/** Mass ratio between P2 + P3 and total*/
		double mu = 0;

		/** Mass ratio between P3 and total */
		double nu = 0;

		/** Scaling constant, non-dim */
		double k = 1/100;

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

		/** Name of P1 */
		std::string P1 = "P1";

		/** Name of P2 */
		std::string P2 = "P2";

		/** Name of P3 */
		std::string P3 = "P3";
};

#endif
//END
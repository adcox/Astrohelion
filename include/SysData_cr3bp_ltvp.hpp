/**
 *  @file SysData_cr3bp_ltvp.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of Astrohelion
 *
 *  Astrohelion is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Astrohelion is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Astrohelion.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "SysData_cr3bp.hpp"

#include "DynamicsModel_cr3bp_ltvp.hpp"

namespace astrohelion{

/**
 *	@brief System data object for the CR3BP Low Thrust, Velocity-Pointing system
 */
class SysData_cr3bp_ltvp : public SysData_cr3bp{

	public:
		SysData_cr3bp_ltvp();
		SysData_cr3bp_ltvp(std::string, std::string, double, double, double);
		SysData_cr3bp_ltvp(const SysData_cr3bp_ltvp&);
		SysData_cr3bp_ltvp(const char*);
		
		SysData_cr3bp_ltvp& operator=(const SysData_cr3bp_ltvp&);
		
		const DynamicsModel* getDynamicsModel() const;
		
		double getIsp() const;
		double getThrust() const;
		double getM0() const;

		void setIsp(double);
		void setIspDim(double);
		void setM0(double);
		void setM0Dim(double);
		void setThrust(double);
		void setThrustDim(double);

		void saveToMat(const char*) const;
		void saveToMat(mat_t*) const;
	private:
		/** The dynamic model that governs motion for this system*/
		DynamicsModel_cr3bp_ltvp model = DynamicsModel_cr3bp_ltvp();
		void readFromMat(mat_t*);
};


}// END of Astrohelion namespace
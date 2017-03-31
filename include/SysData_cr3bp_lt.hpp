/**
 *  \file SysData_cr3bp_lt.hpp
 *	\brief 
 *	
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "ControlLaw_cr3bp_lt.hpp"
#include "DynamicsModel_cr3bp_lt.hpp"

namespace astrohelion{

/**
 *	\ingroup model
 *	\brief System data object for the CR3BP Low Thrust, Velocity-Pointing system
 */
class SysData_cr3bp_lt : public SysData_cr3bp{

	public:
		/**
		 *  \name *structors
		 *  \{
		 */
		SysData_cr3bp_lt();
		SysData_cr3bp_lt(std::string, std::string, double, double, double);
		SysData_cr3bp_lt(const SysData_cr3bp_lt&);
		SysData_cr3bp_lt(const char*);
		//\}

		SysData_cr3bp_lt& operator=(const SysData_cr3bp_lt&);
		
		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		const ControlLaw* getControlLaw() const;
		const DynamicsModel* getDynamicsModel() const;
		
		double getIsp() const;
		double getThrust() const;
		double getThrust_dim() const;
		double getMass() const;

		void setIsp(double);
		void setMass(double);
		void setThrust(double);
		void setThrust_dim(double);
		//\}
		
		void saveToMat(const char*) const override;
		void saveToMat(mat_t*) const;
	private:
		/** The dynamic model that governs motion for this system*/
		DynamicsModel_cr3bp_lt model = DynamicsModel_cr3bp_lt();

		/** Control law object for this system */
		ControlLaw_cr3bp_lt control = ControlLaw_cr3bp_lt();
		
		void readFromMat(mat_t*);
};


}// END of Astrohelion namespace
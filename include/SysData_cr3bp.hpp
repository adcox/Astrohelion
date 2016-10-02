/**
 *  @file SysData_cr3bp.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2016, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "SysData.hpp"

#include "DynamicsModel_cr3bp.hpp"

#include "matio.h"
#include <string>

namespace astrohelion{

/**
 *	@ingroup model cr3bp
 *	@brief A derivative class of the SysData object which
 *	contains information specific to the CR3BP
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 *	@copyright GNU GPL v3.0
 */
class SysData_cr3bp : public SysData{
	public:
		/**
		 *  @name *structors
		 *  @{
		 */
		SysData_cr3bp();
		SysData_cr3bp(std::string, std::string);
		SysData_cr3bp(const SysData_cr3bp&);
		SysData_cr3bp(const char*);
		virtual ~SysData_cr3bp();
		//@}

		SysData_cr3bp& operator=(const SysData_cr3bp&);
		
		/**
		 *  @name Set and Get Functions
		 *  @{
		 */
		const DynamicsModel* getDynamicsModel() const;
		double getMu() const;
		//@}
		
		void saveToMat(const char*) const override;
		void saveToMat(mat_t*) const;
		
	protected:
		void initFromPrimNames(std::string, std::string);
		void readFromMat(mat_t*);
		
	private:
		/** The dynamic model that governs motion for this system*/
		DynamicsModel_cr3bp model = DynamicsModel_cr3bp();
};


}// END of Astrohelion namespace
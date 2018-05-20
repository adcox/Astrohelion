/**
 *  @file SysData_2bp.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version August 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "matio.h"
#include <string>

#include "SysData.hpp"

#include "DynamicsModel_2bp.hpp"

namespace astrohelion{

/**
 *	\ingroup model 2bp
 *	@brief A derivative class of the SysData object which
 *	contains information specific to the 2BP
 *
 *	@author Andrew Cox
 *	@version August 31, 2015
 *	@copyright GNU GPL v3.0
 */
class SysData_2bp : public SysData{
	public:
		/**
		 *  \name *structors
		 *  \{
		 */
		SysData_2bp();
		SysData_2bp(std::string);
		SysData_2bp(const SysData_2bp&);
		void initFromFile(const char*);
		virtual ~SysData_2bp();
		//\}
		
		/**
		 *  \name Operators
		 *  \{
		 */
		SysData_2bp& operator=(const SysData_2bp&);
		//\}

		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		const DynamicsModel* getDynamicsModel() const;

		double getMu() const;
		//\}

		/**
		 *  \name File I/O
		 *  \{
		 */
		void saveToMat(const char*) const override;
		void saveToMat(mat_t*) const;
		//\}

	protected:
		/**
		 *  \name *structors
		 *  \{
		 */
		void initFromPrimNames(std::string);
		//\}

		/**
		 *  \name File I/O
		 *  \{
		 */
		void readFromMat(mat_t*) override;
		//\}
		
	private:
		/** The dynamic model that governs motion for this system*/
		DynamicsModel_2bp model = DynamicsModel_2bp();
};


}// END of Astrohelion namespace
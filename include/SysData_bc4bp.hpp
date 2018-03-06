/**
 *  @file SysData_bc4bp.hpp
 *	@brief 
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
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

#include "SysData.hpp"

#include "Common.hpp" 
#include "DynamicsModel_bc4bp.hpp"



namespace astrohelion{

/**
 *	\ingroup model bc4bp
 *	@brief This derivative class of the SysData super-class
 *	contains information specific to the BCR4BPR
 *
 *	@author Andrew Cox
 *	@version May 18, 2015
 *	@copyright GNU GPL v3.0
 */
class SysData_bc4bp : public SysData{
	public:
		/**
		 *  \name *structors
		 *  \{
		 */
		SysData_bc4bp();
		SysData_bc4bp(std::string, std::string, std::string);
		SysData_bc4bp(const SysData_bc4bp&);
		SysData_bc4bp(const char*);
		//\}

		/**
		 *  \name Operators
		 *  \{
		 */
		SysData_bc4bp& operator=(const SysData_bc4bp&);
		//\}

		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		const DynamicsModel* getDynamicsModel() const;
		double getMu() const;
		double getNu() const;
		double getK() const;
		double getCharLRatio() const;

		double getEpoch0() const;
		double getTheta0() const;
		double getPhi0() const;
		double getGamma() const;

		void setEpoch0(double T);
		void setTheta0(double t);
		void setPhi0(double p);
		void setGamma(double g);
		//\}
		
		/**
		 *  \name File I/O
		 *  \{
		 */
		void saveToMat(const char*) const override;
		void saveToMat(mat_t*) const;
		//\}

		/** Time when geometry is at reference orientation (theta = phi = 0), seconds, J2000, UTC */
		static double REF_EPOCH;	// 2005/06/21 18:21:35

	private:
		/** The dynamic model that governs motion for this system*/
		DynamicsModel_bc4bp model = DynamicsModel_bc4bp();

		/**
		 *  \name *structors
		 *  \{
		 */
		void initFromPrimNames(std::string, std::string, std::string);
		//\}

		/**
		 *  \name File I/O
		 *  \{
		 */
		void readFromMat(mat_t*);
		//\}
};

}// END of Astrohelion namespace
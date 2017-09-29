/**
 *  \file LinMotionEngine_cr3bp.hpp
 *	\brief 
 *	
 *	\author Andrew Cox
 *	\version September 28, 2017
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

#include "LinMotionEngine.hpp"

namespace astrohelion{

// Forward Declarations
class Arcset_cr3bp;

/**
 *  \brief Describe the type of linear motion to construct
 *  \details
 *  
 *  \todo Add description of the different linearization options
 */
class LinMotion_cr3bp_tp : public LinMotion_tp{
public:
	static const unsigned int SPO = 21;		//!< Short-Period Orbit, applies only to Case I linearizations near triangular points
	static const unsigned int LPO = 22;		//!< Long-Period Orbit, applies only to Case I linearizations near triangular points
	static const unsigned int MPO = 23;		//!< Mixed-Period Orbit, applies only to Case I linearizations near triangular points
};

/**
 *  \brief Construct trajectories in the linearized dynamics of the CR3BP
 */
class LinMotionEngine_cr3bp : public LinMotionEngine{
	public:
		LinMotionEngine_cr3bp();

		/**
		 *  \name Set and Get Functions
		 *  \{
		 */
		double getMPORatio() const;

		void setMPORatio(double);
		const char* getTypeStr(unsigned int) const;
		//\}

		/**
		 *  \name Orbit Generation
		 *  \{
		 */
		void getLinear(int, double[3], unsigned int, Arcset_cr3bp*);
		void getLinear(int, double[3], double, double, unsigned int, Arcset_cr3bp*);
		void getLiss(int, double, bool, double, double, double, Arcset_cr3bp*);
		//\}
	
		void reset();
	private:
		/** Ratio between SPO and LPO behavior when constructing an MPO */
		double nu = 1;
};

}
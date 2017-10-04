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

// Forward declarations
class Arcset_cr3bp_lt;
class ControlLaw_cr3bp_lt;

class LinMotionEngine_cr3bp_lt : public LinMotionEngine{
	public:
		LinMotionEngine_cr3bp_lt();

		/**
		 *  \name Orbit Generation
		 *  \{
		 */
		void getLinear(double[3], double, double, double[3], unsigned int, Arcset_cr3bp_lt*, ControlLaw_cr3bp_lt*, unsigned int numNodes = 2);
		//\}
};

}// END of astrohelion namespace
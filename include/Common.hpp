/**
 *	@file Common.hpp
 *	@brief Contains values for physical constants like pi, G, etc. and custom data wrappers
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 *	@copyright GNU GPL v3.0
 */

/*
 *	Astrohelion 
 *	Copyright 2016, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrohelion.
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

#include <complex>

namespace astrohelion{

/**
 * @brief Provides a scale larger than true/false to set verbosity for Core routines
 */
enum class Verbosity_tp{
	NO_MSG 		= 0,	//!< Absolutely no messages
	SOME_MSG	= 1,	//!< Some useful messages are printed to standard output
	ALL_MSG		= 2,	//!< All useful and not-so-useful messages are printed to standard output
	DEBUG 		= 3		//!< Even more messages; helpful for debugging
};

// Custom types
typedef std::complex<double> cdouble;	//!< A complex double
typedef std::complex<int> cint;			//!< A complex integer

/** Universal Gravity Constant, km^3/kg-s^2 */
const double G = 6.67384e-20;

/** 1 G acceleration, km/s^2 */
const double G_GRAV_0 = 9.8065/1000;

/** Pi */
const double PI = 3.14159265358979323846264338327950;

}// END of Astrohelion namespace
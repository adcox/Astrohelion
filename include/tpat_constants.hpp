/**
 *	@file tpat_constants.hpp
 *	@brief Contains values for physical constants like pi, G, etc. and custom data wrappers
 *
 *
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 *	@copyright GNU GPL v3.0
 */

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
#ifndef __H_CONST__
#define __H_CONST__

#include <complex>
 
// Custom types
typedef std::complex<double> cdouble;	//!< A complex double
typedef std::complex<int> cint;			//!< A complex integer

/** Universal Gravity Constant, km^3/kg-s^2 */
const double G = 6.67384e-20;

/** Pi */
const double PI = 3.14159265358979323846264338327950;

#endif
//END
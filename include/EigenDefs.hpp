/**
 *  @file EigenDefs.hpp
 *	@brief Definitions for matrices, vectors, and other objects from the Eigen library.
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
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

#include "Eigen/Core"

namespace astrohelion{
/**
 * A dynamically sized matrix of doubles, with data stored in row-major order
 */
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXRd;

/**
 * A dynamically sized matrix of complex doubles with data stored in row-major order
 */
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXRcd;

/**
 * A 3x3 matrix of doubles with data stored in row-major order
 */
typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> Matrix3Rd;

}// END of Astrohelion namespace
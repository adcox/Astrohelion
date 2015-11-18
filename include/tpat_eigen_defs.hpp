/*
 *	Trajectory Propagation and Analysis Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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
 
 /**
  *  This header file contains definitions for matrices, vectors, and other objects
  *  from the Eigen library. This header does not directly include the Eigen library
  *  to avoid compilation headaches. Any implementing classes must include Eigen libraries
  *  to gain access to Eigen functionality.
  */
#ifndef H_TPAT_EIGEN_DEFS
#define H_TPAT_EIGEN_DEFS

#include "Eigen/Core"

/**
 * A dynamically sized matrix of doubles, with data stored in row-major order
 */
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXRd;

typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXRcd;

typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> Matrix3Rd;

#endif
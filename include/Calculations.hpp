/**
 *	\file Calculations.hpp
 *	\brief Includes miscellaneous calculation functions
 *
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
 */
 
/*
 *	Astrohelion 
 *	Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
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

#include "Common.hpp"
#include "EigenDefs.hpp"
 
#include <vector>

namespace astrohelion{

// Forward declarations
class BaseArcset;
class MultShootEngine;
class MultShootData;
class Node;
class Nodeset;
class Nodeset_bc4bp;
class Nodeset_cr3bp;
class SysData_2bp;
class SysData_bc4bp;
class SysData_cr3bp;
class Traj;
class Traj_2bp;
class Traj_bc4bp;
class Traj_cr3bp;

/**
 *	\brief Describes the plane a periodic orbit can be mirrored across
 */
enum class Mirror_tp{
	MIRROR_XZ,		//!< Mirror over the XZ-Plane; x, z, and y-dot can be fixed if desired
	MIRROR_XY,		//!< Mirror over the XY-Plane; x, y, and z-dot can be fixed if desired
	MIRROR_YZ,		//!< Mirror over the YZ-Plane; y, z, and x-dot can be fixed if desired
	MIRROR_X_AX_H, 	//!< Mirror over the X-Axis, orbit is mostly in xy-plane; x, y-dot, and z-dot can be fixed if desired
	MIRROR_X_AX_V	//!< Mirror over the X-Axis, orbit is mostly in xz-plane; x, y-dot, and z-dot can be fixed if desired
};

/**
 *	\brief Eigenvalue pair types
 *
 *	Eigenvalues come in three different types of pairs. They will
 *	either be complex, real, or exactly equal to 1.0. Since the 
 *	matrices are real, all complex eigenvalues will come in conjucate
 *	pairs. In this type of problem (trajectory design), the other STM
 *	eigenvalues also come in pairs: 2 1.0 eigenvalues, or two real 
 *	eigenvalues that are reciprocals.
 */
enum class EigValSet_tp{
	EIGSET_COMP_CONJ,	//!< Complex conjugate pair
	EIGSET_ONES,		//!< Exactly equal to 1.0
	EIGSET_REAL_RECIP	//!< Real, reciprocal pair
};

/**
 *  \name General Utility Functions
 *  \{
 */
double dateToEphemerisTime(const char*);
double dateToGST(double, double, double, double, double, double);
double gregorianToJulian(double, double, double, double, double, double);
std::vector<double> familyCont_LS(int, double, std::vector<int>, std::vector<double>);
MatrixXRd getMirrorMat(Mirror_tp);
double getStabilityIndex(std::vector<cdouble>);
MatrixXRd solveAX_eq_B(MatrixXRd, MatrixXRd);
std::vector<unsigned int> sortEig(std::vector<cdouble>, std::vector<MatrixXRcd>);
Node interpPointAtTime(const Traj*, double);
/** \} */

/**
 *  \name Orbit Determination
 *  \{
 */
std::vector<double> getSpherical(double, double, double);
std::vector<double> inert2LocalTangent(std::vector<double>, double, double, double);
std::vector<double> localTangent2Inert(std::vector<double>, double, double, double);
std::vector<double> azEl2LocalTangent(double, double, double);
/** \} */

/**
 *  \name 2BP Utility Functions
 *  \{
 */
void r2bp_computeAllKepler(BaseArcset*);
void r2bp_computeKepler(const SysData_2bp*, Node*);
std::vector<double> r2bp_stateFromKepler(const SysData_2bp*, double, double, double, double, double, double);
/** \} */

/**
 *  \name CR3BP Utility Functions
 *  \{
 */
double cr3bp_getVel_withC(const double s[], double, double, int);
Traj_cr3bp cr3bp_EM2SE(Traj_cr3bp, const SysData_cr3bp*, double, double, double);
Nodeset_cr3bp cr3bp_EM2SE(Nodeset_cr3bp, const SysData_cr3bp*, double, double, double);
Traj_cr3bp cr3bp_SE2EM(Traj_cr3bp, const SysData_cr3bp*, double, double, double);
Nodeset_cr3bp cr3bp_SE2EM(Nodeset_cr3bp, const SysData_cr3bp*, double, double, double);
std::vector<double> cr3bp_EM2SE_state(std::vector<double>, double, double, double, double,
	double, double, double, double, double);
std::vector<double> cr3bp_SE2EM_state(std::vector<double>, double, double, double, double,
	double, double, double, double, double);
Nodeset_cr3bp cr3bp_rot2inert(Nodeset_cr3bp, double, int);
Traj_cr3bp cr3bp_rot2inert(Traj_cr3bp, double, int);
std::vector<double> cr3bp_rot2inert_state(std::vector<double>, const SysData_cr3bp*, double, double, int);
Traj_cr3bp cr3bp_getPeriodic(const SysData_cr3bp*, std::vector<double>, double, Mirror_tp, double tol = 1e-12);
Traj_cr3bp cr3bp_getPeriodic(const SysData_cr3bp*, std::vector<double>, double, int, int, Mirror_tp, std::vector<int>, double tol = 1e-12, MultShootData *pItData = nullptr);
/** \} */

/**
 *  \name BC4BP Utility Functions
 *  \{
 */
Traj_bc4bp bcr4bpr_SE2SEM(Traj_cr3bp, const SysData_bc4bp*, int, double);
Nodeset_bc4bp bcr4bpr_SE2SEM(Nodeset_cr3bp, const SysData_bc4bp*, int, double);
Traj_cr3bp bcr4bpr_SEM2SE(Traj_bc4bp, const SysData_cr3bp*);
Nodeset_cr3bp bcr4bpr_SEM2SE(Nodeset_bc4bp, const SysData_cr3bp*);
MatrixXRd bcr4bpr_spLoc_polyFit(const SysData_bc4bp*, double);
Eigen::Vector3d bcr4bpr_getSPLoc(const SysData_bc4bp*, double);
/** \} */
}// END of Astrohelion namespace

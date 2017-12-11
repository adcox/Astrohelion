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
class Arcset;
class Arcset_bc4bp;
class Arcset_cr3bp;
class Arcset_periodic;
class SysData_2bp;
class SysData_bc4bp;
class SysData_cr3bp;
class Arcset;
class Arcset_2bp;
class Arcset_bc4bp;
class Arcset_cr3bp;

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
 *  \name General Analysis Functions
 *  \{
 */
void balanceMat(MatrixXRd&, unsigned int&, unsigned int&, std::vector<double>&);
double dateToEphemerisTime(const char*);
double dateToGST(double, double, double, double, double, double);
void eigVec_backTrans(unsigned int, unsigned int, const std::vector<double>&, MatrixXRcd&);
void exchange(unsigned int, unsigned int, std::vector<double>&, unsigned int, unsigned int, MatrixXRd&);
std::vector<cdouble> getBalancedEigData(MatrixXRd, MatrixXRcd *pVecs = nullptr);
double gregorianToJulian(double, double, double, double, double, double);
MatrixXRd getMirrorMat(Mirror_tp);
double getStabilityIndex(std::vector<cdouble>);
Node interpPointAtTime(const Arcset*, double);
std::vector<unsigned int> sortEig(std::vector<cdouble>, std::vector<MatrixXRcd>);
void reconstructArc(const Arcset*, Arcset*);
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
 *  \name 2BP Analysis Functions
 *  \{
 */
void r2bp_computeAllKepler(BaseArcset*);
void r2bp_computeKepler(const SysData_2bp*, Node*);
std::vector<double> r2bp_stateFromKepler(const SysData_2bp*, double, double, double, double, double, double);
/** \} */

/**
 *  \name CR3BP Analysis Functions
 *  \{
 */
double cr3bp_getVel_withC(const double s[], double, double, int);
Arcset_cr3bp cr3bp_EM2SE(Arcset_cr3bp, const SysData_cr3bp*, double, double, double);
Arcset_cr3bp cr3bp_SE2EM(Arcset_cr3bp, const SysData_cr3bp*, double, double, double);
std::vector<double> cr3bp_EM2SE_state(std::vector<double>, double, double, double, double,
	double, double, double, double, double);
std::vector<double> cr3bp_SE2EM_state(std::vector<double>, double, double, double, double,
	double, double, double, double, double);
Arcset_cr3bp cr3bp_rot2inert(Arcset_cr3bp, double, int);
std::vector<double> cr3bp_rot2inert_state(std::vector<double>, const SysData_cr3bp*, double, double, int);

void cr3bp_addMirrorCons(Arcset_cr3bp*, Mirror_tp, std::vector<unsigned int>);
Arcset_periodic cr3bp_correctHalfPerSymPO(const Arcset_cr3bp*, Arcset_cr3bp*, Mirror_tp, double tol = 1e-12, MultShootData *pItData = nullptr);
void cr3bp_halfPO2fullPO(const Arcset_cr3bp*, Arcset_periodic*, Mirror_tp);
Arcset_cr3bp cr3bp_propHalfPerSymPO(const SysData_cr3bp*, std::vector<double>, double, unsigned int, unsigned int, Mirror_tp, std::vector<unsigned int>, double tol = 1e-12);
/** \} */

/**
 *  \name BC4BP Analysis Functions
 *  \{
 */
Arcset_bc4bp bcr4bpr_SE2SEM(Arcset_cr3bp, const SysData_bc4bp*, int, double);
Arcset_cr3bp bcr4bpr_SEM2SE(Arcset_bc4bp, const SysData_cr3bp*);
MatrixXRd bcr4bpr_spLoc_polyFit(const SysData_bc4bp*, double);
Eigen::Vector3d bcr4bpr_getSPLoc(const SysData_bc4bp*, double);
/** \} */
}// END of Astrohelion namespace

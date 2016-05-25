/**
 *	@brief Primarily contains functions used in numerical integration, like EOMs
 *
 *	This file contains functions used by the simulation_engine object 
 *	(and possibly others), like EOMs and other very "math-y" computations.
 *	Think of this file as a utility library for other classes
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 *	@copyright GNU GPL v3.0
 */
 
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
#ifndef H_CALCULATIONS
#define H_CALCULATIONS

#include "tpat_constants.hpp"
#include "tpat_eigen_defs.hpp"
 
#include <vector>

// Forward declarations
class TPAT_Correction_Engine;
class TPAT_MultShoot_Data;
class TPAT_Nodeset;
class TPAT_Nodeset_BC4BP;
class TPAT_Nodeset_CR3BP;
class TPAT_Sys_Data_BC4BP;
class TPAT_Sys_Data_CR3BP;
class TPAT_Traj;
class TPAT_Traj_BC4BP;
class TPAT_Traj_CR3BP;

/**
 *	@brief Describes the plane a periodic orbit can be mirrored across
 */
enum class TPAT_Mirror_Tp{
	MIRROR_XZ,		//!< Mirror over the XZ-Plane; x, z, and y-dot can be fixed if desired
	MIRROR_XY,		//!< Mirror over the XY-Plane; x, y, and z-dot can be fixed if desired
	MIRROR_YZ,		//!< Mirror over the YZ-Plane; y, z, and x-dot can be fixed if desired
	MIRROR_X_AX_H, 	//!< Mirror over the X-Axis, orbit is mostly in xy-plane; x, y-dot, and z-dot can be fixed if desired
	MIRROR_X_AX_V	//!< Mirror over the X-Axis, orbit is mostly in xz-plane; x, y-dot, and z-dot can be fixed if desired
};

/**
 *	@brief Describes the type of manifold, both stability and direction
 */
enum class TPAT_Manifold_Tp{
	Man_U_P,	//!< Unstable, departing towards +x direction
	Man_U_M,	//!< Unstable, departing towards -x direction
	Man_S_P,	//!< Stable, arriving from +x direction
	Man_S_M		//!< Stable, arriving from -x direction
};

/**
 *	@brief Eigenvalue pair types
 *
 *	Eigenvalues come in three different types of pairs. They will
 *	either be complex, real, or exactly equal to 1.0. Since the 
 *	matrices are real, all complex eigenvalues will come in conjucate
 *	pairs. In this type of problem (trajectory design), the other STM
 *	eigenvalues also come in pairs: 2 1.0 eigenvalues, or two real 
 *	eigenvalues that are reciprocals.
 */
enum class TPAT_EigValSet_Tp{
	EIGSET_COMP_CONJ,	//!< Complex conjugate pair
	EIGSET_ONES,		//!< Exactly equal to 1.0
	EIGSET_REAL_RECIP	//!< Real, reciprocal pair
};

// General Utility Functions
double dateToEpochTime(const char*);
std::vector<double> familyCont_LS(int, double, std::vector<int>, std::vector<double>);
std::vector<TPAT_Traj_CR3BP> getManifolds(TPAT_Manifold_Tp, const TPAT_Traj_CR3BP*, int, double);
MatrixXRd getMirrorMat(TPAT_Mirror_Tp);
double getStabilityIndex(std::vector<cdouble>);
double getTotalDV(const TPAT_MultShoot_Data*);
void finiteDiff_checkMultShoot(const TPAT_Nodeset*);
void finiteDiff_checkMultShoot(const TPAT_Nodeset*, TPAT_Correction_Engine);
MatrixXRd solveAX_eq_B(MatrixXRd, MatrixXRd);
std::vector<cdouble> sortEig(std::vector<cdouble>, std::vector<int>*);

// CR3BP Utility Functions
double cr3bp_getVel_withC(const double s[], double, double, int);
TPAT_Traj_CR3BP cr3bp_EM2SE(TPAT_Traj_CR3BP, const TPAT_Sys_Data_CR3BP*, double, double, double);
TPAT_Nodeset_CR3BP cr3bp_EM2SE(TPAT_Nodeset_CR3BP, const TPAT_Sys_Data_CR3BP*, double, double, double, double);
TPAT_Traj_CR3BP cr3bp_SE2EM(TPAT_Traj_CR3BP, const TPAT_Sys_Data_CR3BP*, double, double, double);
TPAT_Nodeset_CR3BP cr3bp_SE2EM(TPAT_Nodeset_CR3BP, const TPAT_Sys_Data_CR3BP*, double, double, double, double);
std::vector<double> cr3bp_EM2SE_state(std::vector<double>, double, double, double, double,
	double, double, double, double, double);
std::vector<double> cr3bp_SE2EM_state(std::vector<double>, double, double, double, double,
	double, double, double, double, double);
TPAT_Nodeset_CR3BP cr3bp_rot2inert(TPAT_Nodeset_CR3BP, int);
TPAT_Traj_CR3BP cr3bp_rot2inert(TPAT_Traj_CR3BP, int);
std::vector<double> cr3bp_rot2inert_state(std::vector<double>, const TPAT_Sys_Data_CR3BP*, double, int);
TPAT_Traj_CR3BP cr3bp_getPeriodic(const TPAT_Sys_Data_CR3BP*, std::vector<double>, double, TPAT_Mirror_Tp, double);
TPAT_Traj_CR3BP cr3bp_getPeriodic(const TPAT_Sys_Data_CR3BP*, std::vector<double>, double, int, int, TPAT_Mirror_Tp, std::vector<int>, double);
TPAT_Traj_CR3BP cr3bp_getPeriodic(const TPAT_Sys_Data_CR3BP*, std::vector<double>, double, int, int, TPAT_Mirror_Tp, std::vector<int>, double, TPAT_MultShoot_Data*);

// BCR4BPR Utility Functions
TPAT_Traj_BC4BP bcr4bpr_SE2SEM(TPAT_Traj_CR3BP, const TPAT_Sys_Data_BC4BP*, double);
TPAT_Nodeset_BC4BP bcr4bpr_SE2SEM(TPAT_Nodeset_CR3BP, const TPAT_Sys_Data_BC4BP*, double);
TPAT_Traj_BC4BP bcr4bpr_SEM2SE(TPAT_Traj_BC4BP, const TPAT_Sys_Data_CR3BP*);
TPAT_Nodeset_BC4BP bcr4bpr_SEM2SE(TPAT_Nodeset_BC4BP, const TPAT_Sys_Data_CR3BP*);
MatrixXRd bcr4bpr_spLoc_polyFit(const TPAT_Sys_Data_BC4BP*, double);
Eigen::Vector3d bcr4bpr_getSPLoc(const TPAT_Sys_Data_BC4BP*, double);

#endif
//END
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
class tpat_multShoot_data;
class tpat_nodeset;
class tpat_nodeset_bcr4bp;
class tpat_nodeset_cr3bp;
class tpat_sys_data_bcr4bpr;
class tpat_sys_data_cr3bp;
class tpat_traj;
class tpat_traj_bcr4bp;
class tpat_traj_cr3bp;

/**
 *	@brief Describes the plane a periodic orbit can be mirrored across
 */
enum tpat_mirror_tp{
	MIRROR_XZ,		//!< Mirror over the XZ-Plane; x, z, and y-dot can be fixed if desired
	MIRROR_XY,		//!< Mirror over the XY-Plane; x, y, and z-dot can be fixed if desired
	MIRROR_YZ,		//!< Mirror over the YZ-Plane; y, z, and x-dot can be fixed if desired
	MIRROR_X_AX_H, 	//!< Mirror over the X-Axis, orbit is mostly in xy-plane; x, y-dot, and z-dot can be fixed if desired
	MIRROR_X_AX_V	//!< Mirror over the X-Axis, orbit is mostly in xz-plane; x, y-dot, and z-dot can be fixed if desired
};

/**
 *	@brief Describes the type of manifold, both stability and direction
 */
enum tpat_manifold_tp{
	MAN_U_P,	//!< Unstable, departing towards +x direction
	MAN_U_M,	//!< Unstable, departing towards -x direction
	MAN_S_P,	//!< Stable, arriving from +x direction
	MAN_S_M		//!< Stable, arriving from -x direction
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
enum tpat_eigValSet_tp{
	EIGSET_COMP_CONJ,	//!< Complex conjugate pair
	EIGSET_ONES,		//!< Exactly equal to 1.0
	EIGSET_REAL_RECIP	//!< Real, reciprocal pair
};

// General Utility Functions
double dateToEpochTime(const char*);
std::vector<double> familyCont_LS(int, double, std::vector<int>, std::vector<double>);
std::vector<tpat_traj_cr3bp> getManifolds(tpat_manifold_tp, const tpat_traj_cr3bp*, int, double);
MatrixXRd getMirrorMat(tpat_mirror_tp);
double getStabilityIndex(std::vector<cdouble>);
double getTotalDV(const tpat_multShoot_data*);
void finiteDiff_checkMultShoot(const tpat_nodeset*);
MatrixXRd solveAX_eq_B(MatrixXRd, MatrixXRd);
std::vector<cdouble> sortEig(std::vector<cdouble>, std::vector<int>*);

// CR3BP Utility Functions
double cr3bp_getVel_withC(const double s[], double, double, int);
tpat_traj_cr3bp cr3bp_EM2SE(tpat_traj_cr3bp, const tpat_sys_data_cr3bp*, double, double, double);
tpat_nodeset_cr3bp cr3bp_EM2SE(tpat_nodeset_cr3bp, const tpat_sys_data_cr3bp*, double, double, double, double);
tpat_traj_cr3bp cr3bp_SE2EM(tpat_traj_cr3bp, const tpat_sys_data_cr3bp*, double, double, double);
tpat_nodeset_cr3bp cr3bp_SE2EM(tpat_nodeset_cr3bp, const tpat_sys_data_cr3bp*, double, double, double, double);
std::vector<double> cr3bp_EM2SE_state(std::vector<double>, double, double, double, double,
	double, double, double, double, double);
std::vector<double> cr3bp_SE2EM_state(std::vector<double>, double, double, double, double,
	double, double, double, double, double);
tpat_nodeset_cr3bp cr3bp_rot2inert(tpat_nodeset_cr3bp, int);
tpat_traj_cr3bp cr3bp_rot2inert(tpat_traj_cr3bp, int);
std::vector<double> cr3bp_rot2inert_state(std::vector<double>, const tpat_sys_data_cr3bp*, double, int);
tpat_traj_cr3bp cr3bp_getPeriodic(const tpat_sys_data_cr3bp*, std::vector<double>, double, tpat_mirror_tp, double);
tpat_traj_cr3bp cr3bp_getPeriodic(const tpat_sys_data_cr3bp*, std::vector<double>, double, int, int, tpat_mirror_tp, std::vector<int>, double);
tpat_traj_cr3bp cr3bp_getPeriodic(const tpat_sys_data_cr3bp*, std::vector<double>, double, int, int, tpat_mirror_tp, std::vector<int>, double, tpat_multShoot_data*);

// BCR4BPR Utility Functions
tpat_traj_bcr4bp bcr4bpr_SE2SEM(tpat_traj_cr3bp, const tpat_sys_data_bcr4bpr*, double);
tpat_nodeset_bcr4bp bcr4bpr_SE2SEM(tpat_nodeset_cr3bp, const tpat_sys_data_bcr4bpr*, double);
tpat_traj_bcr4bp bcr4bpr_SEM2SE(tpat_traj_bcr4bp, const tpat_sys_data_cr3bp*);
tpat_nodeset_bcr4bp bcr4bpr_SEM2SE(tpat_nodeset_bcr4bp, const tpat_sys_data_cr3bp*);
MatrixXRd bcr4bpr_spLoc_polyFit(const tpat_sys_data_bcr4bpr*, double);
Eigen::Vector3d bcr4bpr_getSPLoc(const tpat_sys_data_bcr4bpr*, double);

#endif
//END
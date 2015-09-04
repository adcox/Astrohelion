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
 
#include <vector>

// Forward declarations
class tpat_nodeset_bcr4bpr;
class tpat_sys_data_bcr4bpr;
class tpat_traj;
class tpat_traj_bcr4bpr;
class tpat_nodeset_cr3bp;
class tpat_sys_data_cr3bp;
class tpat_traj_cr3bp;
class tpat_matrix;

/**
 *	@brief Describes the plane a periodic orbit can be mirrored across
 */
enum mirror_t{
	MIRROR_XZ,		//!< Mirror over the XZ-Plane; x, z, and y-dot can be fixed if desired
	MIRROR_XY,		//!< Mirror over the XY-Plane; x, y, and z-dot can be fixed if desired
	MIRROR_YZ,		//!< Mirror over the YZ-Plane; y, z, and x-dot can be fixed if desired
	MIRROR_X_AX_H, 	//!< Mirror over the X-Axis, orbit is mostly in xy-plane; x, y-dot, and z-dot can be fixed if desired
	MIRROR_X_AX_V	//!< Mirror over the X-Axis, orbit is mostly in xz-plane; x, y-dot, and z-dot can be fixed if desired
};

/**
 *	@brief Describes the type of manifold, both stability and direction
 */
enum manifold_t{
	MAN_U_P,	//!< Unstable, departing towards +x direction
	MAN_U_M,	//!< Unstable, departing towards -x direction
	MAN_S_P,	//!< Stable, arriving from +x direction
	MAN_S_M		//!< Stable, arriving from -x direction
};

// Equations of motion
int cr3bp_EOMs(double, const double[], double[], void*);
int cr3bp_simple_EOMs(double, const double[], double[], void*);

int cr3bp_ltvp_EOMs(double, const double[], double[], void*);
int cr3bp_ltvp_simple_EOMs(double, const double[], double[], void*);

int bcr4bpr_EOMs(double, const double[], double[], void*);
int bcr4bpr_simple_EOMs(double, const double[], double[], void*);

// General Utility Functions
double dateToEpochTime(const char*);
std::vector<double> familyCont_LS(int, double, std::vector<int>, std::vector<double>);
tpat_matrix solveAX_eq_B(tpat_matrix, tpat_matrix);
std::vector<tpat_traj_cr3bp> getManifolds(manifold_t, tpat_traj_cr3bp*, int, double);
std::vector<cdouble> sortEig(std::vector<cdouble>, std::vector<int>*);

// CR3BP Utility Functions
void cr3bp_getUDDots(double, double, double, double, double*);
double cr3bp_getJacobi(const double s[], double);
void cr3bp_getEquilibPt(tpat_sys_data_cr3bp, int, double, double[3]);
tpat_traj_cr3bp cr3bp_EM2SE(tpat_traj_cr3bp, tpat_sys_data_cr3bp*, double, double, double);
tpat_nodeset_cr3bp cr3bp_EM2SE(tpat_nodeset_cr3bp, tpat_sys_data_cr3bp*, double, double, double, double);
tpat_traj_cr3bp cr3bp_SE2EM(tpat_traj_cr3bp, tpat_sys_data_cr3bp*, double, double, double);
tpat_nodeset_cr3bp cr3bp_SE2EM(tpat_nodeset_cr3bp, tpat_sys_data_cr3bp*, double, double, double, double);
std::vector<double> cr3bp_EM2SE_state(std::vector<double>, double, double, double, double,
	double, double, double, double, double);
std::vector<double> cr3bp_SE2EM_state(std::vector<double>, double, double, double, double,
	double, double, double, double, double);
tpat_traj_cr3bp cr3bp_getPeriodic(tpat_sys_data_cr3bp, std::vector<double>, double, mirror_t);
tpat_traj_cr3bp cr3bp_getPeriodic(tpat_sys_data_cr3bp, std::vector<double>, double, int, mirror_t, std::vector<int>);

// BCR4BPR Utility Functions
void bcr4bpr_getPrimaryPos(double, tpat_sys_data_bcr4bpr*, double*);
void bcr4bpr_getPrimaryVel(double, tpat_sys_data_bcr4bpr*, double*);
tpat_traj_bcr4bpr bcr4bpr_SE2SEM(tpat_traj_cr3bp, tpat_sys_data_bcr4bpr*, double);
tpat_nodeset_bcr4bpr bcr4bpr_SE2SEM(tpat_nodeset_cr3bp, tpat_sys_data_bcr4bpr*, double);
void bcr4bpr_orientAtEpoch(double, tpat_sys_data_bcr4bpr*);

#endif
//END
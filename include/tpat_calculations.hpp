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
#ifndef __H_CALCULATIONS__
#define __H_CALCULATIONS__

#include <vector>

// Forward declarations
class tpat_bcr4bpr_nodeset;
class tpat_bcr4bpr_sys_data;
class tpat_bcr4bpr_traj;
class tpat_cr3bp_nodeset;
class tpat_cr3bp_sys_data;
class tpat_cr3bp_traj;

// Equations of motion
int cr3bp_EOMs(double t, const double s[], double sdot[], void *params);
int cr3bp_simple_EOMs(double t, const double s[], double sdot[], void *params);

int bcr4bpr_EOMs(double t, const double s[], double sdot[], void *params);
int bcr4bpr_simple_EOMs(double t, const double s[], double sdot[], void *params);

// General Utility Functions
double dateToEpochTime(const char*);

// CR3BP Utility Functions
void cr3bp_getUDDots(double, double, double, double, double*);
double cr3bp_getJacobi(double s[], double);
void cr3bp_getEquilibPt(tpat_cr3bp_sys_data, int, double, double[3]);
tpat_cr3bp_traj cr3bp_EM2SE(tpat_cr3bp_traj, double, double, double);
tpat_cr3bp_nodeset cr3bp_EM2SE(tpat_cr3bp_nodeset, double, double, double, double);
tpat_cr3bp_traj cr3bp_SE2EM(tpat_cr3bp_traj, double, double, double);
tpat_cr3bp_nodeset cr3bp_SE2EM(tpat_cr3bp_nodeset, double, double, double, double);
std::vector<double> cr3bp_EM2SE_state(std::vector<double>, double, double, double, double,
	double, double, double, double, double);
std::vector<double> cr3bp_SE2EM_state(std::vector<double>, double, double, double, double,
	double, double, double, double, double);

// BCR4BPR Utility Functions
void bcr4bpr_getPrimaryPos(double, tpat_bcr4bpr_sys_data, double*);
void bcr4bpr_getPrimaryVel(double, tpat_bcr4bpr_sys_data, double*);
tpat_bcr4bpr_traj bcr4bpr_SE2SEM(tpat_cr3bp_traj, tpat_bcr4bpr_sys_data, double);
tpat_bcr4bpr_nodeset bcr4bpr_SE2SEM(tpat_cr3bp_nodeset, tpat_bcr4bpr_sys_data, double);
void bcr4bpr_orientAtEpoch(double, tpat_bcr4bpr_sys_data*);

#endif
//END
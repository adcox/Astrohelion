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
 *	This file is part of the Astrodynamics Toolkit (ADTK).
 *
 *  ADTK is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ADTK is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ADTK.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __H_CALCULATIONS__
#define __H_CALCULATIONS__

// #include "adtk_bcr4bpr_sys_data.hpp"

// Forward declarations
class adtk_bcr4bpr_sys_data;
class adtk_cr3bp_sys_data;

// Equations of motion
int cr3bp_EOMs(double t, const double s[], double sdot[], void *params);
int cr3bp_simple_EOMs(double t, const double s[], double sdot[], void *params);

int bcr4bpr_EOMs(double t, const double s[], double sdot[], void *params);
int bcr4bpr_simple_EOMs(double t, const double s[], double sdot[], void *params);

// CR3BP Utility Functions
void cr3bp_getUDDots(double, double, double, double, double*);
double cr3bp_getJacobi(double s[], double);
void cr3bp_getEquilibPt(adtk_cr3bp_sys_data, int, double, double[3]);

// BCR4BPR Utility Functions
void bcr4bpr_getPrimaryPos(double, adtk_bcr4bpr_sys_data, double*);
void bcr4bpr_getPrimaryVel(double, adtk_bcr4bpr_sys_data, double*);

#endif
//END
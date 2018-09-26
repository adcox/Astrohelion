/**
 *	@file Common.hpp
 *	@brief Contains values for physical constants like pi, G, etc. and custom data wrappers
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 *	@copyright GNU GPL v3.0
 */

/*
 *	Astrohelion 
 *	Copyright 2015-2018, Andrew Cox; Protected under the GNU GPL v3.0
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

#include <complex>

namespace astrohelion{

/**
 * @brief Provides a scale larger than true/false to set verbosity for Core routines
 */
enum class Verbosity_tp : int{
	NO_MSG 		= 0,	//!< Absolutely no messages
	SOME_MSG	= 1,	//!< Some useful messages are printed to standard output
	ALL_MSG		= 2,	//!< All useful and not-so-useful messages are printed to standard output
	DEBUG 		= 3		//!< Even more messages; helpful for debugging
};

/**
 *  @brief Describe how much data is saved
 */
enum class Save_tp : int{
	SAVE_ALL = 0, 		//!< Save all data
	SAVE_CURVE = 1,		//!< Do not save STM and extra parameter values along segment evolutions
	SAVE_FRAME = 2		//!< Do not save any states along segment evolutions
};

// Custom types
typedef std::complex<double> cdouble;	//!< A complex double
typedef std::complex<int> cint;			//!< A complex integer

/** Universal Gravity Constant, km^3/kg-s^2 */
const double G = 6.67384e-20;

/** 1 G acceleration, km/s^2 */
const double G_GRAV_0 = 9.80665/1000;

/** Pi */
const double PI = 3.14159265358979323846264338327950;

static const unsigned int BUFFER_SIZE = 256;

// Matlab File Name Variables, Administrative
static const char VARNAME_TIMESTAMP[] = "timestamp";		//!< Matlab file variable name for a timestamp that indicates the veresion

// Matlab File Name Variables for 2BP
static const char VARNAME_2BP_ENERGY[] = "energy";			//!< Matlab file variable name for 2BP energy data
static const char VARNAME_ANGMOM[] = "angMom";				//!< Matlab file variable name for angular momentum data
static const char VARNAME_ECC[] = "ecc";					//!< Matlab file variable name for eccentricity data
static const char VARNAME_SMA[] = "sma";					//!< Matlab file variable name for semimajor axis data

// Matlab File Name Variables for General Arcset
static const char VARNAME_LINKTABLE[] = "LinkTable";		//!< Matlab file variable name for the segment link table
static const char VARNAME_CONSTRAINTS[] = "Constraints";	//!< Matlab file variable name for the constraints
static const char VARNAME_NODECTRL[] = "NodeCtrlState";		//!< Matlab file variable name for node control state
static const char VARNAME_NODESTATE[] = "NodeState";		//!< Matlab file variable name for node state vectors
static const char VARNAME_NODETIME[] = "NodeTime";			//!< Matlab file variable name for node time/epoch values
static const char VARNAME_SEGCTRL[] = "SegCtrl";			//!< Matlab file variable name for Segment control structure
static const char VARNAME_SEGSTATE[] = "SegState";			//!< Matlab file variable name for segment state vectors
static const char VARNAME_SEGTIME[] = "SegTime";			//!< Matlab file variable name for segment time vectors
static const char VARNAME_SEGTOF[] = "SegTOF";				//!< Matlab file variable name for segment time-of-flight values

static const char VARNAME_STATE_DERIV[] = "qdot";			//!< Matlab file variable name for state derivative data
static const char VARNAME_STM[] = "STM";					//!< Matlab file variable name for state transition matrix data
static const char VARNAME_TOF[] = "TOF";					//!< Matlab file variable name for time-of-flight data

// Matlab File Name Variables for Specific Arcsets
static const char VARNAME_CTRL_LAW[] = "CtrlLaw";			//!< Matlab file variable name for control law data
static const char VARNAME_JACOBI[] = "Jacobi";				//!< Matlab file variable name for Jacobi constant data
static const char VARNAME_STATE_EPOCH_DERIV[] = "dqdT";		//!< Matlab file variable name for state derivative w.r.t. epoch

// Deprecated Matlab file variable names
static const char VARNAME_DEP_EPOCH[] = "Epoch";			//!< Matlab file variable name for epoch data
static const char VARNAME_DEP_NODE[] = "Node";				//!< Matlab file variable name for node data
static const char VARNAME_DEP_STATE[] = "State";			//!< Matlab file variable name for state data
static const char VARNAME_DEP_TIME[] = "Time";				//!< Matlab file variable name for time data

// ExtraParam Map Keys
static const char PARAMKEY_2BP_ENERGY[] = "energy";			//!< ExtraParam map key for 2BP energy data
static const char PARAMKEY_ANGMOM[] = "angMom";				//!< ExtraParam map key for angular momentum data
static const char PARAMKEY_CTRL[] = "ctrl";					//!< ExtraParam map key for control vector
static const char PARAMKEY_ECC[] = "ecc";					//!< ExtraParam map key for eccentricity data
static const char PARAMKEY_JACOBI[] = "J";					//!< ExtraParam map key for Jacobi constant data
static const char PARAMKEY_SMA[] = "sma";					//!< ExtraParam map key for semimajor axis data
static const char PARAMKEY_STATE_DERIV[] = "qdot";			//!< ExtraParam map key for state derivative data
static const char PARAMKEY_STATE_EPOCH_DERIV[] = "dqdT";	//!< ExtraParam map key for state derivative w.r.t. epoch 


// Timestamps for file versioning
static const char T_V1_9[] = "19:15:00UTC 02/02/2018";		// Switch all low-thrust magnitudes to store as sqrt(f)
static const char T_V2_0[] = "00:00:00UTC 26/02/2018";		// Introduced timestamps to files
static const char T_V2_1[] = "21:54:00UTC 03/04/2018";		// Changed low-thrust control to new ID with various thrust parameterizations
}// END of Astrohelion namespace
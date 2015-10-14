/** 
 *	@file tpat.hpp
 *
 *	@brief Contains functions and variables that are global for the entire library
 *
 *	Note to end user: You do NOT need to include this header
 *
 *	For developers: This header should be included on ALL library .cpp files to ensure that the library
 *	is fully initialized no matter how the user uses the library. A static <tt>initializer</tt>
 *	object is instantiated as a global variable so that it will be constructed when the 
 *	application starts and deconstructed when the application is finished. The include
 *	guard keeps the initializion from happening more than once.
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

#ifndef H_TPAT_INIT
#define H_TPAT_INIT

#include <cspice/SpiceZfc.h>	// prototypes for functions
#include <cspice/SpiceZdf.h>	// typedefs for SPICE objects, like SpiceDouble
#include <gsl/gsl_errno.h>
 
#include <iostream>

/** 
 * 	What action SPICE should take when it encounters an error
 *	Options are:
 *	'abort' 	- 	Designed for safety; output messages and traceback to
 *					screen or stdout, then stop program and return status
 *					code if possible.
 *
 *	'return'	-	For use in programs that must keep running (our default).
 *					Output messages and traceback are saved to current error
 *					device and can be retrieved for error handling.
 */
static char SPICE_ERR_ACTION[] = "return";

/**
 *	The type of message type SPICE should return to us
 */
static char SPICE_ERR_MSG_TYPE[] = "short,traceback";

/**
 *	Flag to prevent initializations happing lots of times
 */
static bool isInitialized = false;

/**
 *	@brief a structure that contains functions to initialize and unload 
 *	global library stuff
 */
struct initializer {
	/**
	 *	@brief Initialize the library
	 */
	initializer(){
		if(!isInitialized){
			std::cout << "Initializing TPAT System" << std::endl;

			// Tell SPICE how to handle errors
		    erract_c("set", 0, SPICE_ERR_ACTION);
		    errprt_c("set", 0, SPICE_ERR_MSG_TYPE);

		    // Turn GSL's error handler off; we will catch and handle the errors
		    gsl_set_error_handler_off();

		    isInitialized = true;
		}
	}

	/**
	 *	@brief Perform any duties necessary to safely shut down the library
	 */
	~initializer(){
		if(isInitialized){
			std::cout << "Closing out the TPAT System" << std::endl;
			isInitialized = false;
		}
	}
};

static initializer tpat_init;	// Create one to initialize the system

#endif
/** 
 *	@file Core.hpp
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
 *	Astrohelion 
 *	Copyright 2016, Andrew Cox; Protected under the GNU GPL v3.0
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

/**
 * @defgroup engine Engines
 * Contains all engine objects and methods
 */ 

/** @defgroup traj Trajectory Representations
 * Contains objects that represent trajectories in one way or another
 */

/** @defgroup model Model
 * Contains objects that represent models of physical objects and systems
 */

/** @defgroup util Utility
 * Contains all utility objects and methods
 */

/** @defgroup 2bp 2BP
 * Contains objects and methods that handle 2BP trajectories
 */ 

/** @defgroup cr3bp CR3BP
 * Contains objects and methods that handle CR3BP trajectories
 */ 

/** @defgroup bc4bp BC4BP
 * Contains objects and methods that handle BC4BP trajectories
 */

#pragma once

#include <string>

/**
 * @brief A namespace to contain all objects that are part of
 * the Astrohelion library; prevent clashes with other libraries
 */
namespace astrohelion{

/**
 * Type that represents the OS this code is running/compiled on
 */
enum class OS_tp : int {windows=0, apple=1, linux=2};

/**
 * Type that defines the class of operating system
 */
static const OS_tp OS_TYPE =
#ifdef _WIN32
 	// Includes 32-bit and 64-bit windows
	OS_tp::windows;
#elif __APPLE__
	OS_tp::apple;
#elif __linux__
	OS_tp::linux;
#endif

/**
 *  @brief A structure that contains default settings and functions to load
 *  settings from an XML file
 */
struct Core_Settings{
	
	void load(const std::string&);
	void save(const std::string&);

	/** Absolute path to a directory containing SPICE kernels used by this software */
	std::string spice_data_filepath = "~/SPICE/data/";

	/** Time kernel to use */
	std::string spice_time_kernel = "naif0010.tls.pc";
}; // END OF Core_Settings

/**
 *	@brief a structure that contains functions to initialize and unload 
 *	global library stuff
 */
struct Core_Initializer {

	Core_Initializer();
	~Core_Initializer();

	void runInit();

	Core_Settings settings;	//!< A collection of settings that various Core algorithms use
};

/**
 *  @brief An object that serves as a parent to all Core objects.
 *  @details This object's sole purpose is to load settings from
 *  an XML file when any of the Core objects are instantiated. Static
 *  objects are used to ensure that only one copy of the initializer
 *  and settings structures are constructed
 */
class Core{
public:
	static Core_Initializer initializer;				//!< Create one to initialize the system	
	static bool bIsInit;								//!< Flag to prevent initializations happing lots of times

	Core();												//!< Default constructor
	virtual ~Core();									//!< Default destructor
};


}// END of Astrohelion namespace

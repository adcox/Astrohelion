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

#include <string>

/**
 * Type that represents the OS this code is running/compiled on
 */
enum class TPAT_OS_Type : int {windows=0, apple=1, linux=2};

/**
 * Type that defines the class of operating system
 */
static const TPAT_OS_Type TPAT_OS_TYPE =
#ifdef _WIN32
 	// Includes 32-bit and 64-bit windows
	TPAT_OS_Type::windows;
#elif __APPLE__
	TPAT_OS_Type::apple;
#elif __linux__
	TPAT_OS_Type::linux;
#endif

/**
 *  @brief A structure that contains default settings and functions to load
 *  settings from an XML file
 */
struct tpat_settings{
	
	void load(const std::string&);
	void save(const std::string&);

	/** Absolute path to a directory containing SPICE kernels used by this software */
	std::string spice_data_filepath = "~/SPICE/data/";

	/** Time kernel to use */
	std::string spice_time_kernel = "naif0010.tls.pc";
}; // END OF TPAT_SETTINGS

/**
 *	@brief a structure that contains functions to initialize and unload 
 *	global library stuff
 */
struct tpat_initializer {

	tpat_initializer();
	~tpat_initializer();

	void runInit();

	tpat_settings settings;	//!< A collection of settings that various TPAT algorithms use
};

/**
 *  @brief An object that serves as a parent to all TPAT objects.
 *  @details This object's sole purpose is to load settings from
 *  an XML file when any of the TPAT objects are instantiated. Static
 *  objects are used to ensure that only one copy of the initializer
 *  and settings structures are constructed
 */
class TPAT{
public:
	static tpat_initializer initializer;				//!< Create one to initialize the system	
	static bool isInit;									//!< Flag to prevent initializations happing lots of times

	TPAT();												//!< Default constructor
	virtual ~TPAT();									//!< Default destructor
};


#endif
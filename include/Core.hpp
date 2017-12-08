/** 
 *	\file Core.hpp
 *
 *	\brief Contains functions and variables that are global for the entire library
 *
 *	Note to end user: You do NOT need to include this header
 *
 *	For developers: This header should be included on ALL library .cpp files to ensure that the library
 *	is fully initialized no matter how the user uses the library. A static `initializer`
 *	object is instantiated as a global variable so that it will be constructed when the 
 *	application starts and deconstructed when the application is finished. The include
 *	guard keeps the initializion from happening more than once.
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

/**
 * \defgroup engine Engines
 * Contains all engine objects and methods
 */ 

/** \defgroup traj Trajectory Representations
 * Contains objects that represent trajectories in one way or another
 */

/** \defgroup fam Trajectory Families
 * Contains objects the create, store, or operate on families
 */

/** \defgroup model Model
 * Contains objects that represent models of physical objects and systems
 */

/** \defgroup util Utility
 * Contains all utility objects and methods
 */

/** \defgroup 2bp 2BP
 * Contains objects and methods that handle 2BP data
 */ 

/** \defgroup cr3bp CR3BP
 * Contains objects and methods that handle CR3BP data
 */ 

/** \defgroup bc4bp BC4BP
 * Contains objects and methods that handle BC4BP data
 */

/** \defgroup cr3bp_lt CR3BP-LT
 * Contains objects and methods that handle low-thrust CR3BP data
 */

#pragma once

#define NDEBUG		// Dissable all asserts and other debugging features

/* base radix of the computer architecture; almost always 2 (binary),
 * but may be different for some unique systems.
 *
 *	Used in:
 *	- Calculations:balanceMat()
 */
#define RADIX 2	

#include <string>
#include <map>

/**
 * \brief A namespace to contain all objects that are part of
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

struct Body_Data{
	int id = 0;				//!< Unique ID (same as SPICE ID or HORIZONS ID) for this body
	double gravParam = 0;	//!< Gravitational parameter associated with this body, kg^3/s^2
	double bodyRad = 0;		//!< Mean radius of the body, km
	double mass = 0;		//!< Mass of the body, kg
	double orbitRad = 0;	//!< Mean orbital radius (distance from parent body), km

	/** Minmum acceptable fly-by altitude for this body. Altitudes lower than this will
	trigger a crash event in the numerical simulation */
	double minFlyByAlt = 0;

	std::string name = "nullptr";		//!< Name of this body
	std::string parent = "nullptr";	//!< Name of the parent body
}; // END OF BODY_DATA

/**
 *  \brief A structure that contains default settings and functions to load
 *  settings from an XML file
 */
class Core_Settings{
public:
	void load(const std::string&);
	void save(const std::string&);

	/** Absolute path to a directory containing SPICE kernels used by this software */
	std::string spice_data_filepath = "~/SPICE/data/";

	/** Time kernel to use */
	std::string spice_time_kernel = "naif0010.tls.pc";

	/** SPK Kernel to use */
	std::string spice_spk_kernel = "de430.bsp";

	/** Variable to save settings filepath in after it is determined or created */
	std::string settings_filepath = "";

	/** Variable to save body data filepath in after it is determined or created */
	std::string body_data_filepath = "";
	
}; // END OF Core_Settings

/**
 *	\brief a structure that contains functions to initialize and unload 
 *	global library stuff
 */
class Core_Initializer {
public:
	Core_Initializer();
	~Core_Initializer();

	void runInit();

	Core_Settings settings;					//!< A collection of settings that various Core algorithms use
	std::map<int, Body_Data> allBodyData {};	//!< Map of Body_Data objects; key is the SPICE ID
};

/**
 *  \brief An object that serves as a parent to all Core objects.
 *  \details This object's sole purpose is to load settings from
 *  an XML file when any of the Core objects are instantiated. Static
 *  objects are used to ensure that only one copy of the initializer
 *  and settings structures are constructed
 */
class Core{
public:
	Core();												//!< Default constructor
	virtual ~Core();									//!< Default destructor

	/**
	 *  \brief Initialize the core
	 *  \details This function wraps a static Core_Initializer object to avoid multiple
	 *  instances of the initializer in different translation units
	 *  \return a reference to the initializer
	 */
	Core_Initializer& initializer(){
		static Core_Initializer initializer;
		return initializer;
	}

	/**
	 *  \brief Check to see if the core is initialized
	 *  \details This function wraps a static boolean object to avoid multiple 
	 *  instances of the variable in different translation units
	 *  \return whether or not the core is initialized
	 */
	bool& bIsInit(){
		static bool b_isInit;
		return b_isInit;
	}
};


}// END of Astrohelion namespace

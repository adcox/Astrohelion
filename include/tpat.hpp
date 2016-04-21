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

#include <exception> 
#include <iostream>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <cspice/SpiceZfc.h>	// prototypes for functions
#include <cspice/SpiceZdf.h>	// typedefs for SPICE objects, like SpiceDouble
#include <gsl/gsl_errno.h>

#include <pwd.h>

namespace fs = boost::filesystem;

/**
 * Type that represents the OS this code is running/compiled on
 */
enum tpat_os_type {windows, apple, linux};

/**
 * 
 */
static const tpat_os_type TPAT_OS_TYPE =
#ifdef _WIN32
 	// Includes 32-bit and 64-bit windows
	windows;
#elif __APPLE__
	apple;
#elif __linux__
	linux;
#endif

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

/** The type of message type SPICE should return to us */
static char SPICE_ERR_MSG_TYPE[] = "short,traceback";

struct tpat_settings{
	/** Absolute path to a directory containing SPICE kernels used by this software */
	std::string spice_data_filepath = "~/SPICE/data/";

	/** Time kernel to use */
	std::string spice_time_kernel = "naif0010.tls.pc";

	/**
	 *  @brief Load user settings from the prescribed file
	 *  @param filename path to the XML settings file
	 */
	void load(const std::string &filename){
		// Create empty property tree object
		using boost::property_tree::ptree;
		ptree pt;

		try{
			// Load XML file and put its contents in the property tree; if file isn't found, exception is thrown
			read_xml(filename, pt);
		}catch(std::exception &e){
			std::cout << "Error: Could not load tpat.spice.data_filepath from settings file\n";
			throw e;
		}

		try{
			// Read the setting values from the PTree; if they aren't found, an exception is thrown
			spice_data_filepath = pt.get<std::string>("tpat.spice.data_filepath");
			spice_time_kernel = pt.get<std::string>("tpat.spice.time_kernel");
		}catch(std::exception &e){
			std::cout << "Error: Could not load setting values from settings file\n";
			throw e;
		}
	}//==================================================

	/**
	 *  @brief Save user settings to the prescribed file
	 *  @param filename path to the XML settings file
	 */
	void save(const std::string &filename){
		// Create empty property tree object
		using boost::property_tree::ptree;
		ptree pt;

		// Put spice data filepath in property tre
		pt.put("tpat.spice.data_filepath", spice_data_filepath);
		pt.put("tpat.spice.time_kernel", spice_time_kernel);
		
		try{
			// Write property tree to XML file
			write_xml(filename, pt);
		}catch(std::exception &e){
			throw e;
		}
	}//==================================================
}; // END OF TPAT_SETTINGS

/**
 *	@brief a structure that contains functions to initialize and unload 
 *	global library stuff
 */
struct tpat_initializer {
	
	tpat_settings settings;	//!< A collection of settings that various TPAT algorithms use

	tpat_initializer(){
		std::cout << "Constructing tpat_initializer";
	}

	/**
	 *	@brief Initialize the library
	 */
	void runInit(){
		std::cout << "Initializing TPAT System" << std::endl;

		// Get user home folder
		uid_t uid = getuid();
		struct passwd *pw = getpwuid(uid);

		if(pw == NULL){
			std::cout << "Error in tpat_initializer: Could not get system user database\n";
			throw std::exception();
		}

		// Define names of settings files
		const char *defaultSettingsFile = "default_settings.xml";
		const char *userSettingsFile = "user_settings.xml";

		// Construct path to directory that contains the settings file(s)
		char settingsFilepath[128], defaultSettingsFilepath[256], userSettingsFilepath[256];
		switch(TPAT_OS_TYPE){
			case windows: 
				std::cout << "Error: Windows is not currently supported!\n";
				throw std::exception();
				break;
			case apple:
			case linux:
			default:
				sprintf(settingsFilepath, "%s%s", pw->pw_dir, "/.config/tpat/");
				break;
		}

		// Check to see if the directory exists; if it doesn't, create it.
		fs::path p = fs::system_complete(settingsFilepath);	// system_complete interprets any ~ or ../
		if(!fs::exists(p)){
			// Create the directory
			std::cout << "Settings directory does not exist; creating it now...\n";
			fs::create_directory(settingsFilepath);
		}

		// Construct full paths to the settings files
		sprintf(defaultSettingsFilepath, "%s%s", settingsFilepath, defaultSettingsFile);
		sprintf(userSettingsFilepath, "%s%s", settingsFilepath, userSettingsFile);

		// Load user settings from file
		try{
			settings.load(userSettingsFilepath);
			std::cout << "Loaded user settings\n";
		}catch(std::exception &e){
			// Wasn't able to load user settings, so try loading default settings
			std::cout << "Unable to load user settings, attempting to load default_settings\n";
			try{
				settings.load(defaultSettingsFilepath);
				std::cout << "Loaded default settings; creating user settings file at:\n\t" << userSettingsFilepath << std::endl;
			}catch(std::exception &ee){
				// Wasn't able to load default user settings either...
				std::cout << "Unable to load default settings\n";
				std::cout << "Creating default settings file at " << defaultSettingsFilepath << std::endl;
				std::cout << "Creating user settings file at " << userSettingsFilepath << std::endl;
				settings.save(defaultSettingsFilepath);	// Create a default settings file
				settings.save(userSettingsFilepath);	// Also create a user settings file, will be same as defaults
			}
			// If default settings were loaded, save a user settings file that is identical to default for future use/customization
			settings.save(userSettingsFilepath);
		}

		// Tell SPICE how to handle errors
	    erract_c("set", 0, SPICE_ERR_ACTION);
	    errprt_c("set", 0, SPICE_ERR_MSG_TYPE);

	    // Turn GSL's error handler off; we will catch and handle the errors
	    gsl_set_error_handler_off();
	}

	/**
	 *	@brief Perform any duties necessary to safely shut down the library
	 */
	~tpat_initializer(){
		std::cout << "Destructing tpat_initializer; Closing out the TPAT System" << std::endl;
	}
};

class tpat{
public:
	static tpat_initializer initializer;				// Create one to initialize the system	
	static bool isInit;									//!< Flag to prevent initializations happing lots of times

	tpat();
};


#endif
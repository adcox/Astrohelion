/**
 *  @file Core.cpp
 *	@brief Base class for all library objects that can exist independently
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */

#include <exception> 
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <cspice/SpiceZfc.h>	// prototypes for functions
#include <cspice/SpiceZdf.h>	// typedefs for SPICE objects, like SpiceDouble
#include <gsl/gsl_errno.h>		// custom error handling for GSL
#include <pwd.h>				// Get user ID and other system information
#include <thread>				// Hardware thread detection

#ifdef _OPENMP
	#include <omp.h>			// OpenMP header
#endif

#include "Core.hpp"

namespace fs = boost::filesystem;	// shortened for readability

namespace astrohelion{
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


//-----------------------------------------------------
//      CLASS Core FUNCTIONS
//-----------------------------------------------------

Core_Initializer Core::initializer;
bool Core::bIsInit;

/**
 *  @brief Default constructor.
 *  @details Runs an initialization sequence the first time the class is instantiated.
 */
Core::Core(){
	if(!bIsInit){
		initializer.runInit();
		bIsInit = true;
	}
}//======================================================

/**
 *  @brief Default destructor; doesn't do anything :)
 */
Core::~Core(){}


//-----------------------------------------------------
//      CLASS Core_Settings FUNCTIONS
//-----------------------------------------------------


/**
 *  @brief Load user settings from the prescribed file
 *  @param filename path to the XML settings file
 */
void Core_Settings::load(const std::string &filename){
	// Create empty property tree object
	using boost::property_tree::ptree;
	ptree pt;

	try{
		// Load XML file and put its contents in the property tree; if file isn't found, exception is thrown
		read_xml(filename, pt);
	}catch(std::exception &e){
		std::cout << "Error: Could not load the settings file\n";
		throw e;
	}

	try{
		// Read the setting values from the PTree; if they aren't found, an exception is thrown
		spice_data_filepath = pt.get<std::string>("Astrohelion.spice.data_filepath");
		spice_time_kernel = pt.get<std::string>("Astrohelion.spice.time_kernel");
	}catch(std::exception &e){
		std::cout << "Error: Could not load setting values from settings file\n";
		throw e;
	}
}//==================================================

/**
 *  @brief Save user settings to the prescribed file
 *  @param filename path to the XML settings file
 */
void Core_Settings::save(const std::string &filename){
	// Create empty property tree object
	using boost::property_tree::ptree;
	ptree pt;

	// Put spice data filepath in property tre
	pt.put("Astrohelion.spice.data_filepath", spice_data_filepath);
	pt.put("Astrohelion.spice.time_kernel", spice_time_kernel);
	
	try{
		// Write property tree to XML file
		write_xml(filename, pt);
	}catch(std::exception &e){
		throw e;
	}
}//==================================================


//-----------------------------------------------------
//      CLASS Core_Initializer FUNCTIONS
//-----------------------------------------------------


/**
 *  @brief Default constructor
 */
Core_Initializer::Core_Initializer() : settings(){
	std::cout << "Constructing Core_Initializer" << std::endl;
}//================================================

/**
 *	@brief Initialize the library
 */
void Core_Initializer::runInit(){
	std::cout << "Initializing Core System" << std::endl;

	// Get user home folder
	uid_t uid = getuid();
	struct passwd *pw = getpwuid(uid);

	if(pw == NULL){
		std::cout << "Error in Core_Initializer: Could not get system user database\n";
		throw std::exception();
	}

	// Define names of settings files
	const char *defaultSettingsFile = "default_settings.xml";
	const char *userSettingsFile = "user_settings.xml";

	// Construct path to directory that contains the settings file(s)
	char settingsFilepath[128], defaultSettingsFilepath[256], userSettingsFilepath[256];
	switch(OS_TYPE){
		case OS_tp::windows: 
			std::cout << "Error: Windows is not currently supported!\n";
			throw std::exception();
			break;
		case OS_tp::apple:
		case OS_tp::linux:
		default:
			sprintf(settingsFilepath, "%s%s", pw->pw_dir, "/.config/astrohelion/");
			break;
	}

	// Check to see if the directory exists; if it doesn't, create it.
	fs::path p = fs::system_complete(settingsFilepath);	// system_complete interprets any ~ or ../
	if(!fs::exists(p)){
		// Create the directory
		std::cout << "Settings directory does not exist; creating it now...\n";
		std::cout << " mkdir " << settingsFilepath << std::endl;
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

    // Get number of concurrent threads available
    unsigned int nthreads = std::thread::hardware_concurrency();
    std::cout << "Max # Concurrent Threads: " << nthreads << std::endl;

    // Tell OMP to use, at most, the maximum number of concurrent threads available
    #ifdef _OPENMP
    	omp_set_num_threads(nthreads);
    #endif
}//================================================

/**
 *	@brief Perform any duties necessary to safely shut down the library
 */
Core_Initializer::~Core_Initializer(){
	std::cout << "Destructing Core_Initializer; Closing out the Core System" << std::endl;
}

}// END of Astrohelion namespace
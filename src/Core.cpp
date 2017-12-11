/**
 *  \file Core.cpp
 *	\brief Base class for all library objects that can exist independently
 *	
 *	\author Andrew Cox
 *	\version May 25, 2016
 *	\copyright GNU GPL v3.0
 */


#include <exception> 
#include <iostream>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
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
#include "Exceptions.hpp"
#include "Utilities.hpp"

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
static char SPICE_ERR_MSG_TYPE[] = "traceback";		// other options: "short" or "none"


//-----------------------------------------------------
//      CLASS Core FUNCTIONS
//-----------------------------------------------------

/**
 *  \brief Default constructor.
 *  \details Runs an initialization sequence the first time the class is instantiated.
 */
Core::Core(){
	if(!bIsInit()){
		initializer().runInit();
		bIsInit() = true;
	}
}//======================================================

/**
 *  \brief Default destructor; doesn't do anything :)
 */
Core::~Core(){}


//-----------------------------------------------------
//      CLASS Core_Settings FUNCTIONS
//-----------------------------------------------------


/**
 *  \brief Load user settings from the prescribed file
 *  \param filename path to the XML settings file
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
		spice_data_filepath = pt.get<std::string>("astrohelion.spice.data_filepath");
		spice_time_kernel = pt.get<std::string>("astrohelion.spice.time_kernel");
		spice_spk_kernel = pt.get<std::string>("astrohelion.spice.spk_kernel");
	}catch(std::exception &e){
		std::cout << "Error: Could not load setting values from settings file\n";
		throw e;
	}
}//==================================================

/**
 *  \brief Save user settings to the prescribed file
 *  \param filename path to the XML settings file
 */
void Core_Settings::save(const std::string &filename){
	// Create empty property tree object
	using boost::property_tree::ptree;
	ptree pt;

	// Put spice data filepath in property tre
	pt.put("astrohelion.spice.data_filepath", spice_data_filepath);
	pt.put("astrohelion.spice.time_kernel", spice_time_kernel);
	pt.put("astrohelion.spice.spk_kernel", spice_spk_kernel);
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
 *  \brief Default constructor
 */
Core_Initializer::Core_Initializer() : settings(){
	// std::cout << "Constructing Core_Initializer" << std::endl;
}//================================================

/**
 *	\brief Initialize the library
 */
void Core_Initializer::runInit(){
	// std::cout << "Initializing Core System" << std::endl;

	// ************************************************************************
   	// 			Folder, Directory, and File Initialization
    // ************************************************************************
	// Get user home folder
	uid_t uid = getuid();
	struct passwd *pw = getpwuid(uid);

	if(pw == nullptr){
		std::cout << "Error in Core_Initializer: Could not get system user database\n";
		throw std::exception();
	}

	// Define names of settings files
	const char *defaultSettingsFile = "default_settings.xml";
	const char *userSettingsFile = "user_settings.xml";
	const char *bodyDataFile = "body_data.xml";

	// Construct path to directory that contains the settings file(s)
	char settingsFilepath[128], defaultSettingsFilepath[256], userSettingsFilepath[256];
	char bodyDataFilepath[256];
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
		fs::create_directories(settingsFilepath);
	}

	// Construct full paths to the settings files
	sprintf(defaultSettingsFilepath, "%s%s", settingsFilepath, defaultSettingsFile);
	sprintf(userSettingsFilepath, "%s%s", settingsFilepath, userSettingsFile);
	sprintf(bodyDataFilepath, "%s%s", settingsFilepath, bodyDataFile);

	// Save settings filepath to a variable for later
	settings.settings_filepath = std::string(settingsFilepath);
	settings.body_data_filepath = std::string(bodyDataFilepath);

	// Load user settings from file
	try{
		settings.load(userSettingsFilepath);
		// std::cout << "Loaded user settings\n";
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

    // ************************************************************************
   	// 			SPICE Initialization
    // ************************************************************************
    // Tell SPICE how to handle errors
    erract_c("set", 0, SPICE_ERR_ACTION);
    errprt_c("set", 0, SPICE_ERR_MSG_TYPE);

    try{
	    std::string spice_path = settings.spice_data_filepath;
	    std::string time_kernel = settings.spice_time_kernel;
	    std::string spk_kernel = settings.spice_spk_kernel;

	    char timeKernel[512], deKernel[512];
	    sprintf(timeKernel, "%s%s", spice_path.c_str(), time_kernel.c_str());
	    sprintf(deKernel, "%s%s", spice_path.c_str(), spk_kernel.c_str());

	    furnsh_c(timeKernel);
	    checkAndReThrowSpiceErr("Core_Initializer::runInit furnsh_c error");

	    furnsh_c(deKernel);
	    checkAndReThrowSpiceErr("Core_Initializer::runInit furnsh_c error");

	}catch(std::exception &e){}

	// ************************************************************************
   	// 			GSL Initialization
    // ************************************************************************
    // Turn GSL's error handler off; we will catch and handle the errors
    gsl_set_error_handler_off();

    // Get number of concurrent threads available
    unsigned int nthreads = std::thread::hardware_concurrency();
    // std::cout << "Max # Concurrent Threads: " << nthreads << std::endl;

    // ************************************************************************
   	// 			OpenMP Initialization
    // ************************************************************************
    // Tell OMP to use, at most, the maximum number of concurrent threads available
    #ifdef _OPENMP
    	omp_set_num_threads(nthreads);
    #endif

    // ************************************************************************
   	// 			Body Data Initialization
    // ************************************************************************
    // Create empty property tree object
	using boost::property_tree::ptree;
	ptree dataTree;

	try{
		// Load XML file and put its contents in the property tree; if file isn't found, exception is thrown
		read_xml(bodyDataFilepath, dataTree);
	}catch(std::exception &e){
		std::cout << "Could not read the body data file" << std::endl;
		throw e;
	}

	BOOST_FOREACH(ptree::value_type &v, dataTree.get_child("body_data")){
		ptree bodyTree = v.second;	// A tree for each <body> object
		Body_Data bd;

		try{
			bd.name = bodyTree.get<std::string>("name");
		}catch(std::exception &e){
			std::cout << "Error reading body data: name:\n" << e.what() << std::endl;
		}

		try{
			bd.id = bodyTree.get<int>("id");
		}catch(std::exception &e){
			std::cout << "Error reading " << bd.name << " body data: id:\n" << e.what() << std::endl;
		}

		try{
			bd.parent = bodyTree.get<std::string>("parent");
		}catch(std::exception &e){
			std::cout << "Error reading " << bd.name << " body data: parent:\n" << e.what() << std::endl;
		}

		try{
			bd.gravParam = bodyTree.get<double>("gm");
		}catch(std::exception &e){
			std::cout << "Error reading " << bd.name << " body data: gm:\n" << e.what() << std::endl;
		}

		try{
			bd.bodyRad = bodyTree.get<double>("radius");
		}catch(std::exception &e){
			std::cout << "Error reading " << bd.name << " body data: radius:\n" << e.what() << std::endl;
		}

		try{
			bd.orbitRad = bodyTree.get<double>("circ_r");
		}catch(std::exception &e){

			// Report the error unless the body is the sun, which never has a "circ_r" property
			if(!boost::iequals(bd.name, "sun"))
				std::cout << "Error reading " << bd.name << " body data: circ_r:\n" << e.what() << std::endl;
		}

		if(bd.id != 0 && allBodyData.find(bd.id) == allBodyData.end()){
			allBodyData[bd.id] = bd;
		}else{
			std::cout << "Invalid ID = " << bd.id << " and/or body already loaded" << std::endl;
		}
	}
}//================================================

/**
 *	\brief Perform any duties necessary to safely shut down the library
 */
Core_Initializer::~Core_Initializer(){
	// std::cout << "Destructing Core_Initializer; Closing out the Core System" << std::endl;

	SpiceInt count, handle;
	SpiceChar file[128], filetype[32], source[128];
	SpiceBoolean found;

	ktotal_c("all", &count);	// Get number of loaded kernels
	// std::cout << "There are " << count << " SPICE kernels loaded; unloading them..." << std::endl;
	for(SpiceInt i = 0; i < count; i++){
		kdata_c(i, "all", 128, 32, 128, file, filetype, source, &handle, &found);
		if(found){
			// std::cout << "  (" << i << ") unloading " << file << std::endl;
			unload_c(file);

			// Check for errors and report them, but catch the thrown exception
			try{
				checkAndReThrowSpiceErr("Core_Initializer::~Core_Initializer unload_c error");
			}catch(Exception &e){}
		}
	}
}//================================================

}// END of Astrohelion namespace
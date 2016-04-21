#include <iostream>
#include <pwd.h>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

int main(int argc, char *argv[]){
	(void) argc;
	(void) argv;

	uid_t uid = getuid();
	struct passwd *pw = getpwuid(uid);

	if(pw == NULL){
		std::cout << "FAILED\n";
		exit(EXIT_FAILURE);
	}

	char settingsFilepath[256];
	sprintf(settingsFilepath, "%s%s", pw->pw_dir, "/.config/tpat");

	// Check to see if the directory exists; if it doesn't, create it.
	fs::path p(settingsFilepath);
	if(fs::exists(p)){
		std::cout << "Settings directory exists!\n";
	}else{
		// Create the directory
		std::cout << "Settings directory does not exist; creating it now...\n";
		fs::create_directory(settingsFilepath);
	}
}
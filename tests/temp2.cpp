#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

int main(int argc, char **argv) {
	if(argc != 2){
		throw Exception("Expecting filepath...\n");
	}

	const char* file = argv[1];

	printf("Updating %s\n", file);
	SysData_cr3bp_lt sys(file);
	std::vector<ControlLaw*> loadedLaws {};

	Family_PO fam(&sys);
	fam.readFromMat(file, loadedLaws, true);	// reconstruct all arcs as they are loaded

	// Save data again - this will update the STMs to be the sequential guys
	fam.saveToMat(file);

	for(auto &law : loadedLaws){
		delete law;
		law = nullptr;
	}
}
#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;

int main(int argc, char** argv){
	(void) arc;
	(void) argv;

	char filename[128];

	// Load a family, propagate to y = L4, on the right-side of L4
	SysData_cr3bp_lt sys(filename);
	Family_PO_cr3bp_lt fam(&sys);

	std::vector<ControlLaw*> laws {};
	fam.readFromMat(filename, laws);

	
	return 0;
}
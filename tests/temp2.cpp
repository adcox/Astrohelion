#include "AllIncludes.hpp"

#include <exception>
#include <iostream>
#include "matio.h"
#include <vector>

using namespace astrohelion;
using ltlaw = ControlLaw_cr3bp_lt;

// Function declarations
void freeMem(std::vector<ControlLaw*>&);

// Function definitions
void freeMem(std::vector<ControlLaw*>& laws){
	for(ControlLaw* pLaw : laws){
		delete pLaw;
		pLaw = nullptr;
	}
}//====================================================

void fcn(){
	throw LinAlgException("linear algebra exception");
}
int main(){

	try{
		fcn();
	}catch(const Exception &e){
		throw(e);
	}
	
	return EXIT_SUCCESS;
}
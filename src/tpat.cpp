#include "tpat.hpp"

tpat_initializer tpat::initializer;
bool tpat::isInit;

/**
 *  @brief Default constructor.
 *  @details Runs an initialization sequence the first time the class is instantiated.
 */
tpat::tpat(){
	if(!isInit){
		initializer.runInit();
		isInit = true;
	}
}//======================================================

/**
 *  @brief Default destructor; doesn't do anything :)
 */
tpat::~tpat(){}
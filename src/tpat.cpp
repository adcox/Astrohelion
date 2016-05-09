#include "tpat.hpp"

tpat_initializer tpat::initializer;
bool tpat::isInit;

tpat::tpat(){
	if(!isInit){
		initializer.runInit();
		isInit = true;
	}
}//======================================================	
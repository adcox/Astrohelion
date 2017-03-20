/**
* Test the BodyData object
*/

#include <iostream>
#include <math.h>

#include "BodyData.hpp"

using namespace std;

int main(){
	astrohelion::BodyData sun("sun");
	astrohelion::BodyData mer("mercury");

	sun.print();
	mer.print();

	astrohelion::BodyData ear("earth");
	ear.print();
}
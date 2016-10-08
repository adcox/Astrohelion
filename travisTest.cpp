#include "Core.hpp"

#include "AsciiOutput.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"

#include <cstdio>
#include <iostream>

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

using namespace std;
using namespace astrohelion;

int main(){
	std::cout << "Compiled: " << PASS << std::endl;
}
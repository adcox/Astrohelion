/**
 *	Test the Utilities functions
 */

#include "adtk_ascii_output.hpp"
#include "adtk_utilities.hpp"

#include <iostream>

using namespace std;

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

bool testSumInt(){
	int data[] = {1,2,3,4,5};
	return adtk_sum(data, 5) == 15;
}

bool testSumDouble(){
	double data[] = {1.0, 1.1, 1.2, 1.3, 1.4};
	return adtk_sum(data, 5) == 6;
}


int main(void){

	cout << "Test sum<int>: " << (testSumInt() ? PASS : FAIL) << endl;
	cout << "Test sum<double>: " << (testSumDouble() ? PASS : FAIL) << endl;
	return 0;
}
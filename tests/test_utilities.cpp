/**
 *	Test the Utilities functions
 */

#include "tpat_ascii_output.hpp"
#include "tpat_utilities.hpp"

#include <iostream>
#include <vector>

using namespace std;

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

bool testSumInt(){
	int data[] = {1,2,3,4,5};
	return tpat_sum(data, 5) == 15;
}

bool testSumDouble(){
	double data[] = {1.0, 1.1, 1.2, 1.3, 1.4};
	return tpat_sum(data, 5) == 6;
}

bool test_concatVec(){
	vector<double> v1(2,1);
	vector<double> v2(2,2);
	vector<double> concatAns(4,1);
	concatAns[2] = 2;
	concatAns[3] = 2;

	return concatVecs(v1, v2) == concatAns;
}

bool test_concatVecInt(){
	vector<int> v1(2,1);
	vector<int> v2(2,2);
	vector<int> concatAns(4,1);
	concatAns[2] = 2;
	concatAns[3] = 2;

	return concatVecs(v1, v2) == concatAns;
}

int main(void){

	cout << "Test sum<int>: " << (testSumInt() ? PASS : FAIL) << endl;
	cout << "Test sum<double>: " << (testSumDouble() ? PASS : FAIL) << endl;
	cout << "Test concatVecs(vector<int>): " << (test_concatVecInt() ? PASS : FAIL) << endl;
	cout << "Test concatVecs(vector<double>): " << (test_concatVec() ? PASS : FAIL) << endl;
	return 0;
}
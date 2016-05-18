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
	bool check1 = tpat_util::sum(data, 5) == 15;

	std::vector<int> vecData {1,2,3,4,5};
	bool check2 = tpat_util::sum(vecData) == 15;

	return check1 && check2;
}//====================================================

bool testSumDouble(){
	double data[] = {1.0, 1.1, 1.2, 1.3, 1.4};
	bool check1 = tpat_util::sum(data, 5) == 6;

	std::vector<double> vecData {1.0, 1.1, 1.2, 1.3, 1.4};
	bool check2 = tpat_util::sum(vecData) == 6;

	return check1 && check2;
}//====================================================

bool testMeanInt(){
	int data[] = {2,4,5,1};
	bool check1 = tpat_util::mean(data,4) == 3;

	std::vector<int> vecData {2,4,5,1};
	bool check2 = tpat_util::mean(vecData) == 3;

	return check1 && check2;
}//====================================================

bool testMeanDouble(){
	double data[] = {1.1, 2.2, 3.3, 4.4, 5.5};
	bool check1 = tpat_util::mean(data, 5) == 3.3;

	std::vector<double> vecData {1.1, 2.2, 3.3, 4.4, 5.5};
	bool check2 = tpat_util::mean(vecData) == 3.3;

	return check1 && check2;
}//====================================================

bool test_concatVec(){
	vector<double> v1(2,1);
	vector<double> v2(2,2);
	vector<double> concatAns(4,1);
	concatAns[2] = 2;
	concatAns[3] = 2;

	return tpat_util::concatVecs(v1, v2) == concatAns;
}//====================================================

bool test_concatVecInt(){
	vector<int> v1(2,1);
	vector<int> v2(2,2);
	vector<int> concatAns(4,1);
	concatAns[2] = 2;
	concatAns[3] = 2;

	return tpat_util::concatVecs(v1, v2) == concatAns;
}//====================================================

bool test_permute(){
	int numSpots = 3;
	vector<int> values {1,2};
	vector<int> perms = tpat_util::generatePerms<int>(values, numSpots);
	vector<int> sol {1,1,1, 1,1,2, 1,2,1, 1,2,2, 2,1,1, 2,1,2, 2,2,1, 2,2,2};

	vector<int> values2 = {1,2,3};
	vector<int> perms2 = tpat_util::generatePerms<int>(values2);
	vector<int> sol2 = {1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1};

	// printf("permutations:\n%03d: ", 0);
	// int row = 0;
	// for(size_t i = 0; i < perms2.size(); i++){
	// 	if((i+1) % numSpots == 0){
	// 		row++;
	// 		printf("%3d\n%03d: ", perms2[i], row);
	// 	}else
	// 		printf("%3d", perms2[i]);
	// }

	return perms == sol && perms2 == sol2;
}//=================================================

bool test_getSpiceID(){
	bool test1 = getSpiceIDFromName("SUN") == 10;
	bool test2 = getSpiceIDFromName("SOLAR_SYSTEM_BARYCENTER") == 0;
	bool test3 = getSpiceIDFromName("PIONEER 12") == -12;

	return test1 && test2 && test3;
}//=================================================

bool test_getSpiceName(){
	std::string n1 = getNameFromSpiceID(199);
	std::string n2 = getNameFromSpiceID(511);
	std::string n3 = getNameFromSpiceID(-53);

	bool test1 = std::strcmp(n1.c_str(), "MERCURY") == 0;
	bool test2 = std::strcmp(n2.c_str(), "CARME") == 0;
	bool test3 = std::strcmp(n3.c_str(), "MARS SURVEYOR 01 ORBITER") == 0;

	return test1 && test2 && test3;
}//=================================================

int main(void){
	cout << "Test sum<int>: " << (testSumInt() ? PASS : FAIL) << endl;
	cout << "Test sum<double>: " << (testSumDouble() ? PASS : FAIL) << endl;
	cout << "Test mean<int>: " << (testMeanInt() ? PASS : FAIL) << endl;
	cout << "Test mean<double>: " << (testMeanDouble() ? PASS : FAIL) << endl;
	cout << "Test concatVecs(vector<int>): " << (test_concatVecInt() ? PASS : FAIL) << endl;
	cout << "Test concatVecs(vector<double>): " << (test_concatVec() ? PASS : FAIL) << endl;
	cout << "Test permutations(int): " << (test_permute() ? PASS : FAIL) << endl;
	cout << "Test SpiceIDFromName(): " << (test_getSpiceID() ? PASS : FAIL) << endl;
	cout << "Test NameFromSpiceID(): " << (test_getSpiceName() ? PASS : FAIL) << endl;
	return 0;
}//====================================================
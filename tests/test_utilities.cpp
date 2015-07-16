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
	return tpat_util::tpat_sum(data, 5) == 15;
}

bool testSumDouble(){
	double data[] = {1.0, 1.1, 1.2, 1.3, 1.4};
	return tpat_util::tpat_sum(data, 5) == 6;
}

bool test_concatVec(){
	vector<double> v1(2,1);
	vector<double> v2(2,2);
	vector<double> concatAns(4,1);
	concatAns[2] = 2;
	concatAns[3] = 2;

	return tpat_util::concatVecs(v1, v2) == concatAns;
}

bool test_concatVecInt(){
	vector<int> v1(2,1);
	vector<int> v2(2,2);
	vector<int> concatAns(4,1);
	concatAns[2] = 2;
	concatAns[3] = 2;

	return tpat_util::concatVecs(v1, v2) == concatAns;
}

bool test_permute(){
	int numSpots = 3;
	vector<int> values {1,2};
	vector<int> perms = tpat_util::generatePerms<int>(values, numSpots);
	vector<int> sol {1,1,1, 1,1,2, 1,2,1, 1,2,2, 2,1,1, 2,1,2, 2,2,1, 2,2,2};

	// vector<int> values2 = {1,2,3};
	vector<int> values2 = {0,1,2,3,4,5};
	numSpots = 6;
	vector<int> perms2 = tpat_util::generatePerms<int>(values2);
	vector<int> sol2 = {1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1};

	printf("permutations:\n%03d: ", 0);
	int row = 0;
	for(size_t i = 0; i < perms2.size(); i++){
		if((i+1) % numSpots == 0){
			row++;
			printf("%3d\n%03d: ", perms2[i], row);
		}else
			printf("%3d", perms2[i]);
	}

	return perms == sol && perms2 == sol2;
}//=================================================

int main(void){

	cout << "Test sum<int>: " << (testSumInt() ? PASS : FAIL) << endl;
	cout << "Test sum<double>: " << (testSumDouble() ? PASS : FAIL) << endl;
	cout << "Test concatVecs(vector<int>): " << (test_concatVecInt() ? PASS : FAIL) << endl;
	cout << "Test concatVecs(vector<double>): " << (test_concatVec() ? PASS : FAIL) << endl;
	cout << "Test permutations(int): " << (test_permute() ? PASS : FAIL) << endl;
	return 0;
}
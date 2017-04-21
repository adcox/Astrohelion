#define BOOST_TEST_MODULE Utilities

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <typeinfo>
#include <vector>

#include "Common.hpp"
#include "Utilities.hpp"


BOOST_AUTO_TEST_CASE(SumInt){
	int data[] = {1,2,3,4,5};
	BOOST_CHECK(astrohelion::sum(data, 5) == 15);

	std::vector<int> vecData {1,2,3,4,5};
	BOOST_CHECK(astrohelion::sum(vecData) == 15);
}//====================================================

BOOST_AUTO_TEST_CASE(SumDouble){
	double data[] = {1.0, 1.1, 1.2, 1.3, 1.4};
	BOOST_CHECK(astrohelion::sum(data, 5) == 6);

	std::vector<double> vecData {1.0, 1.1, 1.2, 1.3, 1.4};
	BOOST_CHECK(astrohelion::sum(vecData) == 6);
}//====================================================

BOOST_AUTO_TEST_CASE(MeanInt){
	int data[] = {2,4,5,1};
	BOOST_CHECK(astrohelion::mean(data,4) == 3);

	std::vector<int> vecData {2,4,5,1};
	BOOST_CHECK(astrohelion::mean(vecData) == 3);
}//====================================================

BOOST_AUTO_TEST_CASE(MeanDouble){
	double data[] = {1.1, 2.2, 3.3, 4.4, 5.5};
	BOOST_CHECK(astrohelion::mean(data, 5) == 3.3);

	std::vector<double> vecData {1.1, 2.2, 3.3, 4.4, 5.5};
	BOOST_CHECK(astrohelion::mean(vecData) == 3.3);
}//====================================================

BOOST_AUTO_TEST_CASE(ConcatVec){
	std::vector<double> v1(2,1);
	std::vector<double> v2(2,2);
	std::vector<double> concatAns(4,1);
	concatAns[2] = 2;
	concatAns[3] = 2;

	BOOST_CHECK(astrohelion::concatVecs(v1, v2) == concatAns);
}//====================================================

BOOST_AUTO_TEST_CASE(ConcatVecInt){
	std::vector<int> v1(2,1);
	std::vector<int> v2(2,2);
	std::vector<int> concatAns(4,1);
	concatAns[2] = 2;
	concatAns[3] = 2;

	BOOST_CHECK(astrohelion::concatVecs(v1, v2) == concatAns);
}//====================================================

BOOST_AUTO_TEST_CASE(Permute){
	int numSpots = 3;
	std::vector<int> values {1,2};
	std::vector<int> perms = astrohelion::generatePerms<int>(values, numSpots);
	std::vector<int> sol {1,1,1, 1,1,2, 1,2,1, 1,2,2, 2,1,1, 2,1,2, 2,2,1, 2,2,2};

	std::vector<int> values2 = {1,2,3};
	std::vector<int> perms2 = astrohelion::generatePerms<int>(values2);
	std::vector<int> sol2 = {1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1};

	// printf("permutations:\n%03d: ", 0);
	// int row = 0;
	// for(unsigned int i = 0; i < perms2.size(); i++){
	// 	if((i+1) % numSpots == 0){
	// 		row++;
	// 		printf("%3d\n%03d: ", perms2[i], row);
	// 	}else
	// 		printf("%3d", perms2[i]);
	// }

	BOOST_CHECK(perms == sol);
	BOOST_CHECK(perms2 == sol2);
}//=================================================

BOOST_AUTO_TEST_CASE(GetSpiceID){
	BOOST_CHECK(astrohelion::getSpiceIDFromName("SUN") == 10);
	BOOST_CHECK(astrohelion::getSpiceIDFromName("SOLAR_SYSTEM_BARYCENTER") == 0);
	BOOST_CHECK(astrohelion::getSpiceIDFromName("PIONEER 12") == -12);
}//=================================================

BOOST_AUTO_TEST_CASE(GetSpiceName){
	std::string n1 = astrohelion::getNameFromSpiceID(199);
	std::string n2 = astrohelion::getNameFromSpiceID(511);
	std::string n3 = astrohelion::getNameFromSpiceID(-53);

	BOOST_CHECK(std::strcmp(n1.c_str(), "MERCURY") == 0);
	BOOST_CHECK(std::strcmp(n2.c_str(), "CARME") == 0);
	BOOST_CHECK(std::strcmp(n3.c_str(), "MARS SURVEYOR 01 ORBITER") == 0);
}//=================================================

BOOST_AUTO_TEST_CASE(Sign){
	double p = -3.1415926;
	float a = 12345.12345;
	double b = 0;

	BOOST_CHECK(astrohelion::sign(p) == -1);
	BOOST_CHECK(astrohelion::sign(a) == 1);
	BOOST_CHECK(astrohelion::sign(b) == 0);
}//=================================================

BOOST_AUTO_TEST_CASE(CompareMagnitude){
	astrohelion::cdouble a(1,1);
	astrohelion::cdouble b(2,2);
	astrohelion::cdouble c(-1, -1);
	astrohelion::cdouble d(sqrt(2), 0);

	BOOST_CHECK(astrohelion::compareMagnitude(a, b) == true);
	BOOST_CHECK(astrohelion::compareMagnitude(a, c) == false);
	BOOST_CHECK(astrohelion::compareMagnitude(a, d) == false);
}//=================================================

BOOST_AUTO_TEST_CASE(To_Underlying){
	enum class a : int{ VAL1 = 1, VAL2 = 2 };
	enum class b : char{ VAL1 = 'a', VAL2 = 'z' };

	BOOST_CHECK(typeid(astrohelion::to_underlying(a::VAL1)) == typeid(int));
	BOOST_CHECK(typeid(astrohelion::to_underlying(b::VAL2)) == typeid(char));
}//=================================================






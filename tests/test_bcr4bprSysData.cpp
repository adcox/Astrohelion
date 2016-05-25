/**
 *	Test the System Data structure
 */

#include "tpat_sys_data_bc4bp.hpp"
#include "tpat_ascii_output.hpp"
#include <iostream>
#include <cstdio>

using namespace std;

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

int main(void){

	TPAT_Sys_Data_BC4BP semData("sun", "earth", "moon");

	cout << "Data for system:" << endl;
	cout << "  Type: " << semData.getTypeStr() << endl;
	cout << "  P1: " << semData.getPrimary(0) << endl;
	cout << "  P2: " << semData.getPrimary(1) << endl;
	cout << "  P3: " << semData.getPrimary(2) << endl;
	cout << "  charL: " << semData.getCharL() << " km" << endl;
	cout << "  charT: " << semData.getCharT() << " sec" << endl;
	cout << "  charM: " << semData.getCharM() << " kg" << endl;
	printf("  mu: %.5e\n", semData.getMu());
	printf("  nu: %.5e\n", semData.getNu());
	printf("  charLRatio: %.5e\n", semData.getCharLRatio());


	printf("Testing Save/Load functionality...\n");
	semData.saveToMat("semData.mat");
	TPAT_Sys_Data_BC4BP loadedData("semData.mat");

	cout << "Same CharL: " << (semData.getCharL() == loadedData.getCharL() ? PASS : FAIL) << endl;
	cout << "Same CharT: " << (semData.getCharT() == loadedData.getCharT() ? PASS : FAIL) << endl;
	cout << "Same CharM: " << (semData.getCharM() == loadedData.getCharM() ? PASS : FAIL) << endl;
	cout << "Same Mu: " << (semData.getMu() == loadedData.getMu() ? PASS : FAIL) << endl;
	cout << "Same Nu: " << (semData.getNu() == loadedData.getNu() ? PASS : FAIL) << endl;
	cout << "Same CharLRatio: " << (semData.getCharLRatio() == loadedData.getCharLRatio() ? PASS : FAIL) << endl;
	return 0;
}
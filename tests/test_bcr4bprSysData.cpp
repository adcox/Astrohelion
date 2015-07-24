/**
 *	Test the System Data structure
 */

#include "tpat_bcr4bpr_sys_data.hpp"

#include <iostream>
#include <cstdio>

using namespace std;

int main(void){

	tpat_bcr4bpr_sys_data semData("sun", "earth", "moon");

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

	return 0;
}
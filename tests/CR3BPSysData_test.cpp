/**
 *	Test the System Data structure
 */

#include "adtk_cr3bp_sys_data.hpp"

#include <iostream>
#include <cstdio>

using namespace std;

int main(void){

	adtk_cr3bp_sys_data emData("earth", "moon");

	cout << "Data for Sun-Earth System:" << endl;
	cout << "  Type: " << emData.getTypeStr() << endl;
	cout << "  P1: " << emData.getPrimary(0) << endl;
	cout << "  P2: " << emData.getPrimary(1) << endl;
	cout << "  charL: " << emData.getCharL() << " km" << endl;
	cout << "  charT: " << emData.getCharT() << " sec" << endl;
	cout << "  charM: " << emData.getCharM() << " kg" << endl;
	printf("  mu: %.5e\n", emData.getMu());

	adtk_cr3bp_sys_data seData("sun", "earth");
	cout << "Data for Sun-Earth System:" << endl;
	cout << "  Type: " << seData.getTypeStr() << endl;
	cout << "  P1: " << seData.getPrimary(0) << endl;
	cout << "  P2: " << seData.getPrimary(1) << endl;
	cout << "  charL: " << seData.getCharL() << " km" << endl;
	cout << "  charT: " << seData.getCharT() << " sec" << endl;
	cout << "  charM: " << seData.getCharM() << " kg" << endl;
	printf("  mu: %.5e\n", seData.getMu());

	return 0;
}
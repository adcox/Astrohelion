/**
 *	Test functions in the calculations library
 */

#include <iostream>

#include "tpat_ascii_output.hpp"
#include "tpat_calculations.hpp"
#include "tpat_cr3bp_sys_data.hpp"
#include "tpat_utilities.hpp"

using namespace std;

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

bool compareLPts(double *actual, double *computed, double tol){
	for(int i = 0; i < 3; i++){
		if(abs(actual[i] - computed[i]) > tol)
			return false;
	}

	return true;
}

void checkLPts(){
	printf("\nChecking Lagrange Points (Earth-Moon):\n");
	tpat_cr3bp_sys_data emSys("earth", "moon");

	double L1[3] = {0};
	double L2[3] = {0};
	double L3[3] = {0};
	double L4[3] = {0};
	double L5[3] = {0};

	double L1_true[3] = {0.836915132364302, 0, 0};
	double L2_true[3] = {1.15568216029234, 0, 0};
	double L3_true[3] = {-1.00506264525211, 0, 0};
	double L4_true[3] = {0.48784941573006, 0.866025403784439, 0};
	double L5_true[3] = {0.48784941573006, -0.866025403784439, 0};

	cr3bp_getEquilibPt(emSys, 1, 1e-14, L1);
	cr3bp_getEquilibPt(emSys, 2, 1e-14, L2);
	cr3bp_getEquilibPt(emSys, 3, 1e-14, L3);
	cr3bp_getEquilibPt(emSys, 4, 1e-14, L4);
	cr3bp_getEquilibPt(emSys, 5, 1e-14, L5);

	cout << "  L1: " << (compareLPts(L1_true, L1, 1e-14) ? PASS : FAIL) << endl;
	cout << "  L2: " << (compareLPts(L2_true, L2, 1e-14) ? PASS : FAIL) << endl;
	cout << "  L3: " << (compareLPts(L3_true, L3, 1e-14) ? PASS : FAIL) << endl;
	cout << "  L4: " << (compareLPts(L4_true, L4, 1e-14) ? PASS : FAIL) << endl;
	cout << "  L5: " << (compareLPts(L5_true, L5, 1e-14) ? PASS : FAIL) << endl;
}

void checkUDDots(){
	printf("\nChecking Pseudo-Potential DDots (Earth-Moon):\n");
	tpat_cr3bp_sys_data emSys("earth", "moon");
	double ans[] = {1.04726723737184, 0.976366381314081, -0.0236336186859191,
		1.47677525977129, 1.47677525977129, 1.47976912234663};
	double comp[6] = {0};
	cr3bp_getUDDots(emSys.getMu(), 0.5, 0.5, 0.5, comp);

	for(int i = 0; i < 6; i++){
		cout << "  U(" << i << "): " << (abs(comp[i] - ans[i]) < 1e-12 ? PASS : FAIL) << endl;
	}
}

void checkDate2EpochTime(){
	double ET_ans = 172650160;
	double ET_comp = dateToEpochTime("2005/06/21 18:21:35.816");

	cout << "Gregorian Date to Epoch Time Test: " << (floor(ET_ans) == floor(ET_comp) ? PASS : FAIL) << endl;

	double ET_comp2 = dateToEpochTime("jd 2453543.264998");

	cout << "Julian Date to Epoch Time Test: " << (floor(ET_comp) == floor(ET_comp2) ? PASS : FAIL) << endl;

}

int main(void){

	checkLPts();
	checkUDDots();
	checkDate2EpochTime();

	return 0;
}
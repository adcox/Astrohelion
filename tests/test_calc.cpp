/**
 *	Test functions in the calculations library
 */

#include <iostream>

#include "tpat_ascii_output.hpp"
#include "tpat_calculations.hpp"
#include "tpat_sys_data_cr3bp.hpp"
#include "tpat_sys_data_bcr4bpr.hpp"
#include "tpat_eigen_defs.hpp"
#include "tpat_utilities.hpp"

#include <cmath>

using namespace std;

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

bool compareLPts(double *actual, double *computed, double tol){
	for(int i = 0; i < 3; i++){
		if(abs(actual[i] - computed[i]) > tol)
			return false;
	}

	return true;
}//=========================================

void checkLPts(){
	printf("\nChecking Lagrange Points (Earth-Moon):\n");
	tpat_sys_data_cr3bp emSys("earth", "moon");

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
}//=========================================

void checkUDDots(){
	printf("\nChecking Pseudo-Potential DDots (Earth-Moon):\n");
	tpat_sys_data_cr3bp emSys("earth", "moon");
	double ans[] = {1.04726723737184, 0.976366381314081, -0.0236336186859191,
		1.47677525977129, 1.47677525977129, 1.47976912234663};
	double comp[6] = {0};
	cr3bp_getUDDots(emSys.getMu(), 0.5, 0.5, 0.5, comp);

	for(int i = 0; i < 6; i++){
		cout << "  U(" << i << "): " << (abs(comp[i] - ans[i]) < 1e-12 ? PASS : FAIL) << endl;
	}
}//=========================================

void checkDate2EpochTime(){
	double ET_ans = 172650160;
	double ET_comp = dateToEpochTime("2005/06/21 18:21:35.816");

	cout << "Gregorian Date to Epoch Time Test: " << (floor(ET_ans) == floor(ET_comp) ? PASS : FAIL) << endl;

	double ET_comp2 = dateToEpochTime("jd 2453543.264998");

	cout << "Julian Date to Epoch Time Test: " << (floor(ET_comp) == floor(ET_comp2) ? PASS : FAIL) << endl;
}//=========================================

void checkFamilyContLS(){
	double state[] = {	0.86158, 0, 0, 0,  -0.17830, 0, 2.78334, 3.16421,
  						0.86233, 0, 0, 0,  -0.18306, 0, 2.78838, 3.16299, 
  						0.86382, 0, 0, 0,  -0.19251, 0, 2.79881, 3.16050, 
  						0.86449, 0, 0, 0,  -0.19671, 0, 2.80364, 3.15936, 
  						0.86469, 0, 0, 0,  -0.19799, 0, 2.80513, 3.15901};

	std::vector<double> stateVec;
	stateVec.assign(state, state+40);

	std::vector<int> depVars;
	depVars.push_back(4);
	std::vector<double> nextData = familyCont_LS(0, 0.86469+0.0005, depVars, stateVec);

	bool check1 = std::isnan(nextData[0]);
	bool check2 = std::abs(nextData[4] + 0.20112499789815) < 1e-6;

	cout << "Least Squares Test: " << (check1 && check2 ? PASS : FAIL) << endl;

	if(!(check1 && check2)){
		printf("  Element 0 should be NAN, is actually %f\n", nextData[0]);
		printf("  Element 4 should be -0.20112499789815, is actually %.14f\n", nextData[4]);
	}
}//=========================================

void checkSPPos(){
	tpat_sys_data_bcr4bpr sys("Sun", "Earth", "Moon");
	double epoch = 100;
	Eigen::Vector3d trueSPPos(99.8234686455801, 0.000601308121756233, 0.000183265107710598);
	try{
		Eigen::Vector3d calcSPPos = bcr4bpr_getSPLoc(&sys, epoch);
		Eigen::Vector3d diff = trueSPPos - calcSPPos;
		cout << "Saddle Point Position Calc Test: " << (diff.norm() < 1e-10 ? PASS : FAIL) << endl;
	}catch(tpat_diverge &e){
		cout << "Saddle Point Position Calc test: " << FAIL << " (Diverged)" << endl;
	}
}//==========================================

int main(void){

	checkLPts();
	checkUDDots();
	checkDate2EpochTime();
	checkFamilyContLS();
	checkSPPos();

	return 0;
}//=========================================





//
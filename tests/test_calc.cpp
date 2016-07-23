/**
 *	Test functions in the calculations library
 */

#include <iostream>

#include "AsciiOutput.hpp"
#include "Calculations.hpp"
#include "SysData_cr3bp.hpp"
#include "SysData_bc4bp.hpp"
#include "EigenDefs.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"

#include <cmath>

using namespace std;
using namespace astrohelion;

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
	SysData_cr3bp emSys("earth", "moon");

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

	DynamicsModel_cr3bp::getEquilibPt(&emSys, 1, 1e-14, L1);
	DynamicsModel_cr3bp::getEquilibPt(&emSys, 2, 1e-14, L2);
	DynamicsModel_cr3bp::getEquilibPt(&emSys, 3, 1e-14, L3);
	DynamicsModel_cr3bp::getEquilibPt(&emSys, 4, 1e-14, L4);
	DynamicsModel_cr3bp::getEquilibPt(&emSys, 5, 1e-14, L5);

	cout << "  L1: " << (compareLPts(L1_true, L1, 1e-14) ? PASS : FAIL) << endl;
	cout << "  L2: " << (compareLPts(L2_true, L2, 1e-14) ? PASS : FAIL) << endl;
	cout << "  L3: " << (compareLPts(L3_true, L3, 1e-14) ? PASS : FAIL) << endl;
	cout << "  L4: " << (compareLPts(L4_true, L4, 1e-14) ? PASS : FAIL) << endl;
	cout << "  L5: " << (compareLPts(L5_true, L5, 1e-14) ? PASS : FAIL) << endl;
}//=========================================

void checkUDDots(){
	printf("\nChecking Pseudo-Potential DDots (Earth-Moon):\n");
	SysData_cr3bp emSys("earth", "moon");
	double ans[] = {1.04726723737184, 0.976366381314081, -0.0236336186859191,
		1.47677525977129, 1.47677525977129, 1.47976912234663};
	double comp[6] = {0};
	DynamicsModel_cr3bp::getUDDots(emSys.getMu(), 0.5, 0.5, 0.5, comp);

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
	SysData_bc4bp sys("Sun", "Earth", "Moon");
	double epoch = 100;
	Eigen::Vector3d trueSPPos(-0.176227312079519, 0.000601308121756233, 0.000183265107710598);
	try{
		Eigen::Vector3d calcSPPos = bcr4bpr_getSPLoc(&sys, epoch);
		Eigen::Vector3d diff = trueSPPos - calcSPPos;
		cout << "Saddle Point Position Calc Test: " << (diff.norm() < 1e-10 ? PASS : FAIL) << endl;
		if(diff.norm() >= 1e-10){
			printf("  True Position      : [%.15f, %f, %f]\n", trueSPPos(0), trueSPPos(1), trueSPPos(2));
			printf("  Calculated Position: [%.15f, %f, %f]\n", calcSPPos(0), calcSPPos(1), calcSPPos(2));
		}
	}catch(DivergeException &e){
		cout << "Saddle Point Position Calc test: " << FAIL << " (Diverged)" << endl;
	}
}//==========================================

void checkOrientAtEpoch(){
	double leaveHaloEpoch = dateToEpochTime("2016/04/18");

	SysData_bc4bp sys("Sun", "Earth", "moon");
	double T0 = (leaveHaloEpoch - SysData_bc4bp::REF_EPOCH)/sys.getCharT();
	SysData_bc4bp sys_shifted = sys;
	DynamicsModel_bc4bp::orientAtEpoch(leaveHaloEpoch, &sys_shifted);
	
	cout << "BCR4BP Orient At Epoch Test:" << endl;

	double epochs[] = {0, 10, -7};
	for(int i = 0; i < 3; i++){
		double primPos1[9], primVel1[9], primPos2[9], primVel2[9];
		
		DynamicsModel_bc4bp::getPrimaryPos(T0 + epochs[i], &sys, primPos1);
		DynamicsModel_bc4bp::getPrimaryVel(T0 + epochs[i], &sys, primVel1);
		DynamicsModel_bc4bp::getPrimaryPos(epochs[i], &sys_shifted, primPos2);
		DynamicsModel_bc4bp::getPrimaryVel(epochs[i], &sys_shifted, primVel2);

		bool samePos = true, sameVel = true;
		for(int n = 0; n < 9; n++){
			if( !astrohelion::aboutEquals(primPos1[n], primPos2[n], 1e-4) )
				samePos = false;
			if( !astrohelion::aboutEquals(primVel1[n], primVel2[n], 1e-4) )
				sameVel = false;
		}

		cout << "Epoch " << i << ": Position: " << (samePos ? PASS : FAIL) << " Velocity: " << (sameVel ? PASS : FAIL) << endl;
		if(!samePos){
			printf("P1 Pos1 = [%.4f, %.4f, %.4f], P1 Pos 2 (shifted epoch) = [%.4f, %.4f, %.4f]\n", primPos1[0],
				primPos1[1], primPos1[2], primPos2[0], primPos2[1], primPos2[2]);
			printf("P2 Pos1 = [%.4f, %.4f, %.4f], P2 Pos 2 (shifted epoch) = [%.4f, %.4f, %.4f]\n", primPos1[3],
				primPos1[4], primPos1[5], primPos2[3], primPos2[4], primPos2[5]);
			printf("P3 Pos1 = [%.4f, %.4f, %.4f], P3 Pos 2 (shifted epoch) = [%.4f, %.4f, %.4f]\n", primPos1[6],
				primPos1[7], primPos1[8], primPos2[6], primPos2[7], primPos2[8]);
		}
	}

}//==========================================

int main(void){

	checkLPts();
	checkUDDots();
	checkDate2EpochTime();
	checkFamilyContLS();
	checkSPPos();
	checkOrientAtEpoch();

	return 0;
}//=========================================





//
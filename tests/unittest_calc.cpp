/**
 *	Test functions in the calculations library
 */
#define BOOST_TEST_MODULE CalcTests

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <cmath>

#include "Calculations.hpp"
#include "Node.hpp"
#include "SimEngine.hpp"
#include "SysData_2bp.hpp"
#include "SysData_cr3bp.hpp"
#include "SysData_bc4bp.hpp"
#include "Arcset_2bp.hpp"
#include "EigenDefs.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"

#include <Eigen/Dense>

using namespace astrohelion;

BOOST_AUTO_TEST_SUITE(Calculations)

BOOST_AUTO_TEST_CASE(CR3BP_LPTs){
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

	for(unsigned int i = 0; i < 3; i++){
		BOOST_CHECK_SMALL(L1[i] - L1_true[i], 1e-14);
		BOOST_CHECK_SMALL(L2[i] - L2_true[i], 1e-14);
		BOOST_CHECK_SMALL(L3[i] - L3_true[i], 1e-14);
		BOOST_CHECK_SMALL(L4[i] - L4_true[i], 1e-14);
		BOOST_CHECK_SMALL(L5[i] - L5_true[i], 1e-14);
	}
}//=========================================

BOOST_AUTO_TEST_CASE(CR3BP_EM_UDDots){
	SysData_cr3bp emSys("earth", "moon");
	double ans[] = {1.04726723737184, 0.976366381314081, -0.0236336186859191,
		1.47677525977129, 1.47677525977129, 1.47976912234663};
	double comp[6] = {0};
	DynamicsModel_cr3bp::getUDDots(emSys.getMu(), 0.5, 0.5, 0.5, comp);

	for(int i = 0; i < 6; i++){
		BOOST_CHECK_SMALL(comp[i] - ans[i], 1e-12);
	}
}//=========================================

BOOST_AUTO_TEST_CASE(Date2EphemerisTime){
	double ET_ans = 172650160;
	

	// Gregorian date to ephemeris time test
	double ET_comp = dateToEphemerisTime("2005/06/21 18:21:35.816");
	BOOST_CHECK_EQUAL(std::floor(ET_ans), std::floor(ET_comp));

	// Julian date to ephemeris time test
	double ET_comp2 = dateToEphemerisTime("jd 2453543.264998");
	BOOST_CHECK_EQUAL(std::floor(ET_comp), std::floor(ET_comp2));
}//=========================================

// void checkFamilyContLS(){
// 	double state[] = {	0.86158, 0, 0, 0,  -0.17830, 0, 2.78334, 3.16421,
//   						0.86233, 0, 0, 0,  -0.18306, 0, 2.78838, 3.16299, 
//   						0.86382, 0, 0, 0,  -0.19251, 0, 2.79881, 3.16050, 
//   						0.86449, 0, 0, 0,  -0.19671, 0, 2.80364, 3.15936, 
//   						0.86469, 0, 0, 0,  -0.19799, 0, 2.80513, 3.15901};

// 	std::vector<double> stateVec;
// 	stateVec.assign(state, state+40);

// 	std::vector<int> depVars;
// 	depVars.push_back(4);
// 	std::vector<double> nextData = familyCont_LS(0, 0.86469+0.0005, depVars, stateVec);

// 	bool check1 = std::isnan(nextData[0]);
// 	bool check2 = std::abs(nextData[4] + 0.20112499789815) < 1e-6;

// 	cout << "Least Squares Test: " << (check1 && check2 ? PASS : FAIL) << endl;

// 	if(!(check1 && check2)){
// 		printf("  Element 0 should be NAN, is actually %f\n", nextData[0]);
// 		printf("  Element 4 should be -0.20112499789815, is actually %.14f\n", nextData[4]);
// 	}
// }//=========================================

BOOST_AUTO_TEST_CASE(SP_pos){
	SysData_bc4bp sys("Sun", "Earth", "Moon");
	double epoch = 100;
	double trueSPPos[3] = {-0.1762323562648524, 0.0005715303386115, 0.0001812483124701};
	
	Eigen::Vector3d calcSPPos;
	BOOST_REQUIRE_NO_THROW(calcSPPos = bcr4bpr_getSPLoc(&sys, epoch));

	for(unsigned int i = 0; i < calcSPPos.size(); i++){
		BOOST_CHECK_SMALL(trueSPPos[i] - calcSPPos(i), 1e-10);
		// printf("%.16f\n", calcSPPos(i));
	}


}//==========================================

// This unit test does not pass, but is not currently used... Uncomment and fix for future use!
// BOOST_AUTO_TEST_CASE(OrientAtEpoch){
// 	double leaveHaloEpoch = dateToEphemerisTime("2016/04/18");

// 	SysData_bc4bp sys("Sun", "Earth", "moon");
// 	double T0 = (leaveHaloEpoch - SysData_bc4bp::REF_EPOCH)/sys.getCharT();
// 	SysData_bc4bp sys_shifted = sys;
// 	DynamicsModel_bc4bp::orientAtEpoch(leaveHaloEpoch, &sys_shifted);
	

// 	double epochs[] = {0, 10, -7};
// 	for(int i = 0; i < 3; i++){
// 		double 	primPos1[9] = {0},
// 				primVel1[9] = {0},
// 				primPos2[9] = {0},
// 				primVel2[9] = {0};
		
// 		sys.getDynamicsModel()->getPrimPos(T0 + epochs[i], &sys, -1, primPos1);
// 		sys.getDynamicsModel()->getPrimVel(T0 + epochs[i], &sys, -1, primVel1);
// 		sys_shifted.getDynamicsModel()->getPrimPos(epochs[i], &sys_shifted, -1, primPos2);
// 		sys_shifted.getDynamicsModel()->getPrimVel(epochs[i], &sys_shifted, -1, primVel2);

// 		for(int n = 0; n < 9; n++){
// 			BOOST_CHECK_SMALL(primPos1[n] - primPos2[n], 1e-4);
// 			BOOST_CHECK_SMALL(primVel1[n] - primVel2[n], 1e-4);
// 		}

// 		// cout << "Epoch " << i << ": Position: " << (samePos ? PASS : FAIL) << " Velocity: " << (sameVel ? PASS : FAIL) << endl;
// 		// if(!samePos){
// 		// 	printf("P1 Pos1 = [%.4f, %.4f, %.4f], P1 Pos 2 (shifted epoch) = [%.4f, %.4f, %.4f]\n", primPos1[0],
// 		// 		primPos1[1], primPos1[2], primPos2[0], primPos2[1], primPos2[2]);
// 		// 	printf("P2 Pos1 = [%.4f, %.4f, %.4f], P2 Pos 2 (shifted epoch) = [%.4f, %.4f, %.4f]\n", primPos1[3],
// 		// 		primPos1[4], primPos1[5], primPos2[3], primPos2[4], primPos2[5]);
// 		// 	printf("P3 Pos1 = [%.4f, %.4f, %.4f], P3 Pos 2 (shifted epoch) = [%.4f, %.4f, %.4f]\n", primPos1[6],
// 		// 		primPos1[7], primPos1[8], primPos2[6], primPos2[7], primPos2[8]);
// 		// }
// 	}
// }//====================================================

BOOST_AUTO_TEST_CASE(InterpPointAtTime){
	SysData_2bp sys("earth");
	double mu = sys.getMu();
	double r0 = 6378 + 500;
	double v0 = sqrt(mu/r0);
	double ic[] = {r0, 0, 0, 0, v0, 0};
	double period = 2*PI*sqrt(r0*r0*r0/mu);

	Arcset_2bp traj(&sys);
	SimEngine sim;
	sim.setVerbosity(Verbosity_tp::NO_MSG);
	sim.runSim(ic, period, &traj);

	Node node = interpPointAtTime(&traj, period/2.0);
	BOOST_CHECK_SMALL(node.getEpoch() - period/2.0, 1e-10);
}//====================================================

BOOST_AUTO_TEST_CASE(BalanceMatRoutines){
	// Must be using RADIX = 2
	BOOST_REQUIRE_EQUAL(RADIX, 2);

	double data[] = {158.657,  23.9574,   142.03,  27.2438,  48.3038,  15.2807,
					-286.884, -42.9663, -256.583, -48.3038, -87.8345, -27.6321,
 					89.9068,  12.8387,  80.6714,  15.2807,  27.6321,  8.35401,
 					362.388,  54.7925,  325.521,  62.0495,  111.215,  34.6427,
					-458.039, -67.3285, -408.981,  -78.445, -139.574, -43.4001,
 					838.686,  124.921,  753.575,   142.03,  256.583,  80.6714};
 	
 	double perms_ans[] = {0.5, 1, 0.25, 1, 1, 2};
 	unsigned int low_ans = 0, hi_ans = 5;

	MatrixXRd A = Eigen::Map<MatrixXRd>(data, 6, 6);

	// Get eigenvalues of matrix without balancing it
	Eigen::EigenSolver<MatrixXRd> solver(A);
	BOOST_REQUIRE_EQUAL(solver.info(), Eigen::Success);

	Eigen::VectorXcd vals = solver.eigenvalues();
	MatrixXRcd vecs = solver.eigenvectors();

	// Copy matrix; balancing manipulates the input
	MatrixXRd A_bal = A;
	unsigned int low = 0, hi = 0;
	std::vector<double> perms {};

	// Balance the matrix
	BOOST_REQUIRE_NO_THROW(balanceMat(A_bal, low, hi, perms));

	BOOST_CHECK_EQUAL(low, low_ans);
	BOOST_CHECK_EQUAL(hi, hi_ans);
	BOOST_REQUIRE_EQUAL(perms.size(), 6);

	for(unsigned int i = 0; i < perms.size(); i++){
		BOOST_CHECK_SMALL(perms[i] - perms_ans[i], 1e-13);
	}

	// Get eigenvalues and eigenvectors of balanced matrix
	solver.compute(A_bal);
	Eigen::VectorXcd vals_bal = solver.eigenvalues();
	MatrixXRcd vecs_bal = solver.eigenvectors();

	// Check to see if the eigenvalues are the same (should be)
	BOOST_REQUIRE_EQUAL(vals_bal.size(), vals.size());
	for(unsigned int i = 0; i < vals.size(); i++){
		BOOST_CHECK_SMALL(std::real(vals_bal[i] - vals[i]), 1e-10);
		BOOST_CHECK_SMALL(std::imag(vals_bal[i] - vals[i]), 1e-10);
	}

	// Back-transform the eigenvectors
	BOOST_REQUIRE_NO_THROW(eigVec_backTrans(low, hi, perms, vecs_bal));

	// Compute two haves of the eigenvalue problem for the original matrix:
	// AV = VD, A = original matrix, V = back-transformed eigenvectors,
	// D = diagonal matrix of eigenvalues from balanced matrix,
	// M1 = AV, M2 = VD
	MatrixXRcd D = vals_bal.asDiagonal();
	MatrixXRcd M1 = A*vecs_bal, M2 = vecs_bal*D;
	BOOST_REQUIRE_EQUAL(M1.rows(), M2.rows());
	BOOST_REQUIRE_EQUAL(M1.cols(), M2.cols());

	for(unsigned int r = 0; r < M1.rows(); r++){
		for(unsigned int c = 0; c < M1.cols(); c++){
			BOOST_CHECK_SMALL(std::real(M1(r,c) - M2(r,c)), 1e-10);
			BOOST_CHECK_SMALL(std::imag(M1(r,c) - M2(r,c)), 1e-10);
		}
	}
}//====================================================

BOOST_AUTO_TEST_CASE(BalancedMatEigs){
	double data[] = {8934.746, 415.8098, 0, 2381.098, 1748.266, 0, -186.1255,
		-895.9026, -40.66346, 0, -238.5364, -176.2229, 0, 18.75314,
		0, 0, 0.7904959, 0, 0, 0.1851256, 0,
		2.935970e+04, 1366.152, 0, 7824.152, 5745.627, 0, -611.6995,
		-5578.959, -259.6677, 0, -1486.964, -1090.759, 0, 116.2318,
		0, 0, -2.049374, 0, 0, 7.850873e-01, 0,
		0, 0, 0, 0, 0, 0, 1};

	MatrixXRd A = Eigen::Map<MatrixXRd>(data, 7, 7);

	// Get eigenvalues of matrix without balancing it
	Eigen::EigenSolver<MatrixXRd> solver(A);
	BOOST_REQUIRE_EQUAL(solver.info(), Eigen::Success);

	Eigen::VectorXcd vals = solver.eigenvalues();

	// Now get eigenvalues by first balancing the matrix
	MatrixXRcd vecs_bal;
	std::vector<cdouble> vals_bal = getBalancedEigData(A, &vecs_bal);

	BOOST_REQUIRE_EQUAL(vals.size(), vals_bal.size());

	for(unsigned int i = 0; i < vals.size(); i++){
		BOOST_CHECK_SMALL(std::real(vals(i)) - std::real(vals_bal[i]), 1e-6);
		BOOST_CHECK_SMALL(std::imag(vals(i)) - std::imag(vals_bal[i]), 1e-6);
	}
}//====================================================

BOOST_AUTO_TEST_SUITE_END()





//
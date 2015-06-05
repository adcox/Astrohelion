/**
*	Test the matrix object and its operations
*/
#include "adtk_matrix.hpp"

#include "adtk_ascii_output.hpp"
#include "adtk_exceptions.hpp"
#include "adtk_utilities.hpp"

#include <cstdio>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;

static const char* PASS = BOLDGREEN "PASS" RESET;
static const char* FAIL = BOLDRED "FAIL" RESET;

static int I_r = 2, I_c = 2;
// static int B_r = 2, B_c = 2;
static int C_r = 2, C_c = 3;
static double I_data[] = {1,0,0,1};
static double B_data[] = {1,1,0,1};
static double C_data[] = {1,2,3,4,5,6};
static double C_data_trans[] = {1, 4, 2, 5, 3, 6};

bool test_constructor(adtk_matrix C){
	// Check the elements are correct
	for(int r=0; r<C_r; r++){
		for(int c=0; c<C_c; c++){
			if(C.at(r,c) != C_data[r*C_c + c])
				return false;
		}
	}

	return true;
}//===========================================

bool test_matMult(adtk_matrix I, adtk_matrix B, adtk_matrix C){
	// Test to make sure it catches size-mismatch
	try{
		C*C;
		return false;	// if C*C succeeds, there is a problem!
	}
	catch(adtk_exception& e){ /*expect this one, don't do anything*/ }
	catch(...){ throw; }

	// Test simple identity matrix multiplication
	adtk_matrix Q = I*I;
	if(Q != I)
		return false;

	// Test different 2x2 multiplication
	Q = B*B;
	double sol[] = {1,2,0,1};
	adtk_matrix Sol(B.getRows(), B.getCols(), sol);
	if(Q != Sol)
		return false;

	// Test 2x2 * 2x3 multiplication
	adtk_matrix Q2 = B*C;
	double sol2[] = {5,7,9,4,5,6};
	adtk_matrix Sol2(B.getRows(), C.getCols(), sol2);
	if(Q2 != Sol2)
		return false;

	return true;
}//===========================================

bool test_matMult_inPlace(adtk_matrix I, adtk_matrix B, adtk_matrix C){
	try{
		C*=I;
		return false;
	}
	catch(adtk_exception& e){ /*expect this one, don't do anything*/ }
	catch(...){ throw; }

	adtk_matrix tempI = I;
	if((tempI*=tempI) != I){
		return false;
	}

	adtk_matrix tempB = B;
	double sol[] = {1,2,0,1};
	adtk_matrix Sol(B.getRows(), B.getCols(), sol);
	if((tempB*=tempB) != Sol)
		return false;

	return true;
}

bool test_multScalar(adtk_matrix I){
	double sol[] = {5,0,0,5};
	adtk_matrix Sol(I_r, I_c, sol);

	adtk_matrix Q = I*5;
	if(Q != Sol)
		return false;

	return true;
}//===========================================

bool test_multScalar_inPlace(){
	double data[] = {1,0,0,1};
	double sol[] = {-1,0,0,-1};
	adtk_matrix Sol(I_r, I_c, sol);

	adtk_matrix I2(2,2,data);
	I2 *= -1;
	if(I2 != Sol)
		return false;

	return true;
}//===========================================

bool test_matAdd(adtk_matrix B){
	double sol[] = {2,2,0,2};
	adtk_matrix Sol(2,2,sol);
	adtk_matrix Q = B + B;

	if(Q != Sol)
		return false;

	return true;
}//===========================================

bool test_matSubtract(adtk_matrix C){
	adtk_matrix Zeros(2,3);
	if(C - C != Zeros)
		return false;

	return true;
}//===========================================

bool test_plusEquals(){
	double d[] = {1,2,3,4};
	adtk_matrix D = adtk_matrix(2,2,d);
	adtk_matrix I = adtk_matrix::I(2);

	D+=I;

	double sol[] = {2,2,3,5};
	adtk_matrix Sol(2,2,sol);

	return D == Sol;
}//===========================================

bool test_minusEquals(){
	double d[] = {1,2,3,4};
	adtk_matrix D = adtk_matrix(2,2,d);
	adtk_matrix I = adtk_matrix::I(2);

	D-=I;

	double sol[] = {0,2,3,3};
	adtk_matrix Sol(2,2,sol);

	return D == Sol;
}//===========================================

bool test_identity(adtk_matrix I){
	adtk_matrix I2 = adtk_matrix::I(2);
	// cout << endl;
	// I2.print();
	adtk_matrix Q = I - I2;
	// Q.print("%12.3e");

	return I2 == I;
}//===========================================

bool test_getRow(adtk_matrix C){
	double sol[] = {1,2,3};
	adtk_matrix Sol(1,3,sol);

	return Sol == C.getRow(0);
}//===========================================

bool test_getCol(adtk_matrix C){
	double sol[] = {1,4};
	adtk_matrix Sol(2,1,sol);

	return Sol == C.getCol(0);
}//===========================================

bool test_trans(adtk_matrix C){
	adtk_matrix C_trans(3,2, C_data_trans);
	return trans(C) == C_trans;
}//===========================================

bool test_norm(){
	double data[] = {1,2,3,4,5};
	adtk_matrix Q(1,5, data);

	return norm(Q) == sqrt(1 + 4 + 9 + 16 + 25);
}//===========================================

bool test_det(){
	double data[] = {1,2,3,4};
	adtk_matrix Q(2,2,data);
	bool test1 = det(Q) == -2;

	double data3[] = {1,2,0,0};
	adtk_matrix S(2,2,data3);
	bool test3 = det(S) == 0;

	return test1 && test3;
}//===========================================

bool test_diag(){
	double sol_data[] = {1, 0, 0, 0, 2, 0, 0, 0, 3};
	double diag_data[] = {1,2,3};

	adtk_matrix sol(3,3, sol_data);

	return sol == adtk_matrix::diag(diag_data, 3);
}//===========================================

bool test_cross(){
	double sol_data[] = {-3, 6, -3};
	double lhs_data[] = {1,2,3};
	double rhs_data[] = {4,5,6};

	adtk_matrix sol(1,3, sol_data);
	adtk_matrix lhs(1,3, lhs_data);
	adtk_matrix rhs(3,1, rhs_data);

	return sol == cross(lhs, rhs);
}//===========================================

int main(void){

	adtk_matrix I(2, 2, I_data);
	adtk_matrix B(2, 2, B_data);
	adtk_matrix C(2, 3, C_data);

	cout << "Test: getRows()... " << (C.getRows() == C_r ? PASS : FAIL) << endl;
	cout << "Test: getColss()... " << (C.getCols() == C_c ? PASS : FAIL) << endl;
	cout << "Test: Constructor... " << (test_constructor(C) ? PASS : FAIL) << endl;

	cout << "\nIf any of the above tests failed, all the following tests cannot be trusted!\n" <<endl;

	cout << "Test: operator A == B ... " << ( ((I == I) == 1 && (I == B) == 0) ? PASS : FAIL) << endl;
	cout << "Test: operator A != B ... " << ( ((I != I) == 0 && (I != B) == 1) ? PASS : FAIL) << endl;

	cout << "\nIf any of the above tests failed, all the following tests cannot be trusted!\n" <<endl;

	cout << "Test: operator A+B ... " << ( test_matAdd(B) ? PASS : FAIL) << endl;
	cout << "Test: operator A-B ... " << ( test_matSubtract(C) ? PASS : FAIL) << endl;
	cout << "Test: operator A*B ... " << ( test_matMult(I, B, C) ? PASS : FAIL) << endl;
	cout << "Test: operator A*q ... " << ( test_multScalar(I) ? PASS : FAIL) << endl;
	cout << "Test: operator A+=B ... " << ( test_plusEquals() ? PASS : FAIL) << endl;
	cout << "Test: operator A-=B ... " << ( test_minusEquals() ? PASS : FAIL) << endl;
	cout << "Test: operator A*=B ... " << ( test_matMult_inPlace(I, B, C) ? PASS : FAIL) << endl;
	cout << "Test: operator A*=q ... " << ( test_multScalar_inPlace() ? PASS : FAIL) << endl;

	cout << "Test: I()... " << (test_identity(I) ? PASS : FAIL) << endl;
	cout << "Test: diag()... " << ( test_diag() ? PASS : FAIL) << endl;
	cout << "Test: Get Row ... " << ( test_getRow(C) ? PASS : FAIL) << endl;
	cout << "Test: Get Column ... " << ( test_getCol(C) ? PASS : FAIL) << endl;
	cout << "Test: trans()... " << ( test_trans(C) ? PASS : FAIL) << endl;
	cout << "Test: norm()... " << ( test_norm() ? PASS : FAIL) << endl;
	cout << "Test: det()... " << ( test_det() ? PASS : FAIL) << endl;
	cout << "Test: cross()... " << ( test_cross() ? PASS : FAIL) << endl;
	return 0;
}

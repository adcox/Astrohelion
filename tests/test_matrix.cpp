/**
*	Test the matrix object and its operations
*/
#include "tpat_matrix.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_utilities.hpp"

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

bool test_constructor(tpat_matrix C){
	// Check the elements are correct
	for(int r=0; r<C_r; r++){
		for(int c=0; c<C_c; c++){
			if(C.at(r,c) != C_data[r*C_c + c])
				return false;
		}
	}

	return true;
}//===========================================

bool test_matMult(tpat_matrix I, tpat_matrix B, tpat_matrix C){
	// Test to make sure it catches size-mismatch
	try{
		C*C;
		return false;	// if C*C succeeds, there is a problem!
	}
	catch(tpat_exception& e){ /*expect this one, don't do anything*/ }
	catch(...){ throw; }

	// Test simple identity matrix multiplication
	tpat_matrix Q = I*I;
	if(Q != I)
		return false;

	// Test different 2x2 multiplication
	Q = B*B;
	double sol[] = {1,2,0,1};
	tpat_matrix Sol(B.getRows(), B.getCols(), sol);
	if(Q != Sol)
		return false;

	// Test 2x2 * 2x3 multiplication
	tpat_matrix Q2 = B*C;
	double sol2[] = {5,7,9,4,5,6};
	tpat_matrix Sol2(B.getRows(), C.getCols(), sol2);
	if(Q2 != Sol2)
		return false;

	return true;
}//===========================================

bool test_matMult_inPlace(tpat_matrix I, tpat_matrix B, tpat_matrix C){
	try{
		C*=I;
		return false;
	}
	catch(tpat_exception& e){ /*expect this one, don't do anything*/ }
	catch(...){ throw; }

	tpat_matrix tempI = I;
	if((tempI*=tempI) != I){
		return false;
	}

	tpat_matrix tempB = B;
	double sol[] = {1,2,0,1};
	tpat_matrix Sol(B.getRows(), B.getCols(), sol);
	if((tempB*=tempB) != Sol)
		return false;

	return true;
}

bool test_multScalar(tpat_matrix I){
	double sol[] = {5,0,0,5};
	tpat_matrix Sol(I_r, I_c, sol);

	tpat_matrix Q = I*5;
	if(Q != Sol)
		return false;

	return true;
}//===========================================

bool test_multScalar_inPlace(){
	double data[] = {1,0,0,1};
	double sol[] = {-1,0,0,-1};
	tpat_matrix Sol(I_r, I_c, sol);

	tpat_matrix I2(2,2,data);
	I2 *= -1;
	if(I2 != Sol)
		return false;

	return true;
}//===========================================

bool test_matAdd(tpat_matrix B){
	double sol[] = {2,2,0,2};
	tpat_matrix Sol(2,2,sol);
	tpat_matrix Q = B + B;

	if(Q != Sol)
		return false;

	return true;
}//===========================================

bool test_matSubtract(tpat_matrix C){
	tpat_matrix Zeros(2,3);
	if(C - C != Zeros)
		return false;

	return true;
}//===========================================

bool test_plusEquals(){
	double d[] = {1,2,3,4};
	tpat_matrix D = tpat_matrix(2,2,d);
	tpat_matrix I = tpat_matrix::I(2);

	D+=I;

	double sol[] = {2,2,3,5};
	tpat_matrix Sol(2,2,sol);

	return D == Sol;
}//===========================================

bool test_minusEquals(){
	double d[] = {1,2,3,4};
	tpat_matrix D = tpat_matrix(2,2,d);
	tpat_matrix I = tpat_matrix::I(2);

	D-=I;

	double sol[] = {0,2,3,3};
	tpat_matrix Sol(2,2,sol);

	return D == Sol;
}//===========================================

bool test_identity(tpat_matrix I){
	tpat_matrix I2 = tpat_matrix::I(2);
	// cout << endl;
	// I2.print();
	tpat_matrix Q = I - I2;
	// Q.print("%12.3e");

	return I2 == I;
}//===========================================

bool test_getRow(tpat_matrix C){
	double sol[] = {1,2,3};
	tpat_matrix Sol(1,3,sol);

	return Sol == C.getRow(0);
}//===========================================

bool test_getCol(tpat_matrix C){
	double sol[] = {1,4};
	tpat_matrix Sol(2,1,sol);

	return Sol == C.getCol(0);
}//===========================================

bool test_trans(tpat_matrix C){
	tpat_matrix C_trans(3,2, C_data_trans);
	return trans(C) == C_trans;
}//===========================================

bool test_norm(){
	double data[] = {1,2,3,4,5};
	tpat_matrix Q(1,5, data);

	return norm(Q) == sqrt(1 + 4 + 9 + 16 + 25);
}//===========================================

bool test_det(){
	double data[] = {1,2,3,4};
	tpat_matrix Q(2,2,data);
	bool test1 = det(Q) == -2;

	double data3[] = {1,2,0,0};
	tpat_matrix S(2,2,data3);
	bool test3 = det(S) == 0;

	return test1 && test3;
}//===========================================

bool test_diag(){
	double sol_data[] = {1, 0, 0, 0, 2, 0, 0, 0, 3};
	double diag_data[] = {1,2,3};

	tpat_matrix sol(3,3, sol_data);

	return sol == tpat_matrix::diag(diag_data, 3);
}//===========================================

bool test_cross(){
	double sol_data[] = {-3, 6, -3};
	double lhs_data[] = {1,2,3};
	double rhs_data[] = {4,5,6};

	tpat_matrix sol(1,3, sol_data);
	tpat_matrix lhs(1,3, lhs_data);
	tpat_matrix rhs(3,1, rhs_data);

	return sol == cross(lhs, rhs);
}//===========================================

bool test_eig(){
	double D_data[] = {1,2,3,0,4,5,0,0,6};
	tpat_matrix D(3,3,D_data);

	std::vector< std::vector<cdouble> > eigData = eig(D);
	std::vector<cdouble> vals = eigData[0];
	std::vector<cdouble> vecs = eigData[1];

	if(vals.size() != 3)
		return false;

	if(vecs.size() != 9)
		return false;

	bool check1 = real(vals[0]) == 1 && imag(vals[0]) == 0;
	bool check2 = real(vals[1]) == 4 && imag(vals[1]) == 0;
	bool check3 = real(vals[2]) == 6 && imag(vals[2]) == 0;

	bool check_vec1 = real(vecs[0]) == 1 && imag(vecs[0]) == 0;
	check_vec1 = check_vec1 && real(vecs[1]) == 0 && imag(vecs[1]) == 0;
	check_vec1 = check_vec1 && real(vecs[2]) == 0 && imag(vecs[2]) == 0;

	bool check_vec2 = tpat_util::aboutEquals(real(vecs[3]), 0.554700196225229, 1e-13) && imag(vecs[3]) == 0;
	check_vec2 = check_vec2 && tpat_util::aboutEquals(real(vecs[4]), 0.832050294337844, 1e-13) && imag(vecs[4]) == 0;
	check_vec2 = check_vec2 && real(vecs[5]) == 0 && imag(vecs[5]) == 0;

	bool check_vec3 = tpat_util::aboutEquals(real(vecs[6]), 0.510840685451281, 1e-13) && imag(vecs[6]) == 0;
	check_vec3 = check_vec3 && tpat_util::aboutEquals(real(vecs[7]), 0.798188571017626, 1e-13) && imag(vecs[7]) == 0;
	check_vec3 = check_vec3 && tpat_util::aboutEquals(real(vecs[8]), 0.31927542840705, 1e-13) && imag(vecs[8]) == 0;

	bool matrix1Pass = check1 && check2 && check3 && check_vec1 && check_vec2 && check_vec3;

	// Now check a matrix with complex eigenvalues and eigenvectors
	double E_data[] = {1,2,3,4,5,6,1,5,8};
	tpat_matrix E(3,3,E_data);

	vals.clear();
	vecs.clear();
	eigData.clear();
	eigData = eig(E);
	vals = eigData[0];
	vecs = eigData[1];

	if(vals.size() != 3)
		return false;

	if(vecs.size() != 9)
		return false;

	bool check4 = tpat_util::aboutEquals(real(vals[0]), 13.1015555756727, 1e-13) && imag(vals[0]) == 0;
	bool check5 = tpat_util::aboutEquals(real(vals[1]), 0.44922221216366, 1e-13) && 
		tpat_util::aboutEquals(imag(vals[1]), 0.164863116315263, 1e-13);
	bool check6 = real(vals[2]) == real(vals[1]) && imag(vals[2]) == -1*imag(vals[1]);

	bool check_vec4 = tpat_util::aboutEquals(real(vecs[0]), 0.281912833008201, 1e-13) && imag(vecs[0]) == 0;
	check_vec4 = check_vec4 && tpat_util::aboutEquals(real(vecs[1]), 0.657008011457982, 1e-13) && imag(vecs[1]) == 0;
	check_vec4 = check_vec4 && tpat_util::aboutEquals(real(vecs[2]), 0.699189264409372, 1e-13) && imag(vecs[2]) == 0;

	// Matlab and GSL vary a bit here - the vectors are very similar, but not exact.
	bool check_vec5 = tpat_util::aboutEquals(real(vecs[3]), 0.150692642039344, 1e-2) &&
		tpat_util::aboutEquals(imag(vecs[3]), -0.0644436534652211, 1e-2);
	check_vec5 = check_vec5 && tpat_util::aboutEquals(real(vecs[4]), -0.83159310577757, 1e-2) && 
		tpat_util::aboutEquals(imag(vecs[4]), 0.0, 1e-2);
	check_vec5 = check_vec5 && tpat_util::aboutEquals(real(vecs[5]), 0.530270811022171, 1e-2) && 
		tpat_util::aboutEquals(imag(vecs[5]), 0.0201125971560176, 1e-2);

	bool check_vec6 = real(vecs[3]) == real(vecs[6]) && imag(vecs[3]) == -1*imag(vecs[6]);
	check_vec6 = check_vec6 && real(vecs[4]) == real(vecs[7]) && imag(vecs[4]) == -1*imag(vecs[7]);
	check_vec6 = check_vec6 && real(vecs[5]) == real(vecs[8]) && imag(vecs[5]) == -1*imag(vecs[8]);
	
	bool matrix2Pass = check4 && check5 && check6 && check_vec4 && check_vec5 && check_vec6;

	return matrix1Pass && matrix2Pass;
}//===========================================

int main(void){

	tpat_matrix I(2, 2, I_data);
	tpat_matrix B(2, 2, B_data);
	tpat_matrix C(2, 3, C_data);

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
	cout << "Test: eig()... " << (test_eig() ? PASS : FAIL) << endl;
	return 0;
}

/**
 * 	Matrix object that employs GSL's Matrix and CBLAS functionality
 *
 *	The matrix is double precision, rectangular, and real. No other
 *	characteristics are assumed (e.g. banded, symmetric, etc.)
 *
 *	Author: Andrew Cox
 *	Version: May 6, 2015
 */

#include <exception>
#include <iostream>

#include <gsl/gsl_cblas.h>

#include "adtk_matrix.hpp"

using namespace std;

/**
 * 	Exception class for matrices: size mismatch means two matrices cannot
 *	be added, subtracted, or multiplied because they do not have the
 *	appropriate dimensions
 */
class sizeMismatch: public std::exception{
  virtual const char* what() const throw()
  {
    return "Matrix dimensions do not match!";
  }
} adtk_matSizeMismatch;

//---------------------------------------------------------------------------
//    	Constructors and Destructor
//---------------------------------------------------------------------------

/**
 *	Construct a matrix full of zeros
 *	@param r number of rows
 *	@param c number of columns
 */
adtk_matrix::adtk_matrix(const int r, const int c){
	initBasicMatrix(r, c);
}

/**
 *	Construct a matrix from an array of data
 *	@param r number of rows
 *	@param c number of columns
 *	@param data a pointer to an array of data with r*c elements; the array holds the data by row,
 *	e.g. if we create an nx3 matrix, the first 3 elements of data are the first row, the next three
 *	elements are the second row, etc.
 */
adtk_matrix::adtk_matrix(int r, int c, double data[]){
	initBasicMatrix(r, c);
	for (int i = 0; i < rows; i++){
	    for (int j = 0; j < cols; j++){
	    	gsl_matrix_set(a, i, j, data[i*cols+j]);
	    }
	}
}//==============================================

/**
 *	Construct a matrix from an array of data
 *	@param r number of rows
 *	@param c number of columns
 *	@param data a vector containint data with r*c elements; the vector holds the data by row,
 *	e.g. if we create an nx3 matrix, the first 3 elements of data are the first row, the next three
 *	elements are the second row, etc.
 */
adtk_matrix::adtk_matrix(int r, int c, std::vector<double> data){
	initBasicMatrix(r, c);
	for (int i = 0; i < rows; i++){
	    for (int j = 0; j < cols; j++){
	    	gsl_matrix_set(a, i, j, data[i*cols+j]);
	    }
	}
}

/**
 *	Destructor: called when the matrix is deleted or goes out of scope
 *	Primary concern is to de-allocate the memory taken up by the gsl_matrix
 */
adtk_matrix::~adtk_matrix(){
	try{
		gsl_matrix_free(a);
	}catch(...){
		cout << "\nError! cannot free a: my dimensions are " << rows << " x " << cols << endl;
	}
}//==============================================

/**
 *	Copy constructor: creates a new matrix from this one without using the same memory!
 */
adtk_matrix::adtk_matrix(const adtk_matrix &b){
	initBasicMatrix(b.rows, b.cols);
	for(int r=0; r<rows; r++){
		for(int c=0; c<cols; c++){
			gsl_matrix_set(a, r, c, gsl_matrix_get(b.a, r, c));
		}
	}
}//=================================================

/**
 * 	Create an identity matrix
 *	@param size the number of rows or columns the (square) matrix has
 *	@return an identiy matrix
 */
adtk_matrix adtk_matrix::Identity(int size){
	vector<double> data(size*size);
	for(int i = 0; i < size*size; i += size+1){
		data[i] = 1;
	}
	return adtk_matrix(size, size, data);
}//==============================================

//---------------------------------------------------------------------------
//    	Set and Get
//---------------------------------------------------------------------------

/**
 * 	@return the number of rows in this matrix
 */
int adtk_matrix::getRows(){return rows;}

/**
 * 	@return the number of columns in this matrix
 */
int adtk_matrix::getCols(){return cols;}

/**
 *	@return the total number of elements in the matrix
 */
int adtk_matrix::getLength(){return rows*cols;}

/**
 *	@return a pointer to the gsl_matrix object that represents this matrix
 */
gsl_matrix* adtk_matrix::getGSLMat(){return a;}

/**
 *	@param r the row
 *	@param c the column
 *	@return the value in the specified row and column
 */
double adtk_matrix::get(int r, int c){return gsl_matrix_get(a, r, c);}

//---------------------------------------------------------------------------
//   	Overloaded Operators
//---------------------------------------------------------------------------

/**
 *	Copy Operator; copy a matrix into this one
 *	@param b a constant matrix
 *	@return this matrix
 */
adtk_matrix& adtk_matrix::operator =(const adtk_matrix &b){
	rows = b.rows;
	cols = b.cols;

	for(int r = 0; r<rows; r++){
		for(int c = 0; c<cols; c++){
			gsl_matrix_set(a, r, c, gsl_matrix_get(b.a, r, c));
		}
	}
	return *this;
}//================================================

/**
 *	Add two matrices; they must have the same size.
 *	@param b a matrix to add to this one
 *	@return a new matrix that is the sum of this one and b
 */
adtk_matrix adtk_matrix::operator +(const adtk_matrix &b){
	if(rows == b.rows && cols == b.cols){
		vector<double> q(rows*cols);
		for(int r = 0; r<rows; r++){
			for(int c = 0; c<cols; c++){
				q[r*cols+c] = gsl_matrix_get(a, r, c) + gsl_matrix_get(b.a, r, c);
			}
		}
		// Return a new matrix (do not modify this one)
		return adtk_matrix(rows, cols, q);
	}else{
		cout << "Matrices must be the same size to apply addition!" << endl;
		throw adtk_matSizeMismatch;
	}
}//==============================================

/**
 *	Add a matrix to this one; perform operation in place
 *	@param b a matrix, must be the same size as this matrix
 *	@return this matrix
 */
adtk_matrix& adtk_matrix::operator +=(const adtk_matrix &b){
	if(rows == b.rows && cols == b.cols){
		gsl_matrix_add(a, b.a);
		return *this;
	}else{
		cout << "Matrices must be the same size to apply addition!" << endl;
		throw adtk_matSizeMismatch;
	}
}//==============================================

/**
 *	Subtract a matrix from this one; both matrices must be the same size
 *	@param b a matrix to subtract from this one
 *	@return a new matrix
 */
adtk_matrix adtk_matrix::operator -(const adtk_matrix &b){
	if(rows == b.rows && cols == b.cols){
		vector<double> q(rows*cols);
		for(int r = 0; r<rows; r++){
			for(int c = 0; c<cols; c++){
				q[r*cols+c] = gsl_matrix_get(a, r, c) - gsl_matrix_get(b.a, r, c);
			}
		}
		return adtk_matrix(rows, cols, q);
	}else{
		cout << "Matrices must be the same size to apply subtraction!" << endl;
		throw adtk_matSizeMismatch;
	}
}//==============================================

/**
 *	In place subtraction; subtract a matrix from this one
 *	@param b a matrix to subtract from this
 *	@return this matrix
 */
adtk_matrix& adtk_matrix::operator -=(const adtk_matrix &b){
	if(rows == b.rows && cols == b.cols){
		gsl_matrix_sub(a, b.a);
		return *this;
	}else{
		cout << "Matrices must be the same size to apply addition!" << endl;
		throw adtk_matSizeMismatch;
	}
}//==============================================

/**
 *	Multiply this matrix by another.
 *	@param b a matrix
 *	
 *	Let this matrix be A, and the input matrix be B; A is m x k, B is k x n;
 *
 *	@return C = A*B where C is an m x n matrix
 */
adtk_matrix adtk_matrix::operator *(const adtk_matrix &b){
	if(cols != b.rows){
		throw adtk_matSizeMismatch;
	}

	// Get the arrays of doubles from both matrix objects
	double *A = a->data;
	double *B = b.a->data;

	// Create a new array that will hold the data from the mulitplication
	double *C = new double[rows * b.cols];

	/*
	This function multiplies two double-precision matrices: C <- alpha*A*B + beta*C
	  gsl_blas_dgemm(type, transpose A?, transpose B?, m, n, k, alpha, A, k, B, n, beta, C, n)
		- CblasRowMajor indicates the matrices are stored in row major order
		- CblasNoTrans (in both spots) tells the function not to transpose A or B
		- A: array representing m x k matrix
		- B: array representing k x n matrix
		- C: array representing m x n matrix
		- alpha, beta: real scalars
	*/
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, b.cols, cols,
		1.0, A, cols, B, b.cols, 0.0, C, b.cols);

	// Create a new matrix object from the array we just computed
	adtk_matrix ans(rows, b.cols, C);
	delete[] C;
	return ans;
}//==============================================

/**
 * 	Multiplies this matrix, A ( m x m ) with an input matrix B (m x m) and sets A = AB
 *
 *	This operation only works for square matrices (non-square doesn't make any sense,
 *	because you can't assign the new product back to the old matrix)
 *
 *	@param b a square matrix the same size as this one
 *	@return this matrix, modified
 */
adtk_matrix& adtk_matrix::operator *=(const adtk_matrix &b){

	if(cols == rows && b.rows == b.cols && cols == b.cols){
		// Get the arrays of doubles from both matrix objects
		double *A = a->data;
		double *B = b.a->data;
		double *C = new double[rows * b.cols];

		// Do the multiplication, 
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, cols, cols,
			1.0, A, cols, B, cols, 0.0, C, cols);

		// Store data back in this matrix
		for (int i = 0; i < rows; i++){
		    for (int j = 0; j < cols; j++){
		    	gsl_matrix_set(a, i, j, C[i*cols+j]);
		    }
		}

		delete[] C;
		return *this;
	}else{
		throw adtk_matSizeMismatch;
	}
}//=============================================

/**
 *	Multiply the matrix by a scalar and return the product (does not modify this matrix)
 *	@param alpha a real scalar
 *	@return the matrix alpha*A
 */
adtk_matrix adtk_matrix::operator *(const double &alpha){
	vector<double> b(rows*cols);
	for(int r = 0; r<rows; r++){
		for(int c = 0; c<rows; c++){
			b[r*cols + c] = alpha*gsl_matrix_get(a, r, c);
		}
	}

	return adtk_matrix(rows, cols, b);
}//==============================================

/**
 *	Multiply this matrix by a scalar (in place, DOES modify this matrix)
 * 	@param alpha a real scalar
 */
adtk_matrix& adtk_matrix::operator *=(const double &alpha){
	gsl_matrix_scale(a, alpha);
	return *this;
}//==============================================

/**
 *	Compare this matrix to another; note that only EXACT equivalence will
 *	return true.
 *
 *	@param B a matrix
 *	@return whether or not this matrix contains the same elements as B
 */
bool adtk_matrix::operator ==(const adtk_matrix &B){
	return gsl_matrix_equal(a, B.a);
}//==============================================

/**
 *	Compare this matrix to another
 *
 *	@param B a matrix
 *	@return whether or not this matrix is NOT equal to B
 */
bool adtk_matrix::operator !=(const adtk_matrix &B){
	return !gsl_matrix_equal(a, B.a);
}//==============================================

//---------------------------------------------------------------------------
//    	Matrix Operation Functions
//---------------------------------------------------------------------------


/**
 *	Transpose this matrix (does not change the original copy of the matrix)
 *	@return the transpose (NOT the complex conjugate transpose)
 */
adtk_matrix adtk_matrix::trans(){
	gsl_matrix *b = gsl_matrix_alloc(cols, rows);
	gsl_matrix_transpose_memcpy(b, a);
	return adtk_matrix(cols, rows, b->data);
}//==============================================


//---------------------------------------------------------------------------
//    	Utilities
//---------------------------------------------------------------------------

/**
 *	Initialize basic matrix and set all values to zero
 *	@param r number of rows
 *	@param c number of columns
 */
void adtk_matrix::initBasicMatrix(int r, int c){
	a = gsl_matrix_alloc (r, c);
	rows = r;
	cols = c;
	gsl_matrix_set_zero(a);
}//============================================

/**
 * 	Print the matrix to standard output via printf
 */
void adtk_matrix::print(){
	print("%8.4f");		// Use the more advanced function
}//============================================

/**
 *	Print the matrix to standard output with the specified format
 *	@param format a format string for the std::printf function
 */
void adtk_matrix::print(const char *format){
	for (int r = 0; r < rows; r++){
	    for (int c = 0; c < cols; c++){
	    	printf(format, gsl_matrix_get(a, r, c));
	    }
	    printf("\n");
	}
}//==============================================
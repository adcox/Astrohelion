/**
 *	@file adtk_matrix.cpp
 *
 * 	Matrix object that employs GSL's Matrix and CBLAS functionality
 *
 *	The matrix is double precision, rectangular, and real. No other
 *	characteristics are assumed (e.g. banded, symmetric, etc.)
 *
 *	@author Andrew Cox
 *	@version May 6, 2015
 */
#include "adtk_matrix.hpp"

#include <cmath>
#include <exception>
#include <fstream>
#include <gsl/gsl_cblas.h>
#include <iostream>
#include <stdexcept>
 
using namespace std;

/**
 *	@brief Exception class for matrices, fired when the dimension of matrices
 *	makes a function impossible.
 *
 * 	Size mismatch means two matrices cannot
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
 *	@brief Construct a matrix full of zeros
 *	@param r number of rows
 *	@param c number of columns
 */
adtk_matrix::adtk_matrix(const int r, const int c){
	initBasicMatrix(r, c);
}

/**
 *	@brief Construct a matrix from an array of data
 *	@param r number of rows
 *	@param c number of columns
 *	@param data a pointer to an array of data with r*c elements; the array holds the data by row,
 *	e.g. if we create an nx3 matrix, the first 3 elements of data are the first row, the next three
 *	elements are the second row, etc.
 */
adtk_matrix::adtk_matrix(int r, int c, double *data){
	initBasicMatrix(r, c);
	copyDataIntoGSL_Matrix(data);
}//==============================================

/**
 *	@brief Construct a matrix from an array of data
 *	@param r number of rows
 *	@param c number of columns
 *	@param data a vector containint data with r*c elements; the vector holds the data by row,
 *	e.g. if we create an nx3 matrix, the first 3 elements of data are the first row, the next three
 *	elements are the second row, etc.
 */
adtk_matrix::adtk_matrix(int r, int c, std::vector<double> data){
	initBasicMatrix(r, c);
	copyDataIntoGSL_Matrix(data);
}//==================================

/**
 *	@brief Construct a matrix from a <tt>gsl_matrix</tt> object
 *	@param m a pointer to a <tt>gsl_matrix</tt>
 */
adtk_matrix::adtk_matrix(gsl_matrix *m){
	// Copy data from gsl_matrix pointer
	rows = m->size1;
	cols = m->size2;
	a = m;
}//===================================

/**
 *	@brief Construct a 1-D matrix (or vector) from a <tt>gsl_vector</tt>
 *	@param v a pointer to a <tt>gsl_vector</tt>
 *	@param isRow whether or not the matrix is a row-vector. If not, it is a column vector
 */
adtk_matrix::adtk_matrix(gsl_vector *v, bool isRow){
	if(isRow){
		initBasicMatrix(1, v->size);
	}else{
		initBasicMatrix(v->size, 1);
	}

	copyDataIntoGSL_Matrix(v->data);
}//====================================

/**
 *	@brief Destructor: called when the matrix is deleted or goes out of scope
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
 *	@brief Copy constructor: creates a new matrix from this one without using the same memory!
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
 * 	@brief Create an identity matrix
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
int adtk_matrix::getRows() const {return rows;}

/**
 *	@brief Get a row from this matrix
 *	@param i the row index (begins with zero)
 *	@return a 1-D matrix (row vector) that contains the data from row <tt>i</tt>
 */
adtk_matrix adtk_matrix::getRow(int i){
	// Retrieve the row data
	gsl_vector *v = gsl_vector_alloc(static_cast<size_t>(cols));
	gsl_matrix_get_row(v, a, static_cast<size_t>(i));
	adtk_matrix temp(v, true);
	gsl_vector_free(v);
	return temp;
}

/**
 * 	@return the number of columns in this matrix
 */
int adtk_matrix::getCols() const {return cols;}

/**
 *	@brief Get a column from this matrix
 *	@param j the column index (begins with zero)
 *	@return a 1-D matrix (column vector) that contains the data from column <tt>j</tt>
 */
adtk_matrix adtk_matrix::getCol(int j){
	gsl_vector *v = gsl_vector_alloc(static_cast<size_t>(rows));
	gsl_matrix_get_col(v, a, static_cast<size_t>(j));
	adtk_matrix temp(v, false);
	gsl_vector_free(v);
	return temp;
}

/**
 *	@return the total number of elements in the matrix
 */
int adtk_matrix::getLength() const {return rows*cols;}

/**
 *	@return a pointer to the gsl_matrix object that represents this matrix
 */
gsl_matrix* adtk_matrix::getGSLMat(){return a;}

/**
 *	@return a pointer to the array of doubles that represents this matrix;
 *	values are stored in row-major order
 */
double* adtk_matrix::getDataPtr(){ return a->data; }

/**
 *	@param r the row
 *	@param c the column
 *	@return the value in the specified row and column
 */
double adtk_matrix::at(int r, int c) const {
	if(r >= 0 && r < rows && c >= 0 && c < cols){
		return gsl_matrix_get(a, r, c);
	}else{
		throw std::out_of_range("Index out of range");
	}
}//=================================

/**
 *	@brief Retrive a value from a 1D matrix. If the matrix is not 1D, a value from the first row
 *	is returned.
 */
double adtk_matrix::at(int i) const {
	if(i < rows*cols){
		if(rows != 1 && cols != 1){
			cout << "adtk_matrix :: Warning - matrix is not 1D; returning element from storage array" << endl;
		}
		return a->data[i];
	}else{
		throw std::out_of_range("Index out of range");
	}
}//===============================

//---------------------------------------------------------------------------
//   	Overloaded Operators
//---------------------------------------------------------------------------

/**
 *	@brief Copy Operator; copy a matrix into this one
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
 *	@brief Add two matrices; they must have the same size.
 *	@return a new matrix that is the sum of lhs and rhs
 */
adtk_matrix operator +(const adtk_matrix &lhs, const adtk_matrix &rhs){
	if(lhs.rows == rhs.rows && lhs.cols == rhs.cols){
		vector<double> q(lhs.rows * lhs.cols);
		for (int r = 0; r < lhs.rows; r++){
			for (int c = 0; c < lhs.cols; c++){
				q[r * lhs.cols + c] = gsl_matrix_get(lhs.a, r, c) + gsl_matrix_get(rhs.a, r, c);
			}
		}
		// Return a new matrix
		return adtk_matrix(lhs.rows, lhs.cols, q);
	}else{
		cout << "Matrices must be the same size to apply addition!" << endl;
		throw adtk_matSizeMismatch;
	}
}

/**
 *	@brief Add a matrix to this one; perform operation in place
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
 *	@brief Subtract two matrices: A = L - R
 *	@param lhs the left-hand matrix L
 *	@param rhs the right-hand matrix R
 *	@return the difference L - R
 */
adtk_matrix operator -(const adtk_matrix &lhs, const adtk_matrix &rhs){
	return operator +(lhs, -1*rhs);
}//==============================================

/**
 *	@brief In place subtraction; subtract a matrix from this one
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
 *	@brief Multiply this matrix by another.
 *	@param b a matrix
 *	
 *	Let this matrix be A, and the input matrix be B; A is m x k, B is k x n;
 *
 *	@return C = A*B where C is an m x n matrix
 */
adtk_matrix adtk_matrix::operator *(const adtk_matrix &b) const{
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
 * 	@brief Multiplies this matrix, A ( m x m ) with an input matrix B (m x m) and sets A = AB
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
 *	@brief Multiply a scalar and a matrix: A = l*R
 *	@param lhs the scalar l
 *	@param rhs the matrix R
 *	@return the product l*R (same as R*l)
 */
adtk_matrix operator *(const double &lhs, const adtk_matrix &rhs){
	int rows = rhs.rows;
	int cols = rhs.cols;
	vector<double> b(rows * cols);
	for(int r = 0; r< rows; r++){
		for(int c = 0; c< cols; c++){
			b[r * cols + c] = lhs *gsl_matrix_get(rhs.a, r, c);
		}
	}

	return adtk_matrix(rows, cols, b);
}//=============================================

/**
 *	@brief Multiply a scalar and a matrix: A = L*r
 *	@param lhs the matrix L
 *	@param rhs the scalar r
 *	@return the product L*r (same as r*L)
 */
adtk_matrix operator *(const adtk_matrix &lhs, const double &rhs){
	return operator *(rhs, lhs);
}//=============================================

/**
 *	@brief Multiply a scalar and a matrix: A = L/r
 *	@param lhs the matrix L
 *	@param rhs the scalar r
 *	@return the quotient L/r (same as 1/r * L)
 */
adtk_matrix operator /(const adtk_matrix &lhs, const double &rhs){
	return operator *(lhs, 1/rhs);
}

/**
 *	@brief Multiply this matrix by a scalar (in place, DOES modify this matrix)
 * 	@param alpha a real scalar
 */
adtk_matrix& adtk_matrix::operator *=(const double &alpha){
	gsl_matrix_scale(a, alpha);
	return *this;
}//==============================================

/**
 *	@brief Compare two matrices
 *	@param lhs the left-hand matrix L
 *	@param rhs the right-hand matrix R
 *	@return whether or not L and R contain the same values
 */
bool operator ==(const adtk_matrix &lhs, const adtk_matrix &rhs){
	return gsl_matrix_equal(lhs.a, rhs.a);
}

/**
 *	@brief Compare two matrices
 *	@param lhs the left-hand matrix L
 *	@param rhs the right-hand matrix R
 *	@return whether or not L and R do NOT contain the same values
 */
bool operator !=(const adtk_matrix &lhs, const adtk_matrix &rhs){
	return !operator==(lhs, rhs);
}

//---------------------------------------------------------------------------
//    	Matrix Operation Functions
//---------------------------------------------------------------------------


/**
 *	@brief Transpose this matrix (does not change the original copy of the matrix)
 *	@return the transpose (NOT the complex conjugate transpose)
 */
adtk_matrix adtk_matrix::trans() const{
	gsl_matrix *b = gsl_matrix_alloc(cols, rows);
	gsl_matrix_transpose_memcpy(b, a);
	adtk_matrix temp(cols, rows, b->data);
	gsl_matrix_free(b);
	return temp;
}//==============================================

/** 
 *	@brief Get the norm of this "matrix;" only applies for 1-D matrices (vectors)
 *	@return the 2-norm (sqrt of sum of squared elements) of the matrix
 */
double adtk_matrix::norm() const{
	if(rows == 1 || cols == 1){
		double sumSquares = 0;
		for(int i = 0; i < rows*cols; i++){
			sumSquares += (a->data[i])*(a->data[i]);
		}
		return sqrt(sumSquares);
	}else{
		cout << "adtk_matrix :: Cannot take norm of a 2-D matrix" << endl;
		throw;
	}
}//==========================================

//---------------------------------------------------------------------------
//    	Utilities
//---------------------------------------------------------------------------

/**
 *	@brief Initialize basic matrix and set all values to zero
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
 *	@brief Utility to copy data from an array of doubles into the gsl_matrix object that
 *	contains the data for this matrix. 
 *
 *	NOTE: the <tt>rows</tt> and <tt>cols</tt>
 *	variabls MUST be initialized and set to the correct values before running
 *	this function.
 *	@param data an array (row-major order) containing data for this matrix
 */
void adtk_matrix::copyDataIntoGSL_Matrix(double *data){
	for (int i = 0; i < rows; i++){
	    for (int j = 0; j < cols; j++){
	    	gsl_matrix_set(a, i, j, data[i*cols+j]);
	    }
	}
}//============================================

/**
 * 	@brief Utility to copy data from a vector into the gsl_matrix object that
 *	contains the data for this matrix. 
 *
 *	NOTE: the <tt>rows</tt> and <tt>cols</tt>
 *	variabls MUST be initialized and set to the correct values before running
 *	this function.
 *	@param v a vector (row-major order) containing data for this matrix
 */
void adtk_matrix::copyDataIntoGSL_Matrix(std::vector<double> v){
	copyDataIntoGSL_Matrix(&(v[0]));
}//============================================

/**
 * 	@brief Print the matrix to standard output via printf
 */
void adtk_matrix::print() const{
	print("%8.4f");		// Use the more advanced function
}//============================================

/**
 *	@brief Print the matrix to standard output with the specified format
 *	@param format a format string for the std::printf function
 */
void adtk_matrix::print(const char *format) const{
	for (int r = 0; r < rows; r++){
	    for (int c = 0; c < cols; c++){
	    	printf(format, gsl_matrix_get(a, r, c));
	    }
	    printf("\n");
	}
}//==============================================

/**
 *	@brief Write the matrix data out to a CSV file (must specify the extension)
 *	@param filename the filename or filepath, including file extension
 */
void adtk_matrix::toCSV(const char *filename) const{
	ofstream outFile;
	outFile.open(filename, ios::out);
	for (int r = 0; r < rows; r++){
	    for (int c = 0; c < cols; c++){
	    	char buffer[64];
	    	if(c < cols-1)
	    		sprintf(buffer, "%.14f, ", gsl_matrix_get(a,r,c));
	    	else
	    		sprintf(buffer, "%.14f\n", gsl_matrix_get(a,r,c));

	    	outFile << buffer;
	    }
	}

	outFile.close();
}
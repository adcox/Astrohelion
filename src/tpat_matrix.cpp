/**
 *	@file tpat_matrix.cpp
 *
 * 	Matrix object that employs GSL's Matrix and CBLAS functionality
 *
 *	The matrix is double precision, rectangular, and real. No other
 *	characteristics are assumed (e.g. banded, symmetric, etc.)
 *
 *	@author Andrew Cox
 *	@version May 6, 2015
 */
#include "tpat_matrix.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_utilities.hpp"

#include <cmath>
#include <fstream>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>
#include <iostream>
 
using namespace std;

//---------------------------------------------------------------------------
//    	Constructors and Destructor
//---------------------------------------------------------------------------

/**
 *	@brief Construct a matrix full of zeros
 *	@param r number of rows
 *	@param c number of columns
 */
tpat_matrix::tpat_matrix(const int r, const int c){
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
tpat_matrix::tpat_matrix(int r, int c, double *data){
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
tpat_matrix::tpat_matrix(int r, int c, std::vector<double> data){
	initBasicMatrix(r, c);
	copyDataIntoGSL_Matrix(data);
}//==================================

/**
 *	@brief Construct a matrix from a <tt>gsl_matrix</tt> object
 *	@param m a pointer to a <tt>gsl_matrix</tt>
 */
tpat_matrix::tpat_matrix(gsl_matrix *m){
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
tpat_matrix::tpat_matrix(gsl_vector *v, bool isRow){
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
tpat_matrix::~tpat_matrix(){
	try{
		gsl_matrix_free(a);
	}catch(...){
		cout << "\nError! cannot free a: my dimensions are " << rows << " x " << cols << endl;
	}
}//==============================================

/**
 *	@brief Copy constructor: creates a new matrix from this one without using the same memory!
 */
tpat_matrix::tpat_matrix(const tpat_matrix &b){
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
tpat_matrix tpat_matrix::I(int size){
	vector<double> data(size*size, 0);
	for(int i = 0; i < size*size; i += size+1){
		data[i] = 1;
	}
	return tpat_matrix(size, size, data);
}//==============================================

/**
 *	@brief create a simple basis vector \f$ e_j \f$
 *	
 *	For example, e_j(4, 1) will create the row vector {1, 0, 0, 0} and
 *	e_j(3,2) will create the row vector {0, 1, 0}
 *
 *	@param dim the dimension of the vector space (i.e. number of elements in the vector)
 *	@param j the vector number (begins with 1)
 *	@return the row vector
 */
tpat_matrix tpat_matrix::e_j(int dim, int j){
	vector<double> data(dim, 0);
	data[j-1] = 1;
	return tpat_matrix(1, dim, data);
}//===========================================

/**
 *	@brief Create a diagonal, square matrix with the specified elements
 *
 *	@param vals an array of values along the diagonal
 *	@param vals_len the number of values in vals
 *	@return a matrix with the specified values along the diagonal, zeros elsewhere
 */
tpat_matrix tpat_matrix::diag(double* vals, int vals_len){
	vector<double> data(vals_len*vals_len, 0);
	for(int i = 0; i < vals_len; i++){
		data[(vals_len+1)*i] = vals[i];
	}
	return tpat_matrix(vals_len, vals_len, data);
}//===========================================

//---------------------------------------------------------------------------
//    	Set and Get
//---------------------------------------------------------------------------

/**
 * 	@return the number of rows in this matrix
 */
int tpat_matrix::getRows() const {return rows;}

/**
 *	@brief Get a row from this matrix
 *	@param i the row index (begins with zero)
 *	@return a 1-D matrix (row vector) that contains the data from row <tt>i</tt>
 */
tpat_matrix tpat_matrix::getRow(int i){
	// Retrieve the row data
	gsl_vector *v = gsl_vector_alloc(static_cast<size_t>(cols));
	gsl_matrix_get_row(v, a, static_cast<size_t>(i));
	tpat_matrix temp(v, true);
	gsl_vector_free(v);
	return temp;
}

/**
 * 	@return the number of columns in this matrix
 */
int tpat_matrix::getCols() const {return cols;}

/**
 *	@brief Get a column from this matrix
 *	@param j the column index (begins with zero)
 *	@return a 1-D matrix (column vector) that contains the data from column <tt>j</tt>
 */
tpat_matrix tpat_matrix::getCol(int j){
	gsl_vector *v = gsl_vector_alloc(static_cast<size_t>(rows));
	gsl_matrix_get_col(v, a, static_cast<size_t>(j));
	tpat_matrix temp(v, false);
	gsl_vector_free(v);
	return temp;
}

/**
 *	@return the total number of elements in the matrix
 */
int tpat_matrix::getLength() const {return rows*cols;}

/**
 *	@return a pointer to the gsl_matrix object that represents this matrix
 */
gsl_matrix* tpat_matrix::getGSLMat(){return a;}

/**
 *	@return a pointer to the array of doubles that represents this matrix;
 *	values are stored in row-major order
 */
double* tpat_matrix::getDataPtr(){ return a->data; }

/**
 *	@param r the row
 *	@param c the column
 *	@return the value in the specified row and column
 */
double tpat_matrix::at(int r, int c) const {
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
double tpat_matrix::at(int i) const {
	if(i < rows*cols){
		if(rows != 1 && cols != 1){
			cout << "tpat_matrix :: Warning - matrix is not 1D; returning element from storage array" << endl;
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
tpat_matrix& tpat_matrix::operator =(const tpat_matrix &b){
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
tpat_matrix operator +(const tpat_matrix &lhs, const tpat_matrix &rhs){
	if(lhs.rows == rhs.rows && lhs.cols == rhs.cols){
		vector<double> q(lhs.rows * lhs.cols);
		for (int r = 0; r < lhs.rows; r++){
			for (int c = 0; c < lhs.cols; c++){
				q[r * lhs.cols + c] = gsl_matrix_get(lhs.a, r, c) + gsl_matrix_get(rhs.a, r, c);
			}
		}
		// Return a new matrix
		return tpat_matrix(lhs.rows, lhs.cols, q);
	}else{
		cout << "Matrices must be the same size to apply addition!" << endl;
		throw tpat_sizeMismatch();
	}
}

/**
 *	@brief Add a matrix to this one; perform operation in place
 *	@param b a matrix, must be the same size as this matrix
 *	@return this matrix
 */
tpat_matrix& tpat_matrix::operator +=(const tpat_matrix &b){
	if(rows == b.rows && cols == b.cols){
		gsl_matrix_add(a, b.a);
		return *this;
	}else{
		cout << "Matrices must be the same size to apply addition!" << endl;
		throw tpat_sizeMismatch();
	}
}//==============================================

/**
 *	@brief Subtract two matrices: A = L - R
 *	@param lhs the left-hand matrix L
 *	@param rhs the right-hand matrix R
 *	@return the difference L - R
 */
tpat_matrix operator -(const tpat_matrix &lhs, const tpat_matrix &rhs){
	return operator +(lhs, -1*rhs);
}//==============================================

/**
 *	@brief In place subtraction; subtract a matrix from this one
 *	@param b a matrix to subtract from this
 *	@return this matrix
 */
tpat_matrix& tpat_matrix::operator -=(const tpat_matrix &b){
	if(rows == b.rows && cols == b.cols){
		gsl_matrix_sub(a, b.a);
		return *this;
	}else{
		cout << "Matrices must be the same size to apply addition!" << endl;
		throw tpat_sizeMismatch();
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
tpat_matrix tpat_matrix::operator *(const tpat_matrix &b) const{
	if(cols != b.rows){
		throw tpat_sizeMismatch();
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
	tpat_matrix ans(rows, b.cols, C);
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
tpat_matrix& tpat_matrix::operator *=(const tpat_matrix &b){

	if(cols == b.rows && cols == b.cols){
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
		printErr("Operation will not return same size matrix\n");
		throw tpat_sizeMismatch();
	}
}//=============================================

/**
 *	@brief Multiply a scalar and a matrix: A = l*R
 *	@param lhs the scalar l
 *	@param rhs the matrix R
 *	@return the product l*R (same as R*l)
 */
tpat_matrix operator *(const double &lhs, const tpat_matrix &rhs){
	int rows = rhs.rows;
	int cols = rhs.cols;
	vector<double> b(rows * cols);
	for(int r = 0; r< rows; r++){
		for(int c = 0; c< cols; c++){
			b[r * cols + c] = lhs *gsl_matrix_get(rhs.a, r, c);
		}
	}

	return tpat_matrix(rows, cols, b);
}//=============================================

/**
 *	@brief Multiply a scalar and a matrix: A = L*r
 *	@param lhs the matrix L
 *	@param rhs the scalar r
 *	@return the product L*r (same as r*L)
 */
tpat_matrix operator *(const tpat_matrix &lhs, const double &rhs){
	return operator *(rhs, lhs);
}//=============================================

/**
 *	@brief Multiply a scalar and a matrix: A = L/r
 *	@param lhs the matrix L
 *	@param rhs the scalar r
 *	@return the quotient L/r (same as 1/r * L)
 */
tpat_matrix operator /(const tpat_matrix &lhs, const double &rhs){
	return operator *(lhs, 1/rhs);
}

/**
 *	@brief Multiply this matrix by a scalar (in place, DOES modify this matrix)
 * 	@param alpha a real scalar
 */
tpat_matrix& tpat_matrix::operator *=(const double &alpha){
	gsl_matrix_scale(a, alpha);
	return *this;
}//==============================================

/**
 *	@brief Compare two matrices
 *	@param lhs the left-hand matrix L
 *	@param rhs the right-hand matrix R
 *	@return whether or not L and R contain the same values
 */
bool operator ==(const tpat_matrix &lhs, const tpat_matrix &rhs){
	return gsl_matrix_equal(lhs.a, rhs.a);
}

/**
 *	@brief Compare two matrices
 *	@param lhs the left-hand matrix L
 *	@param rhs the right-hand matrix R
 *	@return whether or not L and R do NOT contain the same values
 */
bool operator !=(const tpat_matrix &lhs, const tpat_matrix &rhs){
	return !operator==(lhs, rhs);
}

//---------------------------------------------------------------------------
//    	Matrix Operation Functions
//---------------------------------------------------------------------------

/**
 *	@brief Compute the L2, or Euclidean, norm of a 1D matrix (i.e. a vector)
 *
 *	@param m a 1D matrix
 */
double abs(const tpat_matrix &m){
	return norm(m);
}//==============================

/**
 *	@brief compute the cross product of two, three-element vectors
 *	@param lhs a 3-element vector
 *	@param rhs a 3-element vector
 *	@return \f$ \vec{l} \times \vec{r} \f$ as a 3-element row vector,
 *	regardless of the orientation of the input vectors
 */
tpat_matrix cross(const tpat_matrix &lhs, const tpat_matrix &rhs){
	if(lhs.rows*lhs.cols != 3 || rhs.rows*rhs.cols != 3){
		printErr("Cross product only defined for 3D vectors...\n");
		throw tpat_linalg_err();
	}

	vector<double> data(3,0);
	data[0] = lhs.at(1)*rhs.at(2) - lhs.at(2)*rhs.at(1);
	data[1] = rhs.at(0)*lhs.at(2) - rhs.at(2)*lhs.at(0);
	data[2] = lhs.at(0)*rhs.at(1) - lhs.at(1)*rhs.at(0);

	return tpat_matrix(1, 3, data);
}//============================================

/**
 *	@brief compute the determinant of a square matrix
 *
 *	This function uses GSL's linear algebra library to compute the determinant
 *
 *	@param m a square matrix
 *	@return the determinant of <tt>m</tt>
 */
double det(const tpat_matrix &m){
	if(m.rows != m.cols){
		printErr("Cannot take determinant of non-square matrix!\n");
		throw tpat_linalg_err();
	}
	int permSign;
	gsl_permutation *perm = gsl_permutation_alloc(m.rows);
	int status = gsl_linalg_LU_decomp(m.a, perm, &permSign);
	if(status){
		printErr("Unable to perform LU decomp on matrix; GSL ERR: %s\n", gsl_strerror(status));
		throw tpat_linalg_err();
	}
	double ans = gsl_linalg_LU_det(m.a, permSign);
	gsl_permutation_free(perm);

	return ans;
}//========================

/**
 *	@brief Compute the L2, or Euclidean, norm of a 1D matrix (i.e. a vector)
 *
 *	@param m a 1D matrix
 */
double norm(const tpat_matrix &m){
	if(m.rows == 1 || m.cols == 1){
		double sumSquares = 0;
		for(int i = 0; i < m.rows * m.cols; i++){
			sumSquares += (m.a->data[i])*(m.a->data[i]);
		}
		return sqrt(sumSquares);
	}else{
		cout << "tpat_matrix :: Cannot take norm of a 2-D matrix" << endl;
		throw tpat_exception();
	}
}//==========================================

/**
 *	@brief Transpose a matrix
 *	@param m a matrix
 *	@return the transpose (NOT the complex conjugate transpose)
 */
tpat_matrix trans(const tpat_matrix &m){
	gsl_matrix *b = gsl_matrix_alloc(m.cols, m.rows);
	gsl_matrix_transpose_memcpy(b, m.a);
	tpat_matrix temp(m.cols, m.rows, b->data);
	gsl_matrix_free(b);
	return temp;
}//===============================================

//---------------------------------------------------------------------------
//    	Utilities
//---------------------------------------------------------------------------

/**
 *	@brief Initialize basic matrix and set all values to zero
 *	@param r number of rows
 *	@param c number of columns
 */
void tpat_matrix::initBasicMatrix(int r, int c){
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
void tpat_matrix::copyDataIntoGSL_Matrix(double *data){
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
void tpat_matrix::copyDataIntoGSL_Matrix(std::vector<double> v){
	if(((int)v.size()) != rows*cols){
		printErr("Input vector size does not match matrix size!\n");
		throw tpat_sizeMismatch();
	}
	copyDataIntoGSL_Matrix(&(v[0]));
}//============================================

/**
 * 	@brief Print the matrix to standard output via printf
 */
void tpat_matrix::print() const{
	print("%8.4f");		// Use the more advanced function
}//============================================

/**
 *	@brief Print the matrix to standard output with the specified format
 *	@param format a format string for the std::printf function
 */
void tpat_matrix::print(const char *format) const{
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
void tpat_matrix::toCSV(const char *filename) const{
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
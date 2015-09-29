/**
 *	@file tpat_matrix.cpp
 *
 * 	@brief Matrix object that employs GSL's Matrix and CBLAS functionality
 *
 *	The matrix is double precision and real. No other
 *	characteristics are assumed (e.g. square, banded, symmetric, etc.)
 *
 *	@author Andrew Cox
 *	@version May 6, 2015
 */
/*
 *	Trajectory Propagation and Analysis Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
 *
 *  TPAT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  TPAT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with TPAT.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "tpat.hpp"

#include "tpat_matrix.hpp"

#include "tpat_ascii_output.hpp"
#include "tpat_exceptions.hpp"
#include "tpat_utilities.hpp"

#include <cmath>
#include <fstream>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <Eigen/Dense>

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
	initBasicMatrix(rows, cols);
	gsl_matrix_memcpy(gslMat, m);
	// gslMat = m;
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
		gsl_matrix_free(gslMat);
	}catch(...){
		throw tpat_exception("Cannot free gsl matrix!");
	}
}//==============================================

/**
 *	@brief Copy constructor: creates a new matrix from this one without using the same memory!
 */
tpat_matrix::tpat_matrix(const tpat_matrix &b){
	initBasicMatrix(b.rows, b.cols);
	for(int r=0; r<rows; r++){
		for(int c=0; c<cols; c++){
			gsl_matrix_set(gslMat, r, c, gsl_matrix_get(b.gslMat, r, c));
		}
	}
}//=================================================

/**
 * 	@brief Create an identity matrix
 *	@param size the number of rows or columns the (square) matrix has
 *	@return an identiy matrix
 */
tpat_matrix tpat_matrix::I(int size){
	std::vector<double> data(size*size, 0);
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
	std::vector<double> data(dim, 0);
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
	std::vector<double> data(vals_len*vals_len, 0);
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
	gsl_matrix_get_row(v, gslMat, static_cast<size_t>(i));
	tpat_matrix temp(v, true);
	gsl_vector_free(v);
	return temp;
}//==================================

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
	gsl_matrix_get_col(v, gslMat, static_cast<size_t>(j));
	tpat_matrix temp(v, false);
	gsl_vector_free(v);
	return temp;
}//======================================

/**
 *	@return the total number of elements in the matrix
 */
int tpat_matrix::getLength() const {return rows*cols;}

/**
 *	@return a pointer to the gsl_matrix object that represents this matrix
 */
gsl_matrix* tpat_matrix::getGSLMat(){return gslMat;}

/**
 *	@return a pointer to the array of doubles that represents this matrix;
 *	values are stored in row-major order
 */
double* tpat_matrix::getDataPtr(){ return gslMat->data; }

/**
 *	@param r the row
 *	@param c the column
 *	@return the value in the specified row and column
 */
double tpat_matrix::at(int r, int c) const {
	if(r >= 0 && r < rows && c >= 0 && c < cols){
		return gsl_matrix_get(gslMat, r, c);
	}else{
		throw std::out_of_range("Index out of range");
	}
}//=================================

/**
 *	@brief Retrive a value from a 1D matrix. If the matrix is not 1D, a value from the first row
 *	is returned.
 *	@param i element index
 */
double tpat_matrix::at(int i) const {
	if(i >= 0 && i < rows*cols){
		if(rows != 1 && cols != 1){
			printWarn("tpat_matrix :: matrix is not 1D; returning element from storage array\n");
		}
		return gslMat->data[i];
	}else{
		throw std::out_of_range("Index out of range");
	}
}//===========================================

/**
 *	@brief Set a matrix element equal to a value
 *	@param r row index
 *	@param c column index
 *	@param val the value to store in (r, c)
 */
void tpat_matrix::set(int r, int c, double val){
	if(r >= 0 && r < rows && c >= 0 && c < cols){
		gsl_matrix_set(gslMat, r, c, val);
	}else{
		throw std::out_of_range("Index out of range");
	}
}//===========================================

/**
 *	@brief Set a value in a 1D matrix. If the matrix is not 1D, a value from the first row is set
 *	@param ix matrix/vector index
 *	@param val value to store
 */
void tpat_matrix::set(int ix, double val){
	if(ix >= 0 && ix < rows*cols){
		if(rows != 1 && cols != 1){
			printWarn("tpat_matrix :: matrix is not 1D; setting a row element by default\n");
		}

		if(cols == 1)
			gsl_matrix_set(gslMat, ix, 1, val);
		else{
			if(ix < cols)
				gsl_matrix_set(gslMat, 1, ix, val);
			else
				throw std::out_of_range("Index out of range");		
		}
	}else{
		throw std::out_of_range("Index out of range");
	}
}//===========================================

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
			gsl_matrix_set(gslMat, r, c, gsl_matrix_get(b.gslMat, r, c));
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
		std::vector<double> q(lhs.rows * lhs.cols);
		for (int r = 0; r < lhs.rows; r++){
			for (int c = 0; c < lhs.cols; c++){
				q[r * lhs.cols + c] = gsl_matrix_get(lhs.gslMat, r, c) + gsl_matrix_get(rhs.gslMat, r, c);
			}
		}
		// Return a new matrix
		return tpat_matrix(lhs.rows, lhs.cols, q);
	}else{
		throw tpat_sizeMismatch("Matrices must be the same size to apply addition");
	}
}

/**
 *	@brief Add a matrix to this one; perform operation in place
 *	@param b a matrix, must be the same size as this matrix
 *	@return this matrix
 */
tpat_matrix& tpat_matrix::operator +=(const tpat_matrix &b){
	if(rows == b.rows && cols == b.cols){
		gsl_matrix_add(gslMat, b.gslMat);
		return *this;
	}else{
		throw tpat_sizeMismatch("Matrices must be the same size to apply addition");
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
		gsl_matrix_sub(gslMat, b.gslMat);
		return *this;
	}else{
		throw tpat_sizeMismatch("Cannot perform addition");
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
	double *A = gslMat->data;
	double *B = b.gslMat->data;

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
		double *A = gslMat->data;
		double *B = b.gslMat->data;
		double *C = new double[rows * b.cols];

		// Do the multiplication, 
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, cols, cols,
			1.0, A, cols, B, cols, 0.0, C, cols);

		// Store data back in this matrix
		for (int i = 0; i < rows; i++){
		    for (int j = 0; j < cols; j++){
		    	gsl_matrix_set(gslMat, i, j, C[i*cols+j]);
		    }
		}

		delete[] C;
		return *this;
	}else{
		throw tpat_sizeMismatch("Operation will not return same size matrix");
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
	std::vector<double> b(rows * cols);
	for(int r = 0; r< rows; r++){
		for(int c = 0; c< cols; c++){
			b[r * cols + c] = lhs *gsl_matrix_get(rhs.gslMat, r, c);
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
	gsl_matrix_scale(gslMat, alpha);
	return *this;
}//==============================================

/**
 *	@brief Compare two matrices
 *	@param lhs the left-hand matrix L
 *	@param rhs the right-hand matrix R
 *	@return whether or not L and R contain the same values
 */
bool operator ==(const tpat_matrix &lhs, const tpat_matrix &rhs){
	return gsl_matrix_equal(lhs.gslMat, rhs.gslMat);
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
		throw tpat_linalg_err("Cross product only defined for 3D vectors");
	}

	std::vector<double> data(3,0);
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
		throw tpat_linalg_err("Cannot take determinant of non-square matrix");
	}
	int permSign;
	gsl_permutation *perm = gsl_permutation_alloc(m.rows);
	int status = gsl_linalg_LU_decomp(m.gslMat, perm, &permSign);
	if(status){
		printErr("tpat_matrix::det: GSL ERR: %s\n", gsl_strerror(status));
		throw tpat_linalg_err("Unable to perform LU decomp on matrix");
	}
	double ans = gsl_linalg_LU_det(m.gslMat, permSign);
	gsl_permutation_free(perm);

	return ans;
}//============================================

/**
 *	@brief Compute the eigenvalues of a square matrix
 *	@param m an nxn matrix reference
 *	@return a vector of vectors; the first vector contains n eigenvalues,
 *	while the second vector contains n, n-element, complex eigenvectors in
 *	row-major order. The eigenvalues and eigenvectors are listed in the
 *	same order, so they correspond to one another.
 */	
std::vector< std::vector<cdouble> > eig(const tpat_matrix &m){
	if(m.rows != m.cols)
		throw tpat_sizeMismatch("Cannot compute eigenvalues and vectors of non-square matrix");

	// Space for eigenvalues and eigenvectors
	gsl_vector_complex *eval = gsl_vector_complex_alloc (m.rows);
  	gsl_matrix_complex *evec = gsl_matrix_complex_alloc (m.rows, m.cols);

  	// Allocate space to do computations
  	gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (m.rows);
  
  	// Compute eigenvalues of non-symmetric, real matrix
  	gsl_eigen_nonsymmv (m.gslMat, eval, evec, w);

  	std::vector<cdouble> eigenVals;
  	std::vector<cdouble> eigenVecs;
  	for(int i = 0; i < m.rows; i++){
  		gsl_complex val = gsl_vector_complex_get(eval, i);
  		eigenVals.push_back(cdouble(GSL_REAL(val), GSL_IMAG(val)));

  		for(int j = 0; j < m.rows; j++){
  			gsl_complex vec_el = gsl_matrix_complex_get(evec, j, i);
  			eigenVecs.push_back(cdouble(GSL_REAL(vec_el), GSL_IMAG(vec_el)));
  		}
  	}

  	std::vector< std::vector<cdouble> > eigData;
  	eigData.push_back(eigenVals);
  	eigData.push_back(eigenVecs);

  	// Free memory
  	gsl_eigen_nonsymmv_free (w);
  	gsl_vector_complex_free(eval);
  	gsl_matrix_complex_free(evec);

  	return eigData;
}//=========================================

/**
 *	@brief Invert a matrix
 *
 *	Only square matrices can be inverted
 *	@param m the matrix to invert
 *	@return the inverted matrix
 */
tpat_matrix inv(const tpat_matrix& m){
	if(m.rows != m.cols)
		throw tpat_linalg_err("Cannot invert non-square matrix");

	gsl_permutation *perm = gsl_permutation_alloc(m.rows);
	int permSign = 0;

	int status = gsl_linalg_LU_decomp(m.gslMat, perm, &permSign);
	if(status){
		printErr("tpat_matrix::inv GSL ERR: %s\n", gsl_strerror(status));
		throw tpat_linalg_err("Unable to decompose matrix into L and U");
	}
	gsl_matrix *invMat = gsl_matrix_alloc(m.rows, m.cols);
	status = gsl_linalg_LU_invert(m.gslMat, perm, invMat);
	if(status){
		printErr("tpat_matrix::inv GSL ERR: %s\n", gsl_strerror(status));
		throw tpat_linalg_err("Unable to invert LU decomposition");
	}

	tpat_matrix tpatInvMat(invMat);

	gsl_permutation_free(perm);
	gsl_matrix_free(invMat);
	return tpatInvMat;
}//============================================

/**
 *	@brief Compute the L2, or Euclidean, norm of a 1D matrix (i.e. a vector)
 *
 *	@param m a 1D matrix
 */
double norm(const tpat_matrix &m){
	if(m.rows == 1 || m.cols == 1){
		double sumSquares = 0;
		for(int i = 0; i < m.rows * m.cols; i++){
			sumSquares += (m.gslMat->data[i])*(m.gslMat->data[i]);
		}
		return sqrt(sumSquares);
	}else{
		throw tpat_exception("Cannot take norm of a 2-D matrix");
	}
}//==========================================

/**
 *	@brief Transpose a matrix
 *	@param m a matrix
 *	@return the transpose (NOT the complex conjugate transpose)
 */
tpat_matrix trans(const tpat_matrix &m){
	gsl_matrix *b = gsl_matrix_alloc(m.cols, m.rows);
	gsl_matrix_transpose_memcpy(b, m.gslMat);
	tpat_matrix temp(b);
	gsl_matrix_free(b);
	return temp;
}//===============================================

/**
 *	@brief Compute the rank of a matrix
 *
 *	This algorithm works by computing the QR decomposition and then
 *	counting the number of rows in R that are equal to zero; subtracting
 *	the number of zero rows from the total gives the rank of the matrix.
 *	Note that 'equal to zero' is a numerical computation, so some tolerance
 *	is required. 1e-12 or 1e-14 are generally acceptable values for zero
 *
 *	@param m a matrix reference
 *	@param okErr the acceptable numerical error.
 */
int mat_rank(const tpat_matrix &m, double okErr){
	// It may be prudent to compute a value for okErr from the magnitude of
	// the elements stored in M. R contains elements from the Gram-Schmidt orthogonalization

	// Make a copy so as not to change the input
	tpat_matrix temp(m);

	int minDim = m.rows < m.cols ? m.rows : m.cols;

	gsl_vector *tau = gsl_vector_alloc(minDim);
	int status = gsl_linalg_QR_decomp(temp.getGSLMat(), tau);
	if(status){
		printErr("tpat_matrix::rank QR decomp GSL ERR: %s\n", gsl_strerror(status));
		gsl_vector_free(tau);
		throw tpat_linalg_err("tpat_matrix::rank: Unable to decompose matrix into Q and R");
	}

	gsl_matrix *Q = gsl_matrix_alloc(m.rows, m.rows);
	gsl_matrix *R = gsl_matrix_alloc(m.rows, m.cols);
	status = gsl_linalg_QR_unpack(temp.getGSLMat(), tau, Q, R);
	if(status){
		printErr("tpat_matrix::rank QR unpack GSL ERR: %s\n", gsl_strerror(status));
		gsl_vector_free(tau);
		gsl_matrix_free(Q);
		throw tpat_linalg_err("tpat_matrix::rank: Unable to unpack Q and R");
	}

	// Count number of rows that are zero (within tolerance), starting at bottom
	tpat_matrix R_mat(R);
	int rank = m.rows;
	int r = m.rows-1;
	while(r >= 0 && std::abs(norm(R_mat.getRow(r))) < okErr){
		rank -= 1;
		r -= 1;
	}

	// Free stuff (Freeing R seems to cause an error when temp is destructed... I suspect R points to the temp gslMat)
	gsl_vector_free(tau);
	gsl_matrix_free(Q);
	return rank;
}//==================================================

/**
 *	@brief Compute the right nullspace of the input matrix using QR decomposition
 *	@param m a matrix reference
 *	@retrun the right nullspace of m
 *	@see null_qr(const tpat_matrix&, double)
 */
tpat_matrix null_qr(const tpat_matrix &m){
	return null_qr(m, 1e-14);
}//================================================

/**
 *	@brief Compute an orthonormal basis for the right nullspace of a vector
 *
 *	This method uses the following algorithm to compute the nullspace:
 *	<br>
 *	An \f$ m \times n \f$ matrix \f$ \bf{A} \f$ can be factored into two matrices,
 *	\f$ \bf{Q} \f$ and \f$ \bf{R} \f$ such that \f$ \bf{A} = \bf{QR} \f$. The
 *	range and nullspace of \f$\bf{A}\f$ are stored within \f$ \bf{Q} \f$ as
 *	\f$ \bf{Q} = [\bf{Q_1},~\bf{Q_2}] \f$ where \f$\bf{Q1}\f$ is \f$ m \times m \f$
 *	and \f$\bf{Q_2}\f$ is \f$ m \times (n-r) \f$, \f$r\f$ being the rank of \f$\bf{Q}\f$.
 *	The columns of \f$\bf{Q_1}\f$ span the range of \f$\bf{A}\f$ and the columns of 
 *	\f$\bf{Q_2}\f$ span the nullspace.
 *
 *	@param m a matrix reference
 *	@param okErr the maximum acceptable numerical error when looking for zeros
 *	@return a matrix whose columns span the nullspace of m. If the nullspace is empty,
 *	a matrix with 1 element set equal to 0 will be returned, regardless of the size of
 *	the input matrix.
 */
tpat_matrix null_qr(const tpat_matrix &m, double okErr){
	tpat_matrix temp = trans(m);

	int tempRows = temp.getRows();
	int tempCols = temp.getCols();
	int minDim = tempRows < tempCols ? tempRows : tempCols;
	
	gsl_vector *tau = gsl_vector_alloc(minDim);
	int status = gsl_linalg_QR_decomp(temp.getGSLMat(), tau);
	if(status){
		printErr("tpat_matrix::null QR decomp GSL ERR: %s\n", gsl_strerror(status));
		gsl_vector_free(tau);
		throw tpat_linalg_err("tpat_matrix::null: Unable to decompose matrix into Q and R");
	}

	gsl_matrix *Q = gsl_matrix_alloc(tempRows, tempRows);
	gsl_matrix *R = gsl_matrix_alloc(tempRows, tempCols);
	status = gsl_linalg_QR_unpack(temp.getGSLMat(), tau, Q, R);
	if(status){
		printErr("tpat_matrix::null QR upack GSL ERR: %s\n", gsl_strerror(status));
		gsl_vector_free(tau);
		gsl_matrix_free(Q);
		throw tpat_linalg_err("tpat_matrix::null: Unable to unpack Q and R");
	}

	// Count number of rows that are zero (within tolerance), starting at bottom
	tpat_matrix R_mat(R);
	int rank = tempRows;
	int r = tempRows-1;
	while(r >= 0 && std::abs(norm(R_mat.getRow(r))) < okErr){
		rank -= 1;
		r -= 1;
	}

	// Check to see if a nullspace even exists!
	if(rank == tempRows){
		double data = 0;
		return tpat_matrix(1,1,&data);
	}

	// Extract the columns of Q that correspond to the nullspace
	tpat_matrix N(m.cols, m.cols - rank);
	for(r = 0; r < m.cols; r++){
		for(int c = rank; c < m.cols; c++){
			N.set(r, c - rank, gsl_matrix_get(Q,r,c));
		}
	}

	// Free stuff (Freeing R seems to cause an error when temp is destructed... I suspect R points to the temp gslMat)
	gsl_vector_free(tau);
	gsl_matrix_free(Q);

	return N;
}//================================================

/**
 *	@brief Compute the right nullspace of the input matrix using SVD
 *	@param m a matrix reference
 *	@retrun the right nullspace of m
 *	@see null_svd(const tpat_matrix&, double)
 */
tpat_matrix null_svd(const tpat_matrix &m){
	return null_svd(m, 1e-14);
}//================================================

/**
 *	@brief Compute the right nullspace for an input matrix
 *
 *	This algorithm computes the nullspace from the Singular Value
 *	Decomposition of the input matrix. First, the rank is computed
 *	by counting the number of zero (within 1e-14) singular values
 *	and subtracting them from the total possible rank. The columns
 *	of V that correspond to those zero singular values will span
 *	the nullspace.
 *
 *	@param m a matrix reference
 *	@param okErr the acceptable error when identifying zero-valued singular values
 *	@return the right nullspace of the input matrix
 */
tpat_matrix null_svd(const tpat_matrix &m, double okErr){

	// GSL Algorithm can only handle matrices with rows >= cols
	if(m.rows >= m.cols){
		tpat_matrix U(m);

		gsl_matrix *V = gsl_matrix_alloc(U.getCols(), U.getCols());
		gsl_vector *S = gsl_vector_alloc(U.getCols());
		gsl_vector *w = gsl_vector_alloc(U.getCols());
		int status = 0;

		// The one-sided Jacboi orthogonalization method should give singular values to a higher level of accuracy for S
		// status = gsl_linalg_SV_decomp_jacobi(U.getGSLMat(),V,S);
		status = gsl_linalg_SV_decomp(U.getGSLMat(), V, S, w);

		if(status){
			printErr("tpat_matrix::null_svd: could not perform SVD decomposition; GSL ERR: %s\n", gsl_strerror(status));
			gsl_matrix_free(V);
			gsl_vector_free(S);
			gsl_vector_free(w);
			throw tpat_linalg_err("tpat_matrix::null_svd: Unable to decompose matrix into SVD");
		}

		// Count number of zeros in S
		int rank = S->size;
		for(int r = (int)S->size-1; r >= 0; r--){
			if(std::abs(gsl_vector_get(S, r)) < okErr)
				rank--;
		}
		
		// Check for full-rank case
		if(rank == (int)(S->size)){
			double data = 0;
			return tpat_matrix(1, 1, &data);
		}

		// Extract the right nullspace from V
		tpat_matrix N(U.getCols(), U.getCols()-rank);
		for(int r = 0; r < (int)(V->size1); r++){
			for(int c = rank; c < (int)(V->size2); c++){
				N.set(r, c-rank, gsl_matrix_get(V,r,c));
			}
		}

		if(m.rows < m.cols)
			N = trans(N);

		gsl_matrix_free(V);
		gsl_vector_free(S);
		gsl_vector_free(w);
		return N;
	}else{	// Use Eigen library
		
		// Copy my matrix into an Eigen matrix
		Eigen::MatrixXd A(m.rows, m.cols);
		for(int r = 0; r < m.rows; r++){
			for(int c = 0; c < m.cols; c++){
				A(r,c) = m.at(r,c);
			}
		}

		Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
		svd.setThreshold(okErr);

		Eigen::MatrixXd V = svd.matrixV();
		double rank = svd.rank();

		if(rank == V.cols()){
			double data = 0;
			return tpat_matrix(1, 1, &data);
		}

		// Eigen::MatrixXd N((int)(V.rows()), (int)(V.cols()-rank));
		tpat_matrix N((int)(V.rows()), (int)(V.cols()-rank));
		for(int r = 0; r < V.rows(); r++){
			for(int c = rank; c < V.cols(); c++){
				// N(r,c-rank) = V(r,c);
				N.set(r, c-rank, V(r,c));
			}
		}

		return N;
	}
}//==============================================

//---------------------------------------------------------------------------
//    	Utilities
//---------------------------------------------------------------------------

/**
 *	@brief Initialize basic matrix and set all values to zero
 *	@param r number of rows
 *	@param c number of columns
 */
void tpat_matrix::initBasicMatrix(int r, int c){
	gslMat = gsl_matrix_alloc (r, c);
	rows = r;
	cols = c;
	gsl_matrix_set_zero(gslMat);
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
	    	gsl_matrix_set(gslMat, i, j, data[i*cols+j]);
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
		throw tpat_sizeMismatch("tpat_matrix::copyDataIntoGSL_Matrix: Input vector size does not match matrix size");
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
	    	printf(format, gsl_matrix_get(gslMat, r, c));
	    }
	    printf("\n");
	}
}//==============================================

/**
 *	@brief Write the matrix data out to a CSV file (must specify the extension)
 *	@param filename the filename or filepath, including file extension
 */
void tpat_matrix::toCSV(const char *filename) const{
	std::ofstream outFile(filename, std::ios::out);
	// outFile.open(filename, std::ios::out);
	for (int r = 0; r < rows; r++){
	    for (int c = 0; c < cols; c++){
	    	char buffer[64] = { };
	    	if(c < cols-1)
	    		sprintf(buffer, "%.14f, ", gsl_matrix_get(gslMat,r,c));
	    	else
	    		sprintf(buffer, "%.14f\n", gsl_matrix_get(gslMat,r,c));

	    	outFile << buffer;
	    }
	}

	outFile.close();
}
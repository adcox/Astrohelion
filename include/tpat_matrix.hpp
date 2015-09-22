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

#ifndef H_TPAT_MATRIX
#define H_TPAT_MATRIX

#include "tpat_constants.hpp"	// for complex number typedefs

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <string>
#include <vector>

/**
 *	@brief A data object that represents a matrix. 
 *
 *	This class defines functions that 
 *	allow for easy-to-read matrix operations in code. Some of the functionality
 *	is based on GSL's matrix and CBLAS functions, but some is also hand-coded.
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 *	@copyright GNU GPL v3.0
 */
class tpat_matrix{

	public:

		// Constructors
		tpat_matrix(int r, int c);
		tpat_matrix(int r, int c, double *data);
		tpat_matrix(int r, int c, std::vector<double> data);
		tpat_matrix(gsl_matrix*);
		tpat_matrix(gsl_vector*, bool);
		tpat_matrix(const tpat_matrix &b);	// copy constructor
		~tpat_matrix();

		// Special constructors
		static tpat_matrix I(int);
		static tpat_matrix e_j(int, int);
		static tpat_matrix diag(double*, int);

		// Operator Overloading
		tpat_matrix& operator =(const tpat_matrix&);

		friend tpat_matrix operator +(const tpat_matrix&, const tpat_matrix&);
		tpat_matrix& operator +=(const tpat_matrix&);

		friend tpat_matrix operator -(const tpat_matrix&, const tpat_matrix&);
		tpat_matrix& operator -=(const tpat_matrix&);

		friend tpat_matrix operator *(const double&, const tpat_matrix&);
		friend tpat_matrix operator *(const tpat_matrix&, const double&);
		tpat_matrix operator *(const tpat_matrix&) const;
		tpat_matrix& operator *=(const double&);
		tpat_matrix& operator *=(const tpat_matrix&);

		friend tpat_matrix operator /(const tpat_matrix&, const double&);

		friend bool operator ==(const tpat_matrix&, const tpat_matrix&);
		friend bool operator !=(const tpat_matrix&, const tpat_matrix&);

		// Matrix operations
		friend double abs(const tpat_matrix&);
		friend tpat_matrix cross(const tpat_matrix&, const tpat_matrix&);
		friend double det(const tpat_matrix&);
		friend tpat_matrix inv(const tpat_matrix&);
		friend double norm(const tpat_matrix&);
		friend tpat_matrix null_qr(const tpat_matrix&, double);
		friend tpat_matrix null_svd(const tpat_matrix &);
		friend int mat_rank(const tpat_matrix&, double);
		friend tpat_matrix trans(const tpat_matrix&);
		friend std::vector< std::vector<cdouble> > eig(const tpat_matrix&);

		// Set and Gets
		int getLength() const;
		int getRows() const;
		int getCols() const;
		gsl_matrix* getGSLMat();
		double* getDataPtr();
		double at(int, int) const;
		double at(int) const;
		tpat_matrix getRow(int);
		tpat_matrix getCol(int);

		void set(int, int, double);
		void set(int, double);

		// Utility functions
		void print() const;
		void print(const char *format) const;
		void toCSV(const char *filename) const;
		
		// static tpat_sizeMismatch tpat_matSizeMismatch;
		
	private:
		gsl_matrix *gslMat = 0;		//!< GSL matrix object that stores matrix elements
		int rows = 0;				//!< Number of rows
		int cols = 0;				//!< Number of columns

		void initBasicMatrix(int r, int c);
		void copyDataIntoGSL_Matrix(double*);
		void copyDataIntoGSL_Matrix(std::vector<double>);
};

#endif
//END
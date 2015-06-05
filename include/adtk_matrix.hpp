/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (ADTK).
 *
 *  ADTK is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ADTK is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ADTK.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __ADTK_MATRIX__
#define __ADTK_MATRIX__

#include "adtk_exceptions.hpp"

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
class adtk_matrix{

	public:

		// Constructors
		adtk_matrix(int r, int c);
		adtk_matrix(int r, int c, double *data);
		adtk_matrix(int r, int c, std::vector<double> data);
		adtk_matrix(gsl_matrix*);
		adtk_matrix(gsl_vector*, bool);
		adtk_matrix(const adtk_matrix &b);	// copy constructor
		~adtk_matrix();

		// Special constructors
		static adtk_matrix I(int);
		static adtk_matrix e_j(int, int);
		static adtk_matrix diag(double[], int);

		// Operator Overloading
		adtk_matrix& operator =(const adtk_matrix&);

		friend adtk_matrix operator +(const adtk_matrix&, const adtk_matrix&);
		adtk_matrix& operator +=(const adtk_matrix&);

		friend adtk_matrix operator -(const adtk_matrix&, const adtk_matrix&);
		adtk_matrix& operator -=(const adtk_matrix&);

		friend adtk_matrix operator *(const double&, const adtk_matrix&);
		friend adtk_matrix operator *(const adtk_matrix&, const double&);
		adtk_matrix operator *(const adtk_matrix&) const;
		adtk_matrix& operator *=(const double&);
		adtk_matrix& operator *=(const adtk_matrix&);

		friend adtk_matrix operator /(const adtk_matrix&, const double&);

		friend bool operator ==(const adtk_matrix&, const adtk_matrix&);
		friend bool operator !=(const adtk_matrix&, const adtk_matrix&);

		// Matrix operations
		friend double abs(const adtk_matrix&);
		friend adtk_matrix cross(const adtk_matrix&, const adtk_matrix&);
		friend double det(const adtk_matrix&);
		friend double norm(const adtk_matrix&);
		friend adtk_matrix trans(const adtk_matrix&);
		

		// Set and Gets
		int getLength() const;
		int getRows() const;
		int getCols() const;
		gsl_matrix* getGSLMat();
		double* getDataPtr();
		double at(int, int) const;
		double at(int) const;
		adtk_matrix getRow(int);
		adtk_matrix getCol(int);

		// Utility functions
		void print() const;
		void print(const char *format) const;
		void toCSV(const char *filename) const;
		
		// static adtk_sizeMismatch adtk_matSizeMismatch;
		
	private:
		gsl_matrix *a;
		int rows, cols;

		void initBasicMatrix(int r, int c);
		void copyDataIntoGSL_Matrix(double*);
		void copyDataIntoGSL_Matrix(std::vector<double>);
};

#endif
//END
/**
 *	A data object that represents a matrix. This class defines functions that 
 *	allow for easy-to-read matrix operations in code. Some of the functionality
 *	is based on GSL's matrix and CBLAS functions, but some is also hand-coded.
 *
 *	Author: Andrew Cox
 *
 *	Version: May 15, 2015
 */

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
 *  along with ATDK.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __ADTK_MATRIX__
#define __ADTK_MATRIX__

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <string>
#include <vector>

// Forward Declarations

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
		static adtk_matrix Identity(int size);

		// Operator Overloading
		adtk_matrix& operator =(const adtk_matrix&);
		friend adtk_matrix operator +(const adtk_matrix&, const adtk_matrix&);
		friend adtk_matrix operator -(const adtk_matrix&, const adtk_matrix&);
		friend adtk_matrix operator *(const double&, const adtk_matrix&);
		friend adtk_matrix operator *(const adtk_matrix&, const double&);
		friend adtk_matrix operator /(const adtk_matrix&, const double&);
		
		adtk_matrix operator *(const adtk_matrix&);

		adtk_matrix& operator +=(const adtk_matrix&);
		adtk_matrix& operator -=(const adtk_matrix&);
		adtk_matrix& operator *=(const double&);
		adtk_matrix& operator *=(const adtk_matrix&);
		
		friend bool operator ==(const adtk_matrix&, const adtk_matrix&);
		friend bool operator !=(const adtk_matrix&, const adtk_matrix&);

		// Matrix operations
		adtk_matrix trans();
		double norm();

		// Set and Gets
		int getLength();
		int getRows();
		int getCols();
		gsl_matrix* getGSLMat();
		double get(int, int);
		adtk_matrix getRow(int);
		adtk_matrix getCol(int);

		// Utility functions
		void print();
		void print(const char *format);

	private:
		gsl_matrix *a;
		int rows, cols;

		void initBasicMatrix(int r, int c);
		void copyDataIntoGSL_Matrix(double*);
		void copyDataIntoGSL_Matrix(std::vector<double>);
};

#endif
//END
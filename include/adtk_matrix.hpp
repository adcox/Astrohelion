/**
 *	Headers for matrix object
 *	Author: Andrew Cox
 *	Version: May 6, 2015
 */

#ifndef __ADTK_MATRIX__
#define __ADTK_MATRIX__

#include <gsl/gsl_matrix.h>
#include <string>

class adtk_matrix{

	public:

		// Constructors
		adtk_matrix(int r, int c);
		adtk_matrix(int r, int c, double data[]);
		adtk_matrix(const adtk_matrix &b);	// copy constructor
		~adtk_matrix();

		// Special constructors
		static adtk_matrix Identity(int size);

		// Operator Overloading
		adtk_matrix& operator =(const adtk_matrix&);
		adtk_matrix operator +(const adtk_matrix&);
		adtk_matrix operator -(const adtk_matrix&);
		adtk_matrix operator *(const double&);
		adtk_matrix operator *(const adtk_matrix&);
		adtk_matrix& operator +=(const adtk_matrix&);
		adtk_matrix& operator -=(const adtk_matrix&);
		adtk_matrix& operator *=(const double&);
		adtk_matrix& operator *=(const adtk_matrix&);
		
		bool operator ==(const adtk_matrix&);
		bool operator !=(const adtk_matrix&);

		// Matrix operations
		adtk_matrix trans();

		// Set and Gets
		int getLength();
		int getRows();
		int getCols();
		gsl_matrix* getGSLMat();
		double get(int, int);

		// Utility functions
		void print();
		void print(const char *format);

	private:
		gsl_matrix *a;
		int rows, cols;

		void initBasicMatrix(int, int);
};

#endif
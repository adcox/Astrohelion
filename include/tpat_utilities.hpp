/**
 *	@brief Contains miscellaneous utility functions that make 
 *	coding in C++ easier
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 *	@copyright GNU GPL v3.0
 */

/*
 *	Astrodynamics Toolkit 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrodynamics Toolkit (TPAT).
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

#ifndef __H_UTILITIES_
#define __H_UTILITIES_

#include "tpat_exceptions.hpp"
 
#include "matio.h"
 
#include <complex>
#include <string>
#include <vector>


/**
 *	@brief sum all values in an array
 *
 *	For this function to work, the class must overload +=
 *
 *	@param data an array
 *	@param length the number of elements in the array that can be summed.
 *	@return a single object representing the sum of all the elements of data.
 */
template<class T>
T tpat_sum(T* data, int length){
	T total = 0;
	for(int n = 0; n < length; n++){
		total += data[n];
	}

	return total;
}//=====================================

/**
 *	@brief concatenate two vectors
 *	@param lhs the left-hand-side vector
 *	@param rhs the righ-hand-side vector
 *	@return a new vector that includes both input vectors in order, i.e. [lhs, rhs]
 */
template<class T>
std::vector<T> concatVecs(std::vector<T> lhs, std::vector<T> rhs){
	std::vector<T> tempVec(lhs.begin(), lhs.end());
	tempVec.insert(tempVec.end(), rhs.begin(), rhs.end());
	return tempVec;
}//=======================================

/**
 *	@brief sort a vector and retrieve the indices of the sorted elements
 *
 *	Takes a copy of a vector and sorts it, retaining the original indices
 *	of the elements. For example, sorting {1,3,2} would return the indices 
 *	{0,2,1}.
 *
 *	@param v a vector to sort
 *	@return the indices of the sorted elements
 */
template <typename T>
std::vector<int> getSortedInd(const std::vector<T> &v) {

  	// initialize original index locations
  	std::vector<int> idx(v.size());
  	for (int i = 0; i < ((int)idx.size()); i++){
		idx[i] = i;
  	}

  	// sort indexes based on comparing values in v
  	sort(idx.begin(), idx.end(), 
  		[&v](int i1, int i2) {return v[i1] < v[i2];} );

  	return idx; 
}//===========================================

template <typename T>
bool aboutEquals(T t1, T t2, double tol){
	return std::abs(t1 - t2) < tol;
}//=========================================

std::string complexToStr(std::complex<double> num);
void printErr(const char*, ...);
void printWarn(const char*, ...);
void printVerb(bool, const char*, ...);
void printColor(const char*, const char*, ...);
void printVerbColor(bool, const char*, const char*, ...);
void saveVar(mat_t*, matvar_t*, const char*, matio_compression);
int readIntFromMat(mat_t*, const char*, matio_types, matio_classes);
double readDoubleFromMat(mat_t*, const char*);
std::string readStringFromMat(mat_t*, const char* , matio_types, matio_classes);

void waitForUser();

#endif
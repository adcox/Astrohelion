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
template<class T> T tpat_sum(T* data, int length){
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
template<class T> std::vector<T> concatVecs(std::vector<T> lhs, std::vector<T> rhs){
	std::vector<T> tempVec(lhs.begin(), lhs.end());
	tempVec.insert(tempVec.end(), rhs.begin(), rhs.end());
	return tempVec;
}//=======================================

std::string complexToStr(std::complex<double> num);
void printErr(const char*, ...);
void printWarn(const char*, ...);
void printVerb(bool, const char*, ...);
void printColor(const char*, const char*, ...);
void printVerbColor(bool, const char*, const char*, ...);
void saveVar(mat_t*, matvar_t*, const char*, matio_compression);

void waitForUser();

#endif
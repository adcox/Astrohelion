/**
 *	@file tpat_utilities.hpp
 *	@brief Contains miscellaneous utility functions that make 
 *	coding in C++ easier
 *
 *	@author Andrew Cox
 *	@version May 15, 2015
 *	@copyright GNU GPL v3.0
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

#ifndef H_UTILITIES
#define H_UTILITIES

#include "tpat_constants.hpp"
#include "tpat_exceptions.hpp"
 
#include "matio.h"
 
#include <algorithm>
#include <complex>
#include <string>
#include <vector>

/**
 *	@brief Class containing useful utility-type functions
 *
 *	This class is not an "object", per se, but it stores static template functions
 *	that are useful for other tasks and allows the developer to specify the scope 
 * 	(public/private) of the methods if necessary
 */
class tpat_util{
private:
	/**
	 *	@brief Worker function to compute permutations via recursion
	 *	@param values a vector of all possible values for each element
	 *	@param numSpots the number of elements, or "spots"
	 *	@param ixs the indices of values for each element in the current permuation
	 *	@param perms a pointer to a storage vector for the permutations; stored in 
	 *	row-major order
	 */
	template <typename T>
	static void permute(std::vector<T> values, int numSpots, std::vector<int> ixs, std::vector<T> *perms){
		if((int)(ixs.size()) < numSpots){	// Not done with recursion; continue deeper!
			ixs.push_back(-1);
			for(size_t i = 0; i < values.size(); i++){
				ixs.back() = i;
				permute(values, numSpots, ixs, perms);
			}
		}else{ // We've reached the bottom, compute the permuation and return
			for(int i = 0; i < numSpots; i++){
				perms->push_back(values[ixs[i]]);
			}
			return;
		}
	}//==============================================

	/**
	 *	@brief Worker function to compute permutations via recursion without repeating values
	 *	@param values a vector of all possible values
	 *	@param ixs the indices of values for each element in the current permutation
	 *	@param perms a pointer to a storage vector for the permutations; stored in row-major order
	 */
	template <typename T>
	static void permute(std::vector<T> values, std::vector<int> ixs, std::vector<T> *perms){
		if(ixs.size() < values.size()){
			ixs.push_back(-1);
			for(size_t i = 0; i < values.size(); i++){
				bool alreadyUsed = false;
				for(size_t j = 0; j < ixs.size(); j++){
					if(ixs[j] == (int)i){
						alreadyUsed = true;
						break;
					}
				}

				if(!alreadyUsed){
					ixs.back() = i;
					permute(values, ixs, perms);
				}
			}
		}else{
			for(size_t i = 0; i < ixs.size(); i++){
				perms->push_back(values[ixs[i]]);
			}
			return;
		}
	}//==========================================

public:

	/**
	 *	@brief Check if two numbers are equal within a given tolerance
	 *	@param t1 a number
	 *	@param t2 another number
	 *	@param tol the tolerance
	 *	@return true if the absolute value of the difference between t1 
	 *	and t2 is less than the tolerance
	 */
	template <typename T>
	static bool aboutEquals(T t1, T t2, double tol){
		return std::abs(t1 - t2) < tol;
	}//=========================================

	/**
	 *	@ brief Check if two vectors are equal to a given tolerance
	 *	@param v1 a vector
	 *	@param v2 another vector
	 *	@param tol the tolerance
	 *	@return true if v1 and v2 are the same size and their corresponding
	 *	elements differ by less than the tolerance
	 */
	template <typename T>
	static bool aboutEquals(std::vector<T> v1, std::vector<T> v2, double tol){
		if(v1.size() != v2.size())
			return false;

		for(size_t i = 0; i < v1.size(); i++){
			if(!aboutEquals(v1[i], v2[i], tol))
				return false;
		}

		return true;
	}//==========================================

	/**
	 *	@brief concatenate two vectors
	 *	@param lhs the left-hand-side vector
	 *	@param rhs the righ-hand-side vector
	 *	@return a new vector that includes both input vectors in order, i.e. [lhs, rhs]
	 */
	template<class T>
	static std::vector<T> concatVecs(std::vector<T> lhs, std::vector<T> rhs){
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
	static std::vector<int> getSortedInd(const std::vector<T> &v) {

	  	// initialize original index locations
	  	std::vector<int> idx(v.size());
	  	for (int i = 0; i < ((int)idx.size()); i++){
			idx[i] = i;
	  	}

	  	// sort indexes based on comparing values in v
	  	std::sort(idx.begin(), idx.end(), 
	  		[&v](int i1, int i2) {return v[i1] < v[i2];} );

	  	return idx; 
	}//===========================================

	/**
	 *	@brief Generate all permutations of a set of n elements
	 *
	 *	Each of the n elements can contain any of the values stored in the values vector
	 *	@param values a vector containing all possible values for each element
	 *	@param n the number of elements in the permutation
	 *	
	 *	@return a vector containing all permutations, in row-major order where each row
	 *	is a separate permutation
	 */
	template <typename T>
	static std::vector<T> generatePerms(std::vector<T> values, int n){
		std::vector<int> ixs;
		std::vector<T> perms;
		permute(values, n, ixs, &perms);
		return perms;
	}//=========================================

	/**
	 *	@brief Generate all permutations of a set of values without repeating the values
	 *	@param values the set of values
	 *	@return a vector containing all permutations, in row-major order where each row
	 *	is a seperate permutation of <tt>values</tt>
	 */	
	template <typename T>
	static std::vector<T> generatePerms(std::vector<T> values){
		std::vector<int> ixs;
		std::vector<T> perms;
		permute(values, ixs, &perms);
		return perms;
	}

	/**
	 *	@brief Get the imaginary parts of every element of a vector
	 *	@param compVec a vector of complex numbers
	 *	@return a vector of the imaginary parts of the complex vector
	 */
	template<typename T>
	static std::vector<T> imag(std::vector<std::complex<T> > compVec){
		std::vector<T> imagParts;
		for(size_t i = 0; i < compVec.size(); i++)
			imagParts.push_back(std::imag(compVec[i]));

		return imagParts;
	}//================================================

	/**
	 *	@brief Compute the mean (average) of an array of data
	 *
	 *	@param data an array of numbers
	 *	@param length the length of the array
	 *
	 *	@return the mean
	 */
	template<typename T>
	static T mean(T *data, int length){
		return sum(data, length)/length;
	}//================================================

	/**
	 *  @brief Compute the mean (average) of a vector of data
	 * 
	 *  @param data a vector of numbers
	 *  @tparam T a numerical type, like int, double, long, etc.
	 *  @return the mean
	 */
	template<typename T>
	static T mean(std::vector<T> data){
		return mean(&(data[0]), data.size());
	}//================================================

	/**
	 *	@brief Get the real parts of every element of a vector
	 *	@param compVec a vector of complex numbers
	 *	@return a vector of the real parts of the complex vector
	 */
	template<typename T>
	static std::vector<T> real(std::vector<std::complex<T> > compVec){
		std::vector<T> realParts;
		for(size_t i = 0; i < compVec.size(); i++)
			realParts.push_back(std::real(compVec[i]));

		return realParts;
	}//================================================

	/**
	 *	@brief Get the sign of a number
	 *	@param num a number
	 *	@return +/- 1 for the sign (0 if T = 0)
	 */
	template<typename T>
	static int sign(T num){
		if(num == 0)
			return 0;
		else
			return num < 0 ? -1 : 1;
	}//===========================================
	
	/**
	 *	@brief sum all values in an array
	 *
	 *	For this function to work, the class must overload +=
	 *
	 *	@param data an array
	 *	@param length the number of elements in the array that can be summed.
	 *	@return a single object representing the sum of all the elements of data.
	 */
	template<typename T>
	static T sum(T* data, int length){
		T total = 0;
		for(int n = 0; n < length; n++){
			total += data[n];
		}

		return total;
	}//=====================================

	/**
	 *  @brief Sum all values in a vector
	 * 
	 *  @param data a vector of data
	 *  @tparam T numerical data type, like int, double, long, etc.
	 *  @return the sum
	 */
	template<typename T>
	static T sum(std::vector<T> data){
		return sum(&(data[0]), data.size());
	}//==============================================
};

std::string complexToStr(std::complex<double> num);
void printErr(const char*, ...);
void printWarn(const char*, ...);
void printVerb(bool, const char*, ...);
void printColor(const char*, const char*, ...);
void printVerbColor(bool, const char*, const char*, ...);
void saveMatrixToFile(const char*, const char*, std::vector<double>, size_t, size_t);
void saveVar(mat_t*, matvar_t*, const char*, matio_compression);
int readIntFromMat(mat_t*, const char*, matio_types, matio_classes);
double readDoubleFromMat(mat_t*, const char*);
std::vector<double> readMatrixFromMat(const char*, const char*);
std::string readStringFromMat(mat_t*, const char* , matio_types, matio_classes);

void waitForUser();

#endif
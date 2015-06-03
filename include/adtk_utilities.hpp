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

#ifndef __H_UTILITIES_
#define __H_UTILITIES_

template<class T> T adtk_sum(T* data, int length){
	T total = 0;
	for(int n = 0; n < length; n++){
		total += data[n];
	}

	return total;
}

void waitForUser();

void printErr(const char*, ...);
void printWarn(const char*, ...);
void printVerb(bool, const char*, ...);
void printColor(const char*, const char*, ...);
void printVerbColor(bool, const char*, const char*, ...);

#endif
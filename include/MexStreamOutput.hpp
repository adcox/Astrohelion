/*
 *	Copyright 2018, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of the Astrohelion GUI 
 *
 *  Astrohelion GUI is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Astrohelion GUI is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Astrohelion GUI.  If not, see <http://www.gnu.org/licenses/>.
 *  
 *  *************************************************************************
 *  This file contains functions that are useful for a variety of applications
 *  
 *  Author: Andrew Cox
 *  Version: September 25, 2018
 */
#pragma once

#include <ctype.h>
#include <string.h>
#include <iostream>

namespace astroHelper{

/* *****************************************************************************
 * 	Declarations
 **************************************************************************** */
class mystream;
class scoped_redirect_cout;

bool strcmpi(const char*, const char*);

/* *****************************************************************************
 * 	Definitions
 **************************************************************************** */

/**
 * @brief Compare two character arrays insensitive to case
 * 
 * @param a A character array
 * @param b Another character array
 * 
 * @return Whether or not the two arrays are identical if all the characters are
 * cast to lowercase (or uppercase)
 */
bool strcmpi(const char* a, const char* b){
	unsigned int lenA = strlen(a);
	unsigned int lenB = strlen(b);

	if(lenA != lenB)
		return false;

	for(lenA = 0; lenA < lenB; lenA++){
		if(tolower(a[lenA]) != tolower(b[lenA]))
			return false;
	}

	return true;
}//====================================================

class mystream : public std::streambuf {
	protected:
	virtual std::streamsize xsputn(const char *s, std::streamsize n){
		mexPrintf("%.*s", n, s);
		return n;
	}
	virtual int overflow(int c=EOF){
		if (c != EOF){
			mexPrintf("%.1s", &c);
		}
		return 1;
	}
};

class scoped_redirect_cout {
	public:
		scoped_redirect_cout() {
			old_buf = std::cout.rdbuf();
			std::cout.rdbuf(&mout);
		}
		
		~scoped_redirect_cout() {
			std::cout.rdbuf(old_buf);
		}
	private:
		mystream mout;
		std::streambuf *old_buf;
};

static scoped_redirect_cout mycout_redirect;

}
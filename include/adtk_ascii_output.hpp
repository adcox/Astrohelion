/**
 *	@brief This class contains macros for fancy ASCII outputs
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
 *  along with ATDK.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __H_ASCII_
#define __H_ASCII_

/* 
 *	ASCII Escape Sequences for text coloring. Usage:
 *		printf( RED "Here is some red text!" RESET);
 */


#define RESET   "\033[0m"					/**< Reset ASCII text to default values */
#define BLACK   "\033[30m"					/**< Make ASCII text black */
#define RED     "\033[31m"					/**< Make ASCII text red */
#define GREEN   "\033[32m"					/**< Make ASCII text green */
#define YELLOW  "\033[33m"					/**< Make ASCII text yellow */
#define BLUE    "\033[34m"					/**< Make ASCII text blue */
#define MAGENTA "\033[35m"					/**< Make ASCII text magenta */
#define CYAN    "\033[36m"					/**< Make ASCII text cyan */
#define WHITE   "\033[37m"					/**< Make ASCII text white */
#define BOLDBLACK   "\033[1m\033[30m"		/**< Make ASCII text bold and black */
#define BOLDRED     "\033[1m\033[31m"		/**< Make ASCII text bold and red */
#define BOLDGREEN   "\033[1m\033[32m"		/**< Make ASCII text bold and green */
#define BOLDYELLOW  "\033[1m\033[33m"		/**< Make ASCII text bold and yellow */
#define BOLDBLUE    "\033[1m\033[34m"		/**< Make ASCII text bold and blue */
#define BOLDMAGENTA "\033[1m\033[35m"		/**< Make ASCII text bold and magenta */
#define BOLDCYAN    "\033[1m\033[36m"		/**< Make ASCII text bold and cyan */
#define BOLDWHITE   "\033[1m\033[37m"		/**< Make ASCII text bold and white */

#endif

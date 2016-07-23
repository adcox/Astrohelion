/**
 *  @file Exceptions.hpp
 *	@brief Defines several custom exception classes
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
/*
 *	Astrohelion 
 *	Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *	
 *	This file is part of Astrohelion
 *
 *  Astrohelion is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Astrohelion is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Astrohelion.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <exception>
#include <stdexcept>
#include <string>

namespace astrohelion{
/**
 *	@brief Base Exception class for all my exceptions
 *
 *	By catching this base class, I catch all custom exceptions
 */
class Exception : public std::exception{
	public:
		Exception();
		Exception(const char*);
		Exception(const Exception&);
		virtual ~Exception();
		virtual const char* what() const throw();

	protected:
		std::string msg;	//!< Custom error message
};

/**
 *	@brief Exception class for matrices, fired when the dimension of matrices
 *	makes a function impossible.
 *
 * 	Size mismatch means two matrices cannot be added, subtracted, or multiplied 
 *  because they do not have the appropriate dimensions
 */
class SizeMismatchException : public Exception, std::runtime_error{
	public:
		SizeMismatchException();
		SizeMismatchException(const char*);
		// const char* what() const throw();
};

/**
 *	@brief Exception class for numerical methods that attempt to converge on a
 *	solution.
 *
 *	This exception should be thrown when the solution fails to converge (or diverges)
 *	and provides a way for the calling function to gracefully handle the divergence
 */
class DivergeException : public Exception, std::runtime_error{
	public:
		DivergeException();
		DivergeException(const char*);
		// const char* what() const throw();
};

/** 
 *	@brief Exception class for linear algebra errors
 *
 *	This exception should be thrown for cases like trying to factor singular
 *	matrices, taking the deteriminant of a non-square matrix, etc.
 */
class LinAlgException : public Exception, std::runtime_error{
	public:
		LinAlgException();
		LinAlgException(const char*);
		// const char* what() const throw();
};

}// END of Astrohelion namespace
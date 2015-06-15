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

#ifndef __H_TPAT_EXCEPTIONS_
#define __H_TPAT_EXCEPTIONS_

#include <exception>
#include <stdexcept>

/**
 *	@brief Base Exception class for all my exceptions
 *
 *	By catching this base class, I catch all custom exceptions
 */
class tpat_exception : public std::exception{
	public:
		/** @brief Default constructor */
		tpat_exception() : msg("Custom exception!") {}
		
		/**
		 *	@brief Create an exception with a custom message
		 *	@param m a message
		 */
		tpat_exception(const char* m) : msg(m) {}

		/** @brief describe the exception */
		virtual const char* what() const throw(){
			return msg;
		}
	protected:
		const char *msg;	//<! Custom error message
};

/**
 *	@brief Exception class for matrices, fired when the dimension of matrices
 *	makes a function impossible.
 *
 * 	Size mismatch means two matrices cannot be added, subtracted, or multiplied 
 *  because they do not have the appropriate dimensions
 */
class tpat_sizeMismatch : public tpat_exception, std::runtime_error{
	public:
		/** @brief Default constructor */
		tpat_sizeMismatch() : tpat_exception("Matrix dimensions do not match!"), std::runtime_error("Matrix dimensions do not match!"){}
		
		/**
		 *	@brief Create an exception with a custom message
		 *	@param m a message
		 */
		tpat_sizeMismatch(const char* m) : tpat_exception(m), std::runtime_error(m) {}
		
		/** @brief describe the exception */
		const char* what() const throw(){
	    	return msg;
	  	}
};

/**
 *	@brief Exception class for numerical methods that attempt to converge on a
 *	solution.
 *
 *	This exception should be thrown when the solution fails to converge (or diverges)
 *	and provides a way for the calling function to gracefully handle the divergence
 */
class tpat_diverge : public tpat_exception, std::runtime_error{
	public:
		/** Default constructor */
		tpat_diverge() : tpat_exception("Did not converge!"), std::runtime_error("Did not converge!"){}

		/**
		 *	@brief Create an exception with a custom message
		 *	@param m a message
		 */
		tpat_diverge(const char* m) : tpat_exception(m), std::runtime_error(m) {}
		
		/** @brief describe the exception */
		const char* what() const throw(){
			return msg;
		}
};

/** 
 *	@brief Exception class for linear algebra errors
 *
 *	This exception should be thrown for cases like trying to factor singular
 *	matrices, taking the deteriminant of a non-square matrix, etc.
 */
class tpat_linalg_err : public tpat_exception, std::runtime_error{
	public:
		/** Default constructor */
		tpat_linalg_err() : tpat_exception("Linear algebra error!"), std::runtime_error("Linear algebra error!"){}

		/**
		 *	@brief Create an exception with a custom message
		 *	@param m a message
		 */
		tpat_linalg_err(const char* m) : tpat_exception(m), std::runtime_error(m) {}
		
		/** @brief describe the exception */
		const char* what() const throw(){
			return msg;
		}
};

#endif
/**
 *  @file Exceptions.cpp
 *	@brief Contains member functions for various custom exception classes
 *	
 *	@author Andrew Cox
 *	@version May 25, 2016
 *	@copyright GNU GPL v3.0
 */
#include "Exceptions.hpp"

#include "AsciiOutput.hpp"
#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

namespace astrohelion{
//-----------------------------------------------------------------------------
//  ** Exception Functions **
//-----------------------------------------------------------------------------

/**
 *  @brief Default constructor
 */
Exception::Exception() : msg("Custom exception!") {}

/**
 *	@brief Create an exception with a custom message
 *	@param m a message
 */
Exception::Exception(const char* m) : msg(m) {}

/**
 *  @brief Copy constructor
 * 
 *  @param e reference to another exception object
 */
Exception::Exception(const Exception &e) : msg(e.msg){}

/**
 *  @brief Destructor
 */
Exception::~Exception(){}

/** @brief describe the exception */
const char* Exception::what() const throw(){
	// void *array[10];
	// size_t size = backtrace(array, 10);
	// backtrace_symbols_fd(array, size, STDERR_FILENO);

	void *callstack[128];
	int i, frames = backtrace(callstack, 128);
	char** strs = backtrace_symbols(callstack, frames);
	printf(RED);
	for(i = 0; i < frames; ++i){
		printf("%s\n", strs[i]);
	}
	printf(RESET);
	free(strs);
	return msg.c_str();
}//====================================================

//-----------------------------------------------------------------------------
//  ** SizeMismatch Functions **
//-----------------------------------------------------------------------------

/** @brief Default constructor */
SizeMismatchException::SizeMismatchException() : Exception("Matrix dimensions do not match!"), std::runtime_error("Matrix dimensions do not match!"){}

/**
 *	@brief Create an exception with a custom message
 *	@param m a message
 */
SizeMismatchException::SizeMismatchException(const char* m) : Exception(m), std::runtime_error(m) {}


//-----------------------------------------------------------------------------
//  ** Diverge Functions **
//-----------------------------------------------------------------------------

/** Default constructor */
DivergeException::DivergeException() : Exception("Did not converge!"), std::runtime_error("Did not converge!"){}

/**
 *	@brief Create an exception with a custom message
 *	@param m a message
 */
DivergeException::DivergeException(const char* m) : Exception(m), std::runtime_error(m) {}


//-----------------------------------------------------------------------------
//  ** LinAlg Functions **
//-----------------------------------------------------------------------------

/** Default constructor */
LinAlgException::LinAlgException() : Exception("Linear algebra error!"), std::runtime_error("Linear algebra error!"){}

/**
 *	@brief Create an exception with a custom message
 *	@param m a message
 */
LinAlgException::LinAlgException(const char* m) : Exception(m), std::runtime_error(m) {}



}// END of Astrohelion namespace
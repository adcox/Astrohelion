#include "tpat_exceptions.hpp"

#include "tpat_ascii_output.hpp"
#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

//-----------------------------------------------------------------------------
//  ** TPAT_EXCEPTION Functions **
//-----------------------------------------------------------------------------

/**
 *  @brief Default constructor
 */
TPAT_Exception::TPAT_Exception() : msg("Custom exception!") {}

/**
 *	@brief Create an exception with a custom message
 *	@param m a message
 */
TPAT_Exception::TPAT_Exception(const char* m) : msg(m) {}

/**
 *  @brief Copy constructor
 * 
 *  @param e reference to another exception object
 */
TPAT_Exception::TPAT_Exception(const TPAT_Exception &e) : msg(e.msg){}

/** @brief describe the exception */
const char* TPAT_Exception::what() const throw(){
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
}

//-----------------------------------------------------------------------------
//  ** TPAT_SIZEMISTPAT_Constraint_Tp::MATCH_ Functions **
//-----------------------------------------------------------------------------

/** @brief Default constructor */
TPAT_SizeMismatch::TPAT_SizeMismatch() : TPAT_Exception("Matrix dimensions do not match!"), std::runtime_error("Matrix dimensions do not match!"){}

/**
 *	@brief Create an exception with a custom message
 *	@param m a message
 */
TPAT_SizeMismatch::TPAT_SizeMismatch(const char* m) : TPAT_Exception(m), std::runtime_error(m) {}

/** @brief describe the exception */
const char* TPAT_SizeMismatch::what() const throw(){
	return msg.c_str();
}


//-----------------------------------------------------------------------------
//  ** TPAT_DIVERGE Functions **
//-----------------------------------------------------------------------------

/** Default constructor */
TPAT_Diverge::TPAT_Diverge() : TPAT_Exception("Did not converge!"), std::runtime_error("Did not converge!"){}

/**
 *	@brief Create an exception with a custom message
 *	@param m a message
 */
TPAT_Diverge::TPAT_Diverge(const char* m) : TPAT_Exception(m), std::runtime_error(m) {}

/** @brief describe the exception */
const char* TPAT_Diverge::what() const throw(){
	return msg.c_str();
}

//-----------------------------------------------------------------------------
//  ** TPAT_LINALG_ERR Functions **
//-----------------------------------------------------------------------------

/** Default constructor */
TPAT_LinAlg_Err::TPAT_LinAlg_Err() : TPAT_Exception("Linear algebra error!"), std::runtime_error("Linear algebra error!"){}

/**
 *	@brief Create an exception with a custom message
 *	@param m a message
 */
TPAT_LinAlg_Err::TPAT_LinAlg_Err(const char* m) : TPAT_Exception(m), std::runtime_error(m) {}

/** @brief describe the exception */
const char* TPAT_LinAlg_Err::what() const throw(){
	return msg.c_str();
}
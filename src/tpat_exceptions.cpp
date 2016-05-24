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
tpat_exception::tpat_exception() : msg("Custom exception!") {}

/**
 *	@brief Create an exception with a custom message
 *	@param m a message
 */
tpat_exception::tpat_exception(const char* m) : msg(m) {}

/**
 *  @brief Copy constructor
 * 
 *  @param e reference to another exception object
 */
tpat_exception::tpat_exception(const tpat_exception &e) : msg(e.msg){}

/** @brief describe the exception */
const char* tpat_exception::what() const throw(){
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
//  ** TPAT_SIZEMISMATCH Functions **
//-----------------------------------------------------------------------------

/** @brief Default constructor */
tpat_sizeMismatch::tpat_sizeMismatch() : tpat_exception("Matrix dimensions do not match!"), std::runtime_error("Matrix dimensions do not match!"){}

/**
 *	@brief Create an exception with a custom message
 *	@param m a message
 */
tpat_sizeMismatch::tpat_sizeMismatch(const char* m) : tpat_exception(m), std::runtime_error(m) {}

/** @brief describe the exception */
const char* tpat_sizeMismatch::what() const throw(){
	return msg.c_str();
}


//-----------------------------------------------------------------------------
//  ** TPAT_DIVERGE Functions **
//-----------------------------------------------------------------------------

/** Default constructor */
tpat_diverge::tpat_diverge() : tpat_exception("Did not converge!"), std::runtime_error("Did not converge!"){}

/**
 *	@brief Create an exception with a custom message
 *	@param m a message
 */
tpat_diverge::tpat_diverge(const char* m) : tpat_exception(m), std::runtime_error(m) {}

/** @brief describe the exception */
const char* tpat_diverge::what() const throw(){
	return msg.c_str();
}

//-----------------------------------------------------------------------------
//  ** TPAT_LINALG_ERR Functions **
//-----------------------------------------------------------------------------

/** Default constructor */
tpat_linalg_err::tpat_linalg_err() : tpat_exception("Linear algebra error!"), std::runtime_error("Linear algebra error!"){}

/**
 *	@brief Create an exception with a custom message
 *	@param m a message
 */
tpat_linalg_err::tpat_linalg_err(const char* m) : tpat_exception(m), std::runtime_error(m) {}

/** @brief describe the exception */
const char* tpat_linalg_err::what() const throw(){
	return msg.c_str();
}
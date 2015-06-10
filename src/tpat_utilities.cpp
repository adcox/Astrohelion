
#include "tpat_utilities.hpp"

#include "tpat_ascii_output.hpp"

#include "matio.h"

#include <complex>
#include <iostream>
#include <string>

using namespace std;

std::string complexToStr(std::complex<double> num){
    char buffer[64];
    sprintf(buffer, "%.4e%s%.4ej", real(num), imag(num) < 0 ? " - " : " + ", abs(imag(num)));
    return string(buffer);
}

/**
 *  @brief A wrapper function to print a message
 *	@param verbose whether or not to be verbose; message is not printed if verbose is false
 *	@param format a standard format string literal to pass to <tt>vprintf</tt>
 */
void printVerb(bool verbose, const char * format, ...){
    if(verbose){
        va_list args;
        va_start(args, format);
        vprintf(format, args);
        va_end(args);
    }
}//==========================================

/**
 *	@brief Print an error message to the standard output in red
 *	@param format a standard format string literal to pass to <tt>vprintf</tt>
 */
void printErr(const char * format, ...){
	printf(RED);
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
    printf(RESET);
}//==========================================

/**
 *	@brief Print a warning messate to the standard output in yellow
 *	@param format a standard format string literal to pass to <tt>vprintf</tt>
 */
void printWarn(const char * format, ...){
    printf(YELLOW);
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
    printf(RESET);
}//==========================================

/**
 *	@brief Print a message to the standard output using an ASCII escape-type color
 *	@param color one of the constant color values stored in the ascii-output header
 *	@param format a standard format string literal to pass to <tt>vprintf</tt>
 */
void printColor(const char* color, const char * format, ...){
    printf("%s", color);
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
    printf(RESET);
}//==========================================

/**
 *	@brief Print a message to the standard output using an ASCII escape-type color
 *	@param verbose whether or not to print the string
 *	@param color one of the constant color values stored in the ascii-output header
 *	@param format a standard format string literal to pass to <tt>vprintf</tt>
 */
void printVerbColor(bool verbose, const char* color, const char * format, ...){
	if(verbose){
	    printf("%s", color);
	    va_list args;
	    va_start(args, format);
	    vprintf(format, args);
	    va_end(args);
	    printf(RESET);
	}
}//==========================================

/**
 *  @brief Save a variable to a .mat file, performing error checks along the way. 
 *
 *  Once the variable is written to file, it is freed from memory
 *
 *  @param matFile a pointer to the matlab output file
 *  @param matvar a pointer to the matlab variable object
 *  @param varName a literal string that describes the variable; only used in error message output
 *  @param comp an enum that describes the compression strategy. Options are:
 *      MAT_COMPRESSION_NONE - no compression
 *      MAT_COMPRESSION_ZLIB - zlib compression
 */
void saveVar(mat_t *matFile, matvar_t *matvar, const char* varName, matio_compression comp){
    if(NULL == matvar){
        printErr("Error creating variable for '%s'\n", varName);
    }else{
        Mat_VarWrite(matFile, matvar, comp);
        Mat_VarFree(matvar);
    }
}//==============================

/**
 *  @brief Suspend operation until the user presses a key
 */
void waitForUser(){
    cout << "Press ENTER to continue...";
    cin.ignore();
    cout<<endl;
}
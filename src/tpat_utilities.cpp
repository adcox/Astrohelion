/**
 *  @file tpat_utilities.cpp
 */
/*
 *  Trajectory Propagation and Analysis Toolkit 
 *  Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of the Trajectory Propagation and Analysis Toolkit (TPAT).
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
#include "tpat.hpp"

#include "tpat_utilities.hpp"

#include "tpat_ascii_output.hpp"

#include "matio.h"

#include <complex>
#include <string>

/**
 *  @brief Turn a complex number into a string, e.g. 1.2345 + 0.9876j
 *  @param num a complex number
 *  @return the complex number as a string
 */
std::string complexToStr(std::complex<double> num){
    char buffer[64];
    sprintf(buffer, "%.4e%s%.4ej", real(num), imag(num) < 0 ? " - " : " + ", std::abs(imag(num)));
    return std::string(buffer);
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
 *  @brief Read a double from a mat file
 *  @param matFile a pointer to the matlab file in quesiton
 *  @param varName the name of the variable in the mat file
 *  @return the value of the variable
 *  @throws tpat_exception if there is trouble reading or parsing the variable
 */
double readDoubleFromMat(mat_t *matFile, const char* varName){
    double result = 0;
    
    matvar_t *matvar = Mat_VarRead(matFile, varName);
    if(matvar == NULL){
        throw tpat_exception("tpat_utilities::readDoubleFromMat: Could not read variable from file");
    }else{
        if(matvar->class_type == MAT_C_DOUBLE && matvar->data_type == MAT_T_DOUBLE){
            double *data = static_cast<double *>(matvar->data);

            if(data != NULL){
                result = *data;
            }else{
                throw tpat_exception("tpat_utilities::readDoubleFromMat: No data");
            }
        }else{
            throw tpat_exception("tpat_utilities::readDoubleFromMat: Incompatible data file: unsupported data type/class");
        }
    }

    Mat_VarFree(matvar);
    return result;
}//==============================================

/**
 *  @brief Read an int from a mat file
 *  @param matFile a pointer to the matlab file in quesiton
 *  @param varName the name of the variable in the mat file
 *  @param aType the expected variable type
 *  @param aClass the expected variable class
 *  @return the value of the variable
 *  @throws tpat_exception if there is trouble reading or parsing the variable
 */
int readIntFromMat(mat_t *matFile, const char* varName, matio_types aType,
    matio_classes aClass){

    int result = 0;
    matvar_t *matvar = Mat_VarRead(matFile, varName);
    if(matvar == NULL){
        throw tpat_exception("tpat_utilities::readIntFromMat: Could not read variable from file");
    }else{
        if(matvar->class_type == aClass && matvar->data_type == aType){
            int *data = static_cast<int *>(matvar->data);

            if(data != NULL){
                result = *data;
            }else{
                throw tpat_exception("tpat_utilities::readIntFromMat: No data");
            }
        }else{
            throw tpat_exception("tpat_utilities::readIntFromMat: Incompatible data file: unsupported data type/class");
        }
    }

    Mat_VarFree(matvar);
    return result;
}//==============================================

/**
 *  @brief Read a string from a matlab file
 *
 *  @param matFile a pointer to the matlab file in questions
 *  @param varName the name of the variable in the mat file
 *  @param aType the expected variable type (e.g. MAT_T_UINT8)
 *  @param aClass the expected variable class (e.g. MAT_C_CHAR)
 *  @return the string
 *  @throws tpat_exception if there is trouble reading or parsing the variable.
 */
std::string readStringFromMat(mat_t *matFile, const char* varName, matio_types aType,
    matio_classes aClass){

    std::string result = "";
    matvar_t *matvar = Mat_VarRead(matFile, varName);
    if(matvar == NULL){
        throw tpat_exception("tpat_utilities::readStringFromMat: Could not read variable from file");
    }else{
        if(matvar->class_type == aClass && matvar->data_type == aType){
            char *data = static_cast<char *>(matvar->data);

            if(data != NULL){
                result = std::string(data, matvar->dims[1]);
            }else{
                throw tpat_exception("tpat_utilities::readStringFromMat: No data");
            }
        }else{
            throw tpat_exception("tpat_utilities::readStringFromMat: Incompatible data file: unsupported data type/class");
        }
    }

    Mat_VarFree(matvar);
    return result;
}//=============================================

/**
 *  @brief Suspend operation until the user presses a key
 */
void waitForUser(){
    std::cout << "Press ENTER to continue...";
    std::cin.ignore();
    std::cout << std::endl;
}//========================================================




//
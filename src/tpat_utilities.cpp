/**
 *  @file tpat_utilities.cpp
 *  @brief Similar to tpat_calculations, but these functions are
 *  more general, less math-y
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
    sprintf(buffer, "%.4e%s%.4ej", std::real(num), std::imag(num) < 0 ? " - " : " + ", std::abs(std::imag(num)));
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
 *  @brief Prints out a set of eigenvalues and eigenvectors
 *  @param eigData a set of vectors that contain data about an
 *  eigensystem; such sets are returned from tpat_matrix::eig()
 *  @see tpat_matrix::eig()
 */
void printEigenData(std::vector< std::vector<cdouble> > eigData){
    std::vector<cdouble> vals = eigData[0];
    std::vector<cdouble> vecs = eigData[1];

    printf("Eigenvalues:\n");
    for(size_t v = 0; v < vals.size(); v++){
        printf("%27s", complexToStr(vals[v]).c_str());
    }
    printf("\nEigenvectors:\n");
    for(size_t j = 0; j < vals.size(); j++){
        for(size_t v = 0; v < vals.size(); v++){
            printf("%27s", complexToStr(vecs[v*vals.size() + j]).c_str());
        }
        printf("\n");
    }
    printf("\n");
}//============================================

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
}//=========================================================

/**
 * @brief Save a matrix of data to a Matlab .mat file
 * @details This function will save the transpose of the input matrix, so be 
 * sure to transpose after loading into Matlab
 * 
 * @param filename name/path of the file
 * @param varName variable name
 * @param data vector of data
 * @param rows number of rows in the matrix
 * @param cols number of columns in the matrix
 */
void saveMatrixToFile(const char* filename, const char* varName, std::vector<double> data, size_t rows, size_t cols){
    if(data.size() < rows*cols)
        throw tpat_exception("tpat_utilities::saveMatrixToFile: Input data has fewer elements than specified by the rows and cols arguments");

    mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
    if(NULL != matfp){
        size_t dims[2] = {cols, rows};
        matvar_t *matvar = Mat_VarCreate(varName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(data[0]), MAT_F_DONT_COPY_DATA);
        saveVar(matfp, matvar, varName, MAT_COMPRESSION_NONE);
    }else{
        printErr("tpat_utilities::saveMatrixToFile: Error creating mat file\n");
    }
    Mat_Close(matfp);
}//========================================================

/**
 *  @brief Read a matrix of doubles from a .mat file
 * 
 *  @param filename relative or absolute path to the .mat data file
 *  @param varName the name of the variable/matrix to read from the file
 * 
 *  @return a column-major-order vector containing the data from the desired matrix
 *  @throws tpat_exception if the file cannot be opened, if the variable doesn't exist,
 *  or if the variable contains something other than doubles
 */
std::vector<double> readMatrixFromMat(const char *filename, const char *varName){
 
    mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
    if(matfp ==  NULL)
        throw tpat_exception("tpat_utilities::readMatrixFromFile: Could not open data file\n");

    matvar_t *matvar = Mat_VarRead(matfp, varName);
    if(matvar == NULL){
        throw tpat_exception("tpat_utilities::readMatrixFromFile: Could not read variable data");
    }else{

        int dataSize = (matvar->dims[0])*(matvar->dims[1]);

        if(matvar->class_type == MAT_C_DOUBLE && matvar->data_type == MAT_T_DOUBLE){
            double *data = static_cast<double *>(matvar->data);

            if(data != NULL){
                std::vector<double> vecData(data, data+dataSize);
                
                Mat_VarFree(matvar);
                return vecData;
            }
        }else{
            Mat_VarFree(matvar);
            throw tpat_exception("tpat_utilities::readMatrixFromFile: Incompatible data file: unsupported data type/class");
        }
    }
}//==========================================================

/**
 *  @brief Read a double from a mat file
 *  @param matFile a pointer to the matlab file in quesiton
 *  @param varName the name of the variable in the mat file
 *  @return the value of the variable
 *  @throws tpat_exception if there is trouble reading or parsing the variable
 */
double readDoubleFromMat(mat_t *matFile, const char* varName){
    double result = 2;
    
    matvar_t *matvar = Mat_VarRead(matFile, varName);
    if(matvar == NULL){
        throw tpat_exception("tpat_utilities::readDoubleFromMat: Could not read variable from file");
    }else{
        switch(matvar->data_type){
            case MAT_T_UINT8:
            case MAT_T_UINT16:
            case MAT_T_UINT32:
            case MAT_T_UINT64:
            {
                uint *data = static_cast<uint *>(matvar->data);
                if(data != NULL)
                    result = double(*data);
                else
                    throw tpat_exception("tpat_utilities::readDoubleFromMat: No data");
                break;
            }
            case MAT_T_INT8:
            case MAT_T_INT16:
            case MAT_T_INT32:
            case MAT_T_INT64:
            {
                int *data = static_cast<int *>(matvar->data);
                if(data != NULL)
                    result = double(*data);
                else
                    throw tpat_exception("tpat_utilities::readDoubleFromMat: No data");
                break;
            }
            case MAT_T_SINGLE:
            case MAT_T_DOUBLE:
            {
                double *data = static_cast<double *>(matvar->data);
                if(data != NULL)
                    result = double(*data);
                else
                    throw tpat_exception("tpat_utilities::readDoubleFromMat: No data");
                break;
            }
            default:
                throw tpat_exception("tpat_utilities::readDoubleFromMat: incompatible data-type");
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
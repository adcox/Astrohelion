/**
 *  @file Utilities.cpp
 *  @brief Similar to tpat_calculations, but these functions are
 *  more general, less math-y
 *  
 *  @author Andrew Cox
 *  @version May 25, 2016
 *  @copyright GNU GPL v3.0
 */
/*
 *  Astrohelion 
 *  Copyright 2015, Andrew Cox; Protected under the GNU GPL v3.0
 *  
 *  This file is part of Astrohelion
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

#include "AsciiOutput.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"

#include "cspice/SpiceZfc.h"    // prototypes for functions
#include "matio.h"

#include <complex>
#include <fstream>
#include <iostream>
#include <string>

namespace astrohelion{
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
 *  @brief Save a single double-precision value to a Matlab file
 * 
 *  @param matfp a pointer to the matlab output file
 *  @param varName the name of the variable
 *  @param data data value
 */
void saveDoubleToFile(mat_t *matfp, const char *varName, double data){
    if(NULL != matfp){
        size_t dims[2] = {1, 1};
        matvar_t *matvar = Mat_VarCreate(varName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &data, MAT_F_DONT_COPY_DATA);
        saveVar(matfp, matvar, varName, MAT_COMPRESSION_NONE);
    }else{
        printErr("Utilities::saveDoubleToFile: Error creating mat file\n");
    }
}//========================================================

/**
 * @brief Save a matrix of data to a Matlab .mat file
 * @details The data passed in the <tt>data</tt> vector is transposed 
 *  and saved to the specified Matlab file
 * 
 * @param filename name/path of the file
 * @param varName variable name
 * @param data vector of data in row-major order
 * @param rows number of rows in the matrix
 * @param cols number of columns in the matrix
 */
void saveMatrixToFile(const char* filename, const char* varName, std::vector<double> data, size_t rows, size_t cols){
    mat_t *matfp = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
    saveMatrixToFile(matfp, varName, data, rows, cols);
    Mat_Close(matfp);
}//========================================================

/**
 *  @brief Save a matrix of data to an open matlab .mat file
 *  @details The data passed in the <tt>data</tt> vector is transposed 
 *  and saved to the specified Matlab file 
 * 
 *  @param matfp An open Matlab .mat file
 *  @param varName name of the variable within the .mat file
 *  @param data a vector of data in row-major order
 *  @param rows number of rows in the matrix
 *  @param cols number of columns in the matrix
 *  @throws Exception if <tt>data</tt> does not have enough elements
 *  to construct a matrix with the specified number of rows and columns
 */
void saveMatrixToFile(mat_t *matfp, const char *varName, std::vector<double> data, size_t rows, size_t cols){
    if(data.size() < rows*cols)
        throw Exception("Utilities::saveMatrixToFile: Input data has fewer elements than specified by the rows and cols arguments");

    if(NULL != matfp){
        size_t dims[2] = {rows, cols};

        std::vector<double> data_trans(data.size());
        for(size_t r = 0; r < rows; r++){
            for(size_t c = 0; c < cols; c++){
                data_trans[c*rows + r] = data[r*cols + c];
            }
        }

        matvar_t *matvar = Mat_VarCreate(varName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(data_trans[0]), MAT_F_DONT_COPY_DATA);
        saveVar(matfp, matvar, varName, MAT_COMPRESSION_NONE);
    }else{
        printErr("Utilities::saveMatrixToFile: Error creating mat file\n");
    }
}//=========================================================

/**
 *  @brief Read a matrix of doubles from a .mat file
 * 
 *  @param filename relative or absolute path to the .mat data file
 *  @param varName the name of the variable/matrix to read from the file
 * 
 *  @return a column-major-order vector containing the data from the desired matrix
 *  @throws Exception if the file cannot be opened, if the variable doesn't exist,
 *  or if the variable contains something other than doubles
 *  @throws Exception if the data file cannot be opened or the variable cannot be read
 */
MatrixXRd readMatrixFromMat(const char *filename, const char *varName){
 
    mat_t *matfp = Mat_Open(filename, MAT_ACC_RDONLY);
    if(matfp ==  NULL)
        throw Exception("Utilities::readMatrixFromFile: Could not open data file\n");

    // For debugging, to print a list of all variables in the file
    // matvar_t *tempvar;
    // while( (tempvar = Mat_VarReadNextInfo(matfp)) != NULL){
    //     printf("%s\n", tempvar->name);
    //     Mat_VarFree(tempvar);
    //     tempvar = NULL;
    // }

    matvar_t *matvar = Mat_VarRead(matfp, varName);
    if(matvar == NULL){
        Mat_Close(matfp);
        throw Exception("Utilities::readMatrixFromFile: Could not read variable data");
    }else{

        if(matvar->class_type == MAT_C_DOUBLE && matvar->data_type == MAT_T_DOUBLE){
            double *data = static_cast<double *>(matvar->data);

            MatrixXRd mat;
            
            if(data != NULL){
                mat = Eigen::Map<MatrixXRd>(data, matvar->dims[1], matvar->dims[0]);
                
                Mat_VarFree(matvar);
                return mat;
            }else{
                Mat_Close(matfp);
                return mat;
            }

            Mat_VarFree(matvar);
            Mat_Close(matfp);
            return mat;
        }else{
            Mat_VarFree(matvar);
            Mat_Close(matfp);
            throw Exception("Utilities::readMatrixFromFile: Incompatible data file: unsupported data type/class");
        }
    }
}//==========================================================

/**
 *  @brief Read a double from a mat file
 *  @param matFile a pointer to the matlab file in quesiton
 *  @param varName the name of the variable in the mat file
 *  @return the value of the variable
 *  @throws Exception if there is trouble reading or parsing the variable
 */
double readDoubleFromMat(mat_t *matFile, const char* varName){
    double result = 2;
    
    matvar_t *matvar = Mat_VarRead(matFile, varName);
    if(matvar == NULL){
        char msg[256];
        sprintf(msg, "Utilities::readDoubleFromMat: Could not read %s from file", varName);
        throw Exception(msg);
    }else{
        switch(matvar->data_type){
            case MAT_T_UINT8:
            case MAT_T_UINT16:
            case MAT_T_UINT32:
            case MAT_T_UINT64:
            {
                unsigned int *data = static_cast<unsigned int *>(matvar->data);
                if(data != NULL)
                    result = double(*data);
                else
                    throw Exception("Utilities::readDoubleFromMat: No data");
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
                    throw Exception("Utilities::readDoubleFromMat: No data");
                break;
            }
            case MAT_T_SINGLE:
            case MAT_T_DOUBLE:
            {
                double *data = static_cast<double *>(matvar->data);
                if(data != NULL)
                    result = double(*data);
                else
                    throw Exception("Utilities::readDoubleFromMat: No data");
                break;
            }
            default:
                throw Exception("Utilities::readDoubleFromMat: incompatible data-type");
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
 *  @throws Exception if there is trouble reading or parsing the variable.
 */
std::string readStringFromMat(mat_t *matFile, const char* varName, matio_types aType,
    matio_classes aClass){

    std::string result = "";
    matvar_t *matvar = Mat_VarRead(matFile, varName);
    if(matvar == NULL){
        throw Exception("Utilities::readStringFromMat: Could not read variable from file");
    }else{
        if(matvar->class_type == aClass && matvar->data_type == aType){
            char *data = static_cast<char *>(matvar->data);

            if(data != NULL){
                result = std::string(data, matvar->dims[1]);
            }else{
                throw Exception("Utilities::readStringFromMat: No data");
            }
        }else{
            throw Exception("Utilities::readStringFromMat: Incompatible data file: unsupported data type/class");
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

/**
 *  @brief Retrieve a string representing the name of a spacecraft or celestial body
 *  from its SPICE ID Code
 *  @details For a list of SPICE IDs, see the 
 *  <a href="http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/naif_ids.html">NAIF ID Page</a>
 *  
 *  @param ID Integer ID Code
 *  @return A string representation of the body or satellite
 */
std::string getNameFromSpiceID(int ID){
    SpiceInt lenout = 128;  // maximum length of output name
    SpiceInt code = static_cast<SpiceInt>(ID);
    SpiceChar name[128];
    SpiceBoolean found = false;

    bodc2n_c(code, lenout, name, &found);

    if(!found)
        return "BODY NOT FOUND";
    else
        return std::string(name);
}//============================================

/**
 *  @brief Get a SPICE ID from a body name
 *  @details Body names are case insensitive and leading/trailing spaces are not
 *  important. However, when the name is made up of more than word, they must be
 *  separated by at least one space.
 * 
 *  @param name body name to be translated to a SPICE ID code
 *  @return a SPICE integer ID code
 */
SpiceInt getSpiceIDFromName(const char *name){
    ConstSpiceChar *name_spice = static_cast<ConstSpiceChar*>(name);
    SpiceInt code = 0;
    SpiceBoolean found = false;

    bodn2c_c(name_spice, &code, &found);

    if(!found)
        throw Exception("Utilities::getSpiceIDFromName: Could not find a body name from that SPICE ID");
    else
        return code;
}//============================================

/**
 *  @brief Check to see if SPICE failed; if it did, throw a Exception with
 *  a custom error message
 * 
 *  @param customMsg A message for the Exception object
 *  @throws Exception if an error occured
 */
void checkAndReThrowSpiceErr(const char* customMsg){
    if(failed_c()){
        char errMsg[26];
        getmsg_c("short", 25, errMsg);
        printErr("Spice Error: %s\n", errMsg);
        reset_c();  // reset error status
        throw Exception(customMsg);
    }
}//============================================

/**
 *  @brief Save a matrix as a CSV file to be read by Excel or Matlab
 * 
 *  @param m the matrix
 *  @param filename Filename of the csv file (include the .csv extension!)
 */
void toCSV(MatrixXRd m, const char* filename){
    std::ofstream outFile(filename, std::ios::out);
    
    for (int r = 0; r < m.rows(); r++){
        for (int c = 0; c < m.cols(); c++){
            char buffer[64] = { };
            if(c < m.cols()-1)
                sprintf(buffer, "%.14f, ", m(r,c));
            else
                sprintf(buffer, "%.14f\n", m(r,c));

            outFile << buffer;
        }
    }

    outFile.close();
}//=============================================



} // End of astrohelion namespace
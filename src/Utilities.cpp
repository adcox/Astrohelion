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
 *  Copyright 2015-2017, Andrew Cox; Protected under the GNU GPL v3.0
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

#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>

#include "cspice/SpiceZfc.h"    // prototypes for functions
#include "matio.h"

#include "AsciiOutput.hpp"
#include "Exceptions.hpp"
#include "Utilities.hpp"

namespace astrohelion{

/**
 * @addtogroup util
 * \{
 */

/**
 *  @brief Retrieve the current CPU time in seconds past some reference
 *  @return the current CPU time in seconds past some reference
 */
double getCPUTime(){ return static_cast<double>(clock()) / CLOCKS_PER_SEC; }

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
 *	@param format a standard format string literal to pass to `vprintf`
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
 *	@param format a standard format string literal to pass to `vprintf`
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
 *	@param format a standard format string literal to pass to `vprintf`
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
 *	@param format a standard format string literal to pass to `vprintf`
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
 *	@param format a standard format string literal to pass to `vprintf`
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
 *  @details Once the variable is written to file, it is freed from memory, thus
 *  any further calls to Mat_VarFree on the matvar_t object will result in malloc errors.
 *
 *  @param matFile a pointer to the matlab output file
 *  @param matvar a pointer to the matlab variable object
 *  @param varName a literal string that describes the variable; only used in error message output
 *  @param comp an enum that describes the compression strategy. Options are:
 *      MAT_COMPRESSION_NONE - no compression
 *      MAT_COMPRESSION_ZLIB - zlib compression
 */
void saveVar(mat_t *matFile, matvar_t *matvar, const char* varName, matio_compression comp){
    if(nullptr == matvar){
        printErr("Error creating variable for '%s'\n", varName);
    }else{
        Mat_VarWrite(matFile, matvar, comp);
    }
    Mat_VarFree(matvar);
}//====================================================

/**
 *  @brief Save a single double-precision value to a Matlab file
 * 
 *  @param matfp a pointer to the matlab output file
 *  @param varName the name of the variable
 *  @param data data value
 */
void saveDoubleToFile(mat_t *matfp, const char *varName, double data){
    if(nullptr != matfp){
        size_t dims[2] = {1, 1};
        matvar_t *matvar = Mat_VarCreate(varName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &data, MAT_F_DONT_COPY_DATA);
        saveVar(matfp, matvar, varName, MAT_COMPRESSION_NONE);
    }else{
        printErr("Utilities::saveDoubleToFile: Error creating mat file\n");
    }
}//====================================================

/**
 * @brief Save a matrix of data to a Matlab .mat file
 * @details The data passed in the `data` vector is transposed 
 *  and saved to the specified Matlab file
 * 
 * @param filename name/path of the file
 * @param varName variable name
 * @param data vector of data in row-major order
 * @param rows number of rows in the matrix
 * @param cols number of columns in the matrix
 */
void saveMatrixToFile(const char* filename, const char* varName, std::vector<double> data, size_t rows, size_t cols){
    mat_t *matfp = Mat_CreateVer(filename, nullptr, MAT_FT_DEFAULT);
    saveMatrixToFile(matfp, varName, data, rows, cols);
    Mat_Close(matfp);
}//====================================================

/**
 *  @brief Save a matrix of data to an open matlab .mat file
 *  @details The data passed in the `data` vector is transposed 
 *  and saved to the specified Matlab file 
 * 
 *  @param matfp An open Matlab .mat file
 *  @param varName name of the variable within the .mat file
 *  @param data a vector of data in row-major order
 *  @param rows number of rows in the matrix
 *  @param cols number of columns in the matrix
 *  @throws Exception if `data` does not have enough elements
 *  to construct a matrix with the specified number of rows and columns
 */
void saveMatrixToFile(mat_t *matfp, const char *varName, std::vector<double> data, size_t rows, size_t cols){
    if(data.size() < rows*cols)
        throw Exception("Utilities::saveMatrixToFile: Input data has fewer elements than specified by the rows and cols arguments");

    if(nullptr != matfp){
        size_t dims[2] = {rows, cols};

        std::vector<double> data_trans(data.size());
        for(unsigned int r = 0; r < rows; r++){
            for(unsigned int c = 0; c < cols; c++){
                data_trans[c*rows + r] = data[r*cols + c];
            }
        }

        matvar_t *matvar = Mat_VarCreate(varName, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &(data_trans[0]), MAT_F_DONT_COPY_DATA);
        saveVar(matfp, matvar, varName, MAT_COMPRESSION_NONE);
    }else{
        printErr("Utilities::saveMatrixToFile: Error creating mat file\n");
    }
}//====================================================

/**
 *  @brief Save a string to an open Matlab data file
 *  @details [long description]
 * 
 *  @param matfp pointer to Matlab file
 *  @param varName name of the variable within the Matlab file
 *  @param text value of the variable within the Matlab file
 *  @param strlen number of characters (including the final nullptr char) in `text`
 */
void saveStringToFile(mat_t *matfp, const char *varName, std::string text, const int strlen){
    char *ctxt = (char *)alloca(text.size()+1);
    memcpy(ctxt, text.c_str(), text.size()+1);
    size_t dims[] {1, text.size()+1};
    matvar_t *matvar = Mat_VarCreate(varName, MAT_C_CHAR, MAT_T_UINT8, 2, dims, ctxt, 0);
    astrohelion::saveVar(matfp, matvar, varName, MAT_COMPRESSION_NONE);

    (void) strlen;
    // char text_chars[strlen];
    // strcpy(text_chars, text.c_str());
    // size_t dims[] {1, text.length()};
    // matvar_t *matvar = Mat_VarCreate(varName, MAT_C_CHAR, MAT_T_UINT8, 2, dims, &(text_chars[0]), 0);
    // astrohelion::saveVar(matfp, matvar, varName, MAT_COMPRESSION_NONE);
}//====================================================

/**
 * @brief Save a timestamp to a file
 * @details The timestamp is saved as "hh:mm:ssUTC dd/mm/yyyy"
 * 
 * @param matfp pointer to Matlab file
 * @param varName name of the variable
 */
void saveTimestampToFile(mat_t *matfp, const char *varName){
    time_t timer;

    time(&timer);   // Get the current time
    struct tm now_utc = *gmtime(&timer);    // Convert current time to UTC
    char txt_now[22];
    sprintf(txt_now, "%02d:%02d:%02dUTC %02d/%02d/%4d", now_utc.tm_hour,
        now_utc.tm_min, now_utc.tm_sec, now_utc.tm_mday,
        now_utc.tm_mon + 1, now_utc.tm_year + 1900);

    size_t dims[] {1, 22};
    matvar_t *matvar = Mat_VarCreate(varName, MAT_C_CHAR, MAT_T_UINT8, 2, dims, txt_now, 0);
    saveVar(matfp, matvar, varName, MAT_COMPRESSION_NONE);
}//====================================================

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
    if(matfp ==  nullptr)
        throw Exception("Utilities::readMatrixFromFile: Could not open data file\n");

    // For debugging, to print a list of all variables in the file
    // matvar_t *tempvar;
    // while( (tempvar = Mat_VarReadNextInfo(matfp)) != nullptr){
    //     printf("%s\n", tempvar->name);
    //     Mat_VarFree(tempvar);
    //     tempvar = nullptr;
    // }

    matvar_t *matvar = Mat_VarRead(matfp, varName);
    if(matvar == nullptr){
        Mat_Close(matfp);
        throw Exception("Utilities::readMatrixFromFile: Could not read variable data");
    }else{

        if(matvar->class_type == MAT_C_DOUBLE && matvar->data_type == MAT_T_DOUBLE){
            double *data = static_cast<double *>(matvar->data);

            MatrixXRd mat;
            
            if(data != nullptr){
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
    if(matvar == nullptr){
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
                if(data != nullptr)
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
                if(data != nullptr)
                    result = double(*data);
                else
                    throw Exception("Utilities::readDoubleFromMat: No data");
                break;
            }
            case MAT_T_SINGLE:
            case MAT_T_DOUBLE:
            {
                double *data = static_cast<double *>(matvar->data);
                if(data != nullptr)
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
    if(matvar == nullptr){
        throw Exception("Utilities::readStringFromMat: Could not read variable from file");
    }else{
        if(matvar->class_type == aClass && matvar->data_type == aType){
            char *data = static_cast<char *>(matvar->data);

            if(data != nullptr){
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
    checkAndReThrowSpiceErr("getNameFromSpiceID error");
    
    if(!found){
        return "BODY NOT FOUND";
    }
    else{
        return std::string(name);
    }
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
    checkAndReThrowSpiceErr("getSpiceIDFromName error");

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
 *  @brief Resolve double angle ambiquity from inverse trig
 *  functions
 *  @details [long description]
 * 
 *  @param asinVal The result of arcsin(value)
 *  @param acosVal The result of arccos(value)
 * 
 *  @return The common value between the double angles returned by arcsin(value)
 *  and arccos(value)
 */
double resolveAngle(double asinVal, double acosVal){
    if(std::isnan(asinVal) || std::isnan(acosVal))
        throw Exception("resolveAngle: one of the inputs is NAN; cannot resolve angle");

    double asinVals[] = {asinVal, PI - asinVal};
    double acosVals[] = {acosVal, -acosVal};

    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            if(std::abs(sin(acosVals[i]) - sin(asinVals[j])) < 1e-8)
                return acosVals[i];
            // if(std::abs(acosVals[i] - asinVals[j]) < 1e-8)
            //     return acosVals[i];
        }
    }

    // Could not find a value
    printErr("No Common Angle:\n  arcsin vals = %.4f, %.4f\n", asinVals[0], asinVals[1]);
    printErr("  arccos vals = %.4f, %.4f\n", acosVals[0], acosVals[1]);
    throw Exception("resolveAngle: Could not resolve angle");
}//=============================================

/**
 *  @brief Determine the value of a number within some bounds
 *  @details This is particularly useful for arguments for 
 *  inverse trig functions: small numerical errors may result 
 *  in arccos(1 + eps) or arsin(-1 - eps) where eps is some small
 *  value. The resulting inverse trig function will be invalide
 *  and return NAN. In this case, use boundValue(value, -1, 1)
 *  to ensure that the argument inside the trig function is between
 *  -1 and 1.
 * 
 *  @param val the value to bound
 *  @param min minimum allowable value
 *  @param max maximum allowable value
 *  @return The value if it is within the bounds, or the min/max
 *  bound if val is outside the bounds.
 */
double boundValue(double val, double min, double max){
    if(val < min)
        return min;

    if(val > max)
        return max;

    return val;
}//====================================================

bool fileExists (const char *filename) {
    struct stat buffer;   
    return (stat (filename, &buffer) == 0); 
}//====================================================

/**
 *  @brief Convert and angle to its equivalent between -pi and pi
 * 
 *  @param val value in radians
 *  @return An equivalent value between -pi and pi
 */
double wrapToPi(double val){
    while(abs(val) > PI)
        val -= astrohelion::sign(val)*2*PI;

    return val;
}//====================================================

/**
 *  @brief Convert an Eigen::ComputationInfo enumerated type into a human-readable string
 * 
 *  @param ci enumerated type
 *  @return a human-readable string representation of the enumerated type
 */
std::string eigenCompInfo2Str(Eigen::ComputationInfo ci){
    switch(ci){
        case Eigen::Success: return std::string("Computation was successful");
        case Eigen::NumericalIssue: return std::string("the provided data did not satisfy the prerequisites");
        case Eigen::NoConvergence: return std::string("Iterative procedure did not converge");
        case Eigen::InvalidInput: return std::string("The inputs are invalid, or the algorithm has been improperly called");
        default: return std::string("Unknown ComputationInfo enumerator");
    }
}//====================================================

/** \} */ // END of util group
} // End of astrohelion namespace
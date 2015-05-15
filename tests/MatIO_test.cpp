/**
 *  Test out the matio tools
 */

#include <stdlib.h>
#include <stdio.h>
#include "matio.h"

int main(int argc, char **argv){
	mat_t *matfp;
	matvar_t *matvar;
	size_t dims[2] = {5, 2};
	double 	x[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
			y[10] = {11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
	struct mat_complex_split_t z = {x, y};

	// Create the .mat file
	matfp = Mat_CreateVer("Test.mat", NULL, MAT_FT_DEFAULT);
	if(NULL == matfp){
		fprintf(stderr, "Error creating MAT File\n");
		return EXIT_FAILURE;
	}

	// Save the x vector to the file
	matvar = Mat_VarCreate("x", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, x, 0);
	if(NULL == matvar){
		fprintf(stderr, "Error creating variable for 'x'\n");
	}else{
		Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
		Mat_VarFree(matvar);
	}

	// Save the y vector to the file
	matvar = Mat_VarCreate("y", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, y, 0);
	if(NULL == matvar){
		fprintf(stderr, "Error creating variable for 'y'\n");
	}else{
		Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
		Mat_VarFree(matvar);
	}

	// Save a variable z that has real parts x and imaginary parts y
	matvar = Mat_VarCreate("z", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &z, MAT_F_COMPLEX);
	if(NULL == matvar){
		fprintf(stderr, "Error creating variable for 'z'\n");
	}else{
		Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
		Mat_VarFree(matvar);
	}

	// Save a 3D array containing three 2x2 matrices. The first should have a column of ones in the first col.
	double q[12] = {1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1};
	size_t dims3D[3] = {2,2,3};
	matvar = Mat_VarCreate("I", MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims3D, q, 0);
	if(NULL == matvar){
		fprintf(stderr, "Error creating variable for 'q'\n");
	}else{
		Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
		Mat_VarFree(matvar);
	}

	Mat_Close(matfp);
	return EXIT_SUCCESS;
}
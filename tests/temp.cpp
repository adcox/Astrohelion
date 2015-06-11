#include <iostream>
#include <gsl/gsl_matrix.h>

using namespace std;

int main(void){
	
	double data[] = {1,2,3,4};
	gsl_matrix *m = gsl_matrix_alloc(2,2);
	copy(data, data+4, m->data);

	printf("Matrix:\n%f %f\n%f %f\n", gsl_matrix_get(m, 0, 0), gsl_matrix_get(m, 0, 1),
		gsl_matrix_get(m, 1, 0), gsl_matrix_get(m, 1, 1));
	return 0;
}
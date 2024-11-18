/*
 * DESC: Module for naive matrix multiplication.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */
#include <math.h>  // for pow()
#include <stdlib.h>
#include <time.h>
#include <cblas.h> // for CBLAS

#include "../include/naive_matmat.h"
#include "../include/IO.h"

/* generate a random floating point number from min to max */
double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

void gen_rand_matrix(double *A, const size_t m, const size_t n){
	for(size_t i=0; i<m*n; i++){
		A[i] = randfrom(0.,1.);    
	}
}

int main(int argc, char *argv[]) {
	size_t N = 4;	// max power dimension of matrix
	for (size_t i = 1; i < N; i++) {
		const size_t n = pow(2, i);
		const size_t m = n+1;
		const size_t k = n-1;
		double *A = malloc(m * n * sizeof(double));
		double *B = malloc(n * k * sizeof(double));
		double *C = malloc(m * k * sizeof(double));
		gen_rand_matrix(A,m,n);
		gen_rand_matrix(B,n,k);

		printf("# naive_matmat() # \n");
		naive_matmat(A,B,C,m,n,k);
		print_mat_double(C,m,k);
		
		printf("# cblas_dgemm() (matmat) # \n");
		double alpha = 1.;
		double beta = 1.;
		bli_dgemm(BLIS_TRANSPOSE, BLIS_TRANSPOSE, m, k, n, &alpha, A, 1, m, B, 1, n, &beta, C, 1, m);

		print_mat_double(C,m,k);


		// Test 1 and write time to file.
		// Test 2 and write time to same file, next column.
		// Test ...

		free(A);
		free(B);
		free(C);
	}
}

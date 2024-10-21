/*
 * DESC: Module for naive matrix multiplication.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */
#include <math.h>  // for pow()

int main(int argc, char *argv[]) {
	size_t N = 10;	// max power dimension of matrix
	for (size_t i = 1; i < N, i++) {
		size_t n = pow(2, i);
		double *A = malloc(n * n * sizeof(double));
		double *B = malloc(n * n * sizeof(double));
		double *C = malloc(n * n * sizeof(double));
		gen_rand_matrix(A);
		// Test 1 and write time to file.
		// Test 2 and write time to same file, next column.
		// Test ...

		free(A);
		free(B);
		free(C);
	}
}

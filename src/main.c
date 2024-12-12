/*
 * DESC: Result of file compilation to be executed. Runs tests on all
 * implementations on different matrix sizes.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */
#include <math.h>    // for pow()
#include <stddef.h>  // for size_t
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../include/IO.h"
#include "../include/test.h"

int main(int argc, char *argv[]) {
	const size_t N = 12;  // max power dimension of matrix
	srand(time(NULL));
	for (size_t i = 2; i < N; i++) {
		const size_t n = pow(2, i) - 3;
		printf("Size of test: %zu \n", i);
		flush_cache();
		printf("test_naive_matmat time: %lf \n",
		       test_naive_matmat(i, 0.001));
		flush_cache();
		printf("test_strassen_matmat time: %lf \n",
		       test_strassen_matmat(i, 0.001));

		// Test invert functions

		// Initialize a random invertible A
		double *A = malloc(n * n * sizeof(double));
		gen_rand_matrix(A, n, n);
		int invertible = 0;
		double *A_copy = malloc(n * n * sizeof(double));
		while (!invertible) {
			gen_rand_matrix(A, n, n);
			memcpy(A_copy, A, n * n * sizeof(double));
			invertible = is_invertible(A_copy, n);
		}
		free(A_copy);

		// Test the invert functions
		flush_cache();
		printf("test_strassen_invert_naive_matmat time: %lf\n",
		       test_strassen_invert_naive_matmat(&A, n, 0.001));
		flush_cache();
		printf("test_strassen_invert_strassen_matmat time: %lf \n",
		       test_strassen_invert_strassen_matmat(&A, n, 0.001));
		flush_cache();
		printf("test_lu_invert time: %lf \n",
		       test_lu_invert(A, n, 0.001));
		free(A);
	}
}

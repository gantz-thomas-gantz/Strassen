/*
 * DESC: Module for naive matrix multiplication.
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

// Done: changes in main, test_naive_matmat, test_strassen_matmat
// TODO: Test LU invert

int main(int argc, char *argv[]) {
    
	size_t N = 5; // default max power dimension of matrix
    // Check if an argument is provided for N
    if (argc > 1) {
        N = strtoul(argv[1], NULL, 10);
        if (N == 0) {
            fprintf(stderr, "Invalid input for N. Please provide a positive integer.\n");
            return EXIT_FAILURE;
        }
    }

	double tolerance = 1e-3;
	srand(time(NULL));	 // random seed

	// Files to save results for plotting
	FILE *file_matmat = fopen("data/matmat.txt", "w");
	FILE *file_matinv = fopen("data/matinv.txt", "w");
	
	for (size_t i = 1; i <= N; i++) {
		
		const size_t n = pow(2, i) - 1; // Test matrices how are not of size (2^n,2^n)

		printf("# ######################\n");
		printf("# Size of test: %zu \n", i);
		printf("# ######################\n");
		printf("\n");

		/* ####################################################### */

		printf("### TEST 1 : Matrix Multiplication\n");

		// Initialize data
		const size_t m = n + 1;
		const size_t k = n - 1;
		double *A_mul = malloc(m * n * sizeof(double));
		double *B_mul = malloc(n * k * sizeof(double));
		gen_rand_matrix(A_mul, m, n);
		gen_rand_matrix(B_mul, n, k);

		// Perform tests
		flush_cache();
        double naive_time = test_naive_matmat(&A_mul, &B_mul, m, n, k, tolerance);
		flush_cache();
        double strassen_time = test_strassen_matmat(&A_mul, &B_mul, m, n, k, tolerance);

		// Write to console
		printf("- naive_matmat :    %.5lf\n", naive_time);
		printf("- strassen_matmat : %.5lf\n", strassen_time);
		printf("\n");

		// Write to file
		fprintf(file_matmat, "%zu %.5lf %.5lf\n", i, naive_time, strassen_time);

		// Free data
		free(A_mul);
		free(B_mul);
		
		/* ####################################################### */
		
		printf("### TEST 2 : Matrix Inversion\n");

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

		// Perform tests
		flush_cache();
		double time_strassen_invert_naive_matmat = test_strassen_invert_naive_matmat(&A, n, tolerance);
		flush_cache();
        double time_strassen_invert_strassen_matmat = test_strassen_invert_strassen_matmat(&A, n, tolerance);
        

		// Write to console
		printf("- strassen_invert_naive_matmat :    %.5lf\n", time_strassen_invert_naive_matmat);
		printf("- strassen_invert_strassen_matmat : %.5lf\n", time_strassen_invert_strassen_matmat);
		printf("\n");

		// Write to file
		fprintf(file_matinv, "%zu %.5lf %.5lf\n", i, time_strassen_invert_naive_matmat, time_strassen_invert_strassen_matmat);
	
		// free data
		free(A);

	}

	fclose(file_matinv);
	fclose(file_matmat);
}
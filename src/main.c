/*
 * DESC: Main file to execute all tests.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../include/IO.h"
#include "../include/test.h"

int main(int argc, char *argv[]) {
	size_t N = 5;  // default max power dimension of matrix

	// Check if an argument is provided for N
	if (argc > 1) {
		N = strtoul(argv[1], NULL, 10);
		if (N == 0) {  // Handle invalid input
			fprintf(stderr,
				"Invalid input for N. Please provide a "
				"positive integer.\n");
			return EXIT_FAILURE;
		}
	}

	double tolerance = 1e-3;  // Set test tolerance level
	srand(time(
	    NULL));  // Seed the random number generator with the current time

	// Open files to save test results for plotting
	FILE *file_matmat = fopen("matmat.txt", "w");
	FILE *file_matinv = fopen("matinv.txt", "w");

	for (size_t i = 2; i <= N; i++) {
		const size_t n =
		    pow(2, i) - 1;  // Determine the size of test matrices (not
				    // standard 2^n,2^n)

		printf("# ######################\n");
		printf("# Size of test: %zu \n", i);
		printf("# ######################\n");
		printf("\n");

		/* ####################################################### */
		printf("### TEST 1 : Matrix Multiplication\n");

		// Initialize data for matrix multiplication
		const size_t m = n + 1;
		const size_t k = n - 1;
		double *A_mul = malloc(m * n * sizeof(double));
		double *B_mul = malloc(n * k * sizeof(double));

		// Generate random matrices
		gen_rand_matrix(A_mul, m, n);
		gen_rand_matrix(B_mul, n, k);

		// Flush cache to ensure fair timing
		flush_cache();

		// Perform naive matrix multiplication test
		double naive_time =
		    test_naive_matmat(&A_mul, &B_mul, m, n, k, tolerance);

		// Flush cache to ensure fair timing
		flush_cache();

		// Perform Strassen matrix multiplication test
		double strassen_time =
		    test_strassen_matmat(&A_mul, &B_mul, m, n, k, tolerance);

		// Output the test results to console
		printf("- naive_matmat :    %.5lf\n", naive_time);
		printf("- strassen_matmat : %.5lf\n", strassen_time);
		printf("\n");

		// Write test results to the corresponding file
		fprintf(file_matmat, "%zu %lf %lf\n", i, naive_time,
			strassen_time);

		// Free allocated memory for matrix multiplication
		free(A_mul);
		free(B_mul);

		/* ####################################################### */
		printf("### TEST 2 : Matrix Inversion\n");

		// Initialize a random invertible matrix A
		double *A = malloc(n * n * sizeof(double));
		gen_rand_matrix(A, n, n);
		int invertible = 0;
		double *A_copy = malloc(n * n * sizeof(double));

		// Ensure matrix A is invertible
		while (!invertible) {
			gen_rand_matrix(A, n, n);
			memcpy(A_copy, A, n * n * sizeof(double));
			invertible = is_invertible(A_copy, n);
		}
		free(A_copy);

		// Flush cache to ensure fair timing
		flush_cache();

		// Perform Strassen inversion with naive matrix multiplication
		double time_strassen_invert_naive_matmat =
		    test_strassen_invert_naive_matmat(&A, n, tolerance);

		flush_cache();

		// Perform Strassen inversion with Strassen's matrix
		// multiplication
		double time_strassen_invert_strassen_matmat =
		    test_strassen_invert_strassen_matmat(&A, n, tolerance);

		flush_cache();

		// Perform LU-based inversion
		double time_lu_invert = test_lu_invert(A, n, tolerance);

		// Output results to console
		printf("- strassen_invert_naive_matmat :    %.5lf\n",
		       time_strassen_invert_naive_matmat);
		printf("- strassen_invert_strassen_matmat : %.5lf\n",
		       time_strassen_invert_strassen_matmat);
		printf("- lu_invert :                       %.5lf\n",
		       time_lu_invert);
		printf("\n");

		// Write test results to file
		fprintf(file_matinv, "%zu %lf %lf %lf\n", i, time_lu_invert,
			time_strassen_invert_naive_matmat,
			time_strassen_invert_strassen_matmat);

		free(A);
	}

	// Close the opened files
	fclose(file_matinv);
	fclose(file_matmat);
}


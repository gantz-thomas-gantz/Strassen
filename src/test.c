/*
 * DESC: Module for testing all project implementations.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "../include/IO.h"
#include "../include/naive_lu.h"
#include "../include/naive_matmat.h"
#include "../include/strassen_inv.h"
#include "../include/strassen_matmat.h"

void flush_cache() {
	const size_t cache_size = 32 * 1024 * 1024;  // 32 MB (adjust if needed)
	char *dummy = (char *)malloc(cache_size);
	if (dummy) {
		for (size_t i = 0; i < cache_size;
		     i += 64) {	 // Step through memory in cache line-sized
				 // increments
			dummy[i] = i;
		}
		free(dummy);
	}
}

static double randfrom() {
	double a =
	    (double)rand() / RAND_MAX;	// Generate a random double in [0, 1]
	int p_o_m =
	    a < 0.5 ? 1
		    : 0;  // Decide if the number will be positive or negative
	return p_o_m * (-1) +
	       (double)rand() /
		   RAND_MAX;  // Return a random double with the chosen sign
}

void gen_rand_matrix(double *A, const size_t m, const size_t n) {
	for (size_t i = 0; i < m * n;
	     i++) {  // Fill the matrix with random values
		A[i] = randfrom();
	}
}

int compare_mat(const double *const A, const double *const B, const size_t m,
		const size_t n, const double eps) {
	for (size_t i = 0; i < m * n; i++)  // Iterate over all elements
		if (fabs(A[i] - B[i]) > eps)
			return 0;  // Return 0 if elements differ by more than
				   // eps
	return 1;		   // Matrices are equal within tolerance
}

double test_naive_matmat(double **A, double **B, const size_t m, const size_t n,
			 const size_t k, const double eps) {
	double *C_gt = malloc(
	    m * k * sizeof(double));  // Ground truth matrix (using cblas)
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, k, n, 1., *A,
		    n, *B, k, 0., C_gt, k);

	double *C = malloc(m * k * sizeof(double));  // Result matrix
	clock_t start = clock();		     // Record start time
	naive_matmat(*A, *B, C, m, n,
		     k);	// Perform naive matrix multiplication
	clock_t end = clock();	// Record end time
	double time_spent =
	    (double)(end - start) / CLOCKS_PER_SEC;  // Calculate elapsed time

	double result = -1.0;
	if (compare_mat(C, C_gt, m, k, eps))
		result = time_spent;  // Validate result

	free(C);
	free(C_gt);

	return result;
}

double test_strassen_matmat(double **A, double **B, const size_t m,
			    const size_t n, const size_t k, const double eps) {
	double *C_gt = malloc(m * k * sizeof(double));	// Ground truth matrix
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, k, n, 1., *A,
		    n, *B, k, 0., C_gt, k);

	double *C = malloc(m * k * sizeof(double));  // Result matrix
	clock_t start = clock();		     // Record start time
	strassen_matmat(A, B, &C, m, n,
			k);	// Perform Strassen's matrix multiplication
	clock_t end = clock();	// Record end time
	double time_spent =
	    (double)(end - start) / CLOCKS_PER_SEC;  // Calculate elapsed time

	double result = -1.0;
	if (compare_mat(C, C_gt, m, k, eps))
		result = time_spent;  // Validate result

	free(C);
	free(C_gt);

	return result;
}

int is_invertible(double *A, int n) {
	int *ipiv = (int *)malloc(n * sizeof(int));  // Pivot indices
	int info;

	// Perform LU decomposition
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A, n, ipiv);

	free(ipiv);

	// If info > 0, the matrix is singular
	return info == 0;  // Returns 1 if invertible, 0 if not
}

double test_strassen_invert_strassen_matmat(double **A, const size_t n,
					    const double eps) {
	double *inverse_A = calloc(
	    n * n,
	    sizeof(double));  // Allocate memory for Strassen's inverted matrix
	double *inverse_A_gt = calloc(
	    n * n, sizeof(double));  // Allocate memory for ground truth inverse
	int *ipiv = malloc(
	    n * sizeof(int));  // Pivot indices for ground truth inversion

	// Compute ground truth inverse using LAPACK
	memcpy(inverse_A_gt, *A,
	       n * n * sizeof(double));	 // Copy input matrix to ground truth
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, inverse_A_gt, n,
		       ipiv);  // Perform LU decomposition
	LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, inverse_A_gt, n,
		       ipiv);  // Compute inverse from LU factors

	clock_t start = clock();  // Record start time
	strassen_invert_strassen_matmat(A, &inverse_A,
					n);  // Perform Strassen's inversion
	clock_t end = clock();		     // Record end time
	double time_spent =
	    (double)(end - start) / CLOCKS_PER_SEC;  // Calculate elapsed time

	double result = -1.0;
	if (compare_mat(inverse_A, inverse_A_gt, n, n,
			eps))  // Validate result against ground truth
		result = time_spent;

	free(ipiv);
	free(inverse_A);
	free(inverse_A_gt);

	return result;
}

double test_strassen_invert_naive_matmat(double **A, const size_t n,
					 const double eps) {
	double *inverse_A = calloc(
	    n * n, sizeof(double));  // Allocate memory for Naive inversion
	double *inverse_A_gt = calloc(
	    n * n, sizeof(double));  // Allocate memory for ground truth inverse
	int *ipiv = malloc(
	    n * sizeof(int));  // Pivot indices for ground truth inversion

	// Compute ground truth inverse using LAPACK
	memcpy(inverse_A_gt, *A,
	       n * n * sizeof(double));	 // Copy input matrix to ground truth
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, inverse_A_gt, n,
		       ipiv);  // Perform LU decomposition
	LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, inverse_A_gt, n,
		       ipiv);  // Compute inverse from LU factors

	clock_t start = clock();  // Record start time
	strassen_invert_naive_matmat(A, &inverse_A,
				     n);  // Perform Naive inversion
	clock_t end = clock();		  // Record end time
	double time_spent =
	    (double)(end - start) / CLOCKS_PER_SEC;  // Calculate elapsed time

	double result = -1.0;
	if (compare_mat(inverse_A, inverse_A_gt, n, n,
			eps))  // Validate result against ground truth
		result = time_spent;

	free(ipiv);
	free(inverse_A);
	free(inverse_A_gt);

	return result;
}

double test_lu_invert(const double *const A, const size_t n, const double eps) {
	double *inverse_A =
	    calloc(n * n, sizeof(double));  // Allocate memory for LU inversion
	double *inverse_A_gt = calloc(
	    n * n, sizeof(double));  // Allocate memory for ground truth inverse
	int *ipiv = malloc(
	    n * sizeof(int));  // Pivot indices for ground truth inversion

	// Compute ground truth inverse using LAPACK
	memcpy(inverse_A_gt, A,
	       n * n * sizeof(double));	 // Copy input matrix to ground truth
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, inverse_A_gt, n,
		       ipiv);  // Perform LU decomposition
	LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, inverse_A_gt, n,
		       ipiv);  // Compute inverse from LU factors

	clock_t start = clock();     // Record start time
	lu_invert(A, inverse_A, n);  // Perform LU-based inversion
	clock_t end = clock();	     // Record end time
	double time_spent =
	    (double)(end - start) / CLOCKS_PER_SEC;  // Calculate elapsed time

	double result = -1.0;
	if (compare_mat(inverse_A, inverse_A_gt, n, n,
			eps))  // Validate result against ground truth
		result = time_spent;

	free(ipiv);
	free(inverse_A);
	free(inverse_A_gt);

	return result;
}


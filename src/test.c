/*
 * DESC: Module for naive matrix multiplication.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */
// TODO: generate same random matrix for same functionality such that tests are
// comparable (done for inversion, todo for multiplication)
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <stddef.h>  // for size_t
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "../include/IO.h"
#include "../include/lu_naive.h"
#include "../include/naive_matmat.h"
#include "../include/strassen_inv.h"
#include "../include/strassen_matmat.h"

void flush_cache() {
	const size_t cache_size = 32 * 1024 * 1024;  // 32 MB (adjust if needed)
	char *dummy = (char *)malloc(cache_size);
	if (dummy) {
		for (size_t i = 0; i < cache_size; i += 64) {
			dummy[i] = i;
		}
		free(dummy);
	}
}

/* generate a random floating point number from min to max */
double randfrom(double min, double max) {
	return (double)rand() / RAND_MAX;  // * (max - min);
}

void gen_rand_matrix(double *A, const size_t m, const size_t n) {
	for (size_t i = 0; i < m * n; i++) {
		A[i] = randfrom(0., 1.);
	}
}

int compare_mat(const double *const A, const double *const B, const size_t m,
		const size_t n, const double eps) {
	for (size_t i = 0; i < m * n; i++)
		if (fabs(A[i] - B[i]) > eps) return 0;
	return 1;
}

double test_naive_matmat(const size_t N, const double eps) {
	const size_t n = pow(2, N);
	const size_t m = n + 1;
	const size_t k = n - 1;
	double *A = malloc(m * n * sizeof(double));
	double *B = malloc(n * k * sizeof(double));
	double *C = malloc(m * k * sizeof(double));
	double *C_gt = malloc(m * k * sizeof(double));
	gen_rand_matrix(A, m, n);
	gen_rand_matrix(B, n, k);
	clock_t start = clock();  // Record start time
	naive_matmat(A, B, C, m, n, k);
	clock_t end = clock();	// Record end time
	double time_spent =
	    (double)(end - start) / CLOCKS_PER_SEC;  // Calculate elapsed time
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, k, n, 1., A,
		    n, B, k, 0., C_gt, k);
	double result = -1.0;
	if (compare_mat(C, C_gt, m, k, eps)) result = time_spent;
	free(A);
	free(B);
	free(C);
	free(C_gt);

	return result;
}
double test_strassen_matmat(const size_t N, const double eps) {
	const size_t n = pow(2, N);
	const size_t m = n + 1;
	const size_t k = n - 1;
	double *A = malloc(m * n * sizeof(double));
	double *B = malloc(n * k * sizeof(double));
	double *C = malloc(m * k * sizeof(double));
	double *C_gt = malloc(m * k * sizeof(double));
	gen_rand_matrix(A, m, n);
	gen_rand_matrix(B, n, k);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, k, n, 1., A,
		    n, B, k, 0., C_gt, k);

	clock_t start = clock();  // Record start time
	strassen_matmat(&A, &B, &C, m, n, k);
	clock_t end = clock();	// Record end time
	double time_spent =
	    (double)(end - start) / CLOCKS_PER_SEC;  // Calculate elapsed time
	double result = -1.0;
	if (compare_mat(C, C_gt, m, k, eps)) result = time_spent;
	free(A);
	free(B);
	free(C);
	free(C_gt);

	return result;
};

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
	double *inverse_A = calloc(n * n, sizeof(double));
	double *inverse_A_gt = calloc(n * n, sizeof(double));
	int *ipiv = malloc(n * sizeof(int));  // Pivot indices
	memcpy(inverse_A_gt, *A, n * n * sizeof(double));
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, inverse_A_gt, n, ipiv);
	LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, inverse_A_gt, n, ipiv);

	clock_t start = clock();  // Record start time
	strassen_invert_strassen_matmat(A, &inverse_A, n);
	clock_t end = clock();	// Record end time
	double time_spent =
	    (double)(end - start) / CLOCKS_PER_SEC;  // Calculate elapsed time
	double result = -1.0;
	if (compare_mat(inverse_A, inverse_A_gt, n, n, eps))
		result = time_spent;
	free(ipiv);
	free(inverse_A);
	free(inverse_A_gt);
	return result;
};
double test_strassen_invert_naive_matmat(double **A, const size_t n,
					 const double eps) {
	double *inverse_A = calloc(n * n, sizeof(double));
	double *inverse_A_gt = calloc(n * n, sizeof(double));
	int *ipiv = malloc(n * sizeof(int));  // Pivot indices
	memcpy(inverse_A_gt, *A, n * n * sizeof(double));
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, inverse_A_gt, n, ipiv);
	LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, inverse_A_gt, n, ipiv);

	clock_t start = clock();  // Record start time
	strassen_invert_naive_matmat(A, &inverse_A, n);
	clock_t end = clock();	// Record end time
	double time_spent =
	    (double)(end - start) / CLOCKS_PER_SEC;  // Calculate elapsed time
	double result = -1.0;
	if (compare_mat(inverse_A, inverse_A_gt, n, n, eps))
		result = time_spent;
	free(ipiv);
	free(inverse_A);
	free(inverse_A_gt);
	return result;
};

double test_lu_invert(const double *const A, const size_t n, const double eps) {
	double *inverse_A = calloc(n * n, sizeof(double));
	double *inverse_A_gt = calloc(n * n, sizeof(double));
	int *ipiv = malloc(n * sizeof(int));  // Pivot indices
	memcpy(inverse_A_gt, *A, n * n * sizeof(double));
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, inverse_A_gt, n, ipiv);
	LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, inverse_A_gt, n, ipiv);

	clock_t start = clock();  // Record start time
	lu_invert(A, inverse_A, n);
	clock_t end = clock();	// Record end time
	double time_spent =
	    (double)(end - start) / CLOCKS_PER_SEC;  // Calculate elapsed time
	double result = -1.0;
	if (compare_mat(inverse_A, inverse_A_gt, n, n, eps))
		result = time_spent;
	free(ipiv);
	free(inverse_A);
	free(inverse_A_gt);
	return result;
};


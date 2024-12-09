/*
 * DESC: Module for naive matrix multiplication.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */
// TODO: generate same random matrix for same functionality such that tests are
// comparable
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <stddef.h>  // for size_t
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "../include/IO.h"
#include "../include/naive_matmat.h"
#include "../include/strassen_inv.h"
#include "../include/strassen_matmat.h"

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

// TODO: make sure matrix is invertible, else gen new matrix
double test_strassen_invert_strassen_matmat(const size_t N, const double eps) {
	const size_t n = pow(2, N);
	double *A = malloc(n * n * sizeof(double));
	double *inverse_A = malloc(n * n * sizeof(double));
	double *inverse_A_gt = malloc(n * n * sizeof(double));
	gen_rand_matrix(A, n, n);
	int *ipiv = malloc(n * sizeof(int));  // Pivot indices
	memcpy(inverse_A_gt, A, n * n * sizeof(double));
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, inverse_A_gt, n, ipiv);
	LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, inverse_A_gt, n, ipiv);

	clock_t start = clock();  // Record start time
	strassen_invert_strassen_matmat(A, inverse_A, n);
	clock_t end = clock();	// Record end time
	double time_spent =
	    (double)(end - start) / CLOCKS_PER_SEC;  // Calculate elapsed time
	double result = -1.0;
	if (compare_mat(inverse_A, inverse_A_gt, n, n, eps))
		result = time_spent;

	free(A);
	free(ipiv);
	free(inverse_A);
	free(inverse_A_gt);

	return result;
};
double test_strassen_invert_naive_matmat(const size_t N, const double eps) {
	const size_t n = pow(2, N);
	double *A = malloc(n * n * sizeof(double));
	double *inverse_A = malloc(n * n * sizeof(double));
	double *inverse_A_gt = malloc(n * n * sizeof(double));
	gen_rand_matrix(A, n, n);
	int *ipiv = malloc(n * sizeof(int));  // Pivot indices
	memcpy(inverse_A_gt, A, n * n * sizeof(double));
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, inverse_A_gt, n, ipiv);
	LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, inverse_A_gt, n, ipiv);

	clock_t start = clock();  // Record start time
	strassen_invert_naive_matmat(A, inverse_A, n);
	clock_t end = clock();	// Record end time
	double time_spent =
	    (double)(end - start) / CLOCKS_PER_SEC;  // Calculate elapsed time
	double result = -1.0;
	if (compare_mat(inverse_A, inverse_A_gt, n, n, eps))
		result = time_spent;

	free(A);
	free(ipiv);
	free(inverse_A);
	free(inverse_A_gt);

	return result;
};


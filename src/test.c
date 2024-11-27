/*
 * DESC: Module for naive matrix multiplication.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */

#include <cblas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../include/naive_matmat.h"

/* generate a random floating point number from min to max */
double randfrom(double min, double max) {
	double range = (max - min);
	double div = RAND_MAX / range;
	return min + (rand() / div);
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
	if (compare_mat(C, C_gt, m, k, eps))
		return time_spent;
	else
		return -1.0;
	free(A);
	free(B);
	free(C);
	free(C_gt);
}
double test_strassen_matmat(i, 0.001)){
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
	if (compare_mat(C, C_gt, m, k, eps))
		return time_spent;
	else
		return -1.0;
	free(A);
	free(B);
	free(C);
	free(C_gt);
};


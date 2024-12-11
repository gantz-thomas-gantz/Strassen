#include "block_utilities.h"

#include <assert.h>
#include <stddef.h>  // for size_t
#include <stdlib.h>

double *darray_add(const double *const A, const double *const B,
		   const size_t size, const double alpha) {
	double *C = (double *)malloc(size * sizeof(double));
	for (size_t i = 0; i < size; i++) C[i] = A[i] + alpha * B[i];
	return C;
}

void darray_block_add(const double *const A, double *C, const size_t start1,
		      const size_t start2, const size_t m, const size_t n,
		      const double alpha) {
	assert((start1 + (m / 2 - 1) * n + (n / 2 - 1) < m * n) ||
	       (start2 + (m / 2 - 1) * n + (n / 2 - 1) < m * n));
	for (size_t i = 0; i < m / 2; i++)
		for (size_t j = 0; j < n / 2; j++)
			C[i * n / 2 + j] = A[start1 + i * n + j] +
					   alpha * A[start2 + i * n + j];
}

double *create_block(const double *const A, const size_t start, const size_t m,
		     const size_t n) {
	assert((start + (m / 2 - 1) * n + (n / 2 - 1) < m * n));
	double *a = (double *)malloc(m / 2 * n / 2 * sizeof(double));
	for (size_t i = 0; i < m / 2; i++)
		for (size_t j = 0; j < n / 2; j++)
			a[i * n / 2 + j] = A[start + i * n + j];
	return a;
}

void mat_inplace_block_add(double *C, const double *const a,
			   const double *const b, const size_t start,
			   const size_t m, const size_t n, const double alpha,
			   const double beta) {
	assert((start + (m / 2 - 1) * n + (n / 2 - 1) < m * n));
	if (b == NULL) {
		for (size_t i = 0; i < m / 2; i++) {
			for (size_t j = 0; j < n / 2; j++) {
				C[start + i * n + j] +=
				    alpha * a[i * n / 2 + j];
			}
		}
	} else {
		for (size_t i = 0; i < m / 2; i++) {
			for (size_t j = 0; j < n / 2; j++) {
				C[start + i * n + j] +=
				    alpha * a[i * n / 2 + j] +
				    beta * b[i * n / 2 + j];
			}
		}
	}
}


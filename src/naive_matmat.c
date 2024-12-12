/*
 * DESC: Module for naive matrix multiplication.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */

#include <stddef.h>  // for size_t

void naive_matmat(double *A, double *B, double *C, const size_t m,
		  const size_t n, const size_t k) {
	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < k; j++) {
			C[i * k + j] = 0;
			for (size_t l = 0; l < n; l++) {
				C[i * k + j] += A[i * n + l] * B[l * k + j];
			}
		}
	}
}


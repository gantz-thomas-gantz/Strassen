/*
 * DESC: Module for naive matrix multiplication.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */

#include <stddef.h>

void naive_matmat(double *A, double *B, double *C, const size_t m,
		  const size_t n, const size_t k) {
	// Iterate over rows of A (m rows)
	for (size_t i = 0; i < m; i++) {
		// Iterate over columns of B (k columns)
		for (size_t j = 0; j < k; j++) {
			// Initialize the result matrix C at row i, column j
			C[i * k + j] = 0;
			// Perform dot product of row i from A and column j from
			// B
			for (size_t l = 0; l < n; l++) {
				C[i * k + j] += A[i * n + l] * B[l * k + j];
			}
		}
	}
}


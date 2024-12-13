/*
 * DESC: Module for strassen matrix multiplication.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */
#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/block_utilities.h"
#include "../include/naive_matmat.h"

static void pad_matrix(double **A, size_t m, size_t n, size_t *new_m,
		       size_t *new_n) {
	*new_m = (m % 2 == 0) ? m : m + 1;  // Ensure rows are even
	*new_n = (n % 2 == 0) ? n : n + 1;  // Ensure columns are even

	if (!(m % 2 == 0) || !(n % 2 == 0)) {
		// Allocate new memory for the padded matrix
		double *padded_A =
		    (double *)calloc((*new_m) * (*new_n), sizeof(double));

		// Copy the original data into the padded matrix
		for (size_t i = 0; i < m; i++)
			for (size_t j = 0; j < n; j++)
				padded_A[i * (*new_n) + j] = (*A)[i * n + j];

		free(*A);	// Free the old memory
		*A = padded_A;	// Point to the new padded matrix
	}
}

static void depad_matrix(double **padded_A, size_t m, size_t n, size_t og_m,
			 size_t og_n) {
	// If no padding was applied, no action is needed
	if (m == og_m && n == og_n) {
		return;
	}

	// Allocate memory for the original matrix without padded extra space
	double *A = (double *)malloc(og_m * og_n * sizeof(double));

	// Copy valid elements back from the padded matrix
	for (size_t i = 0; i < og_m; i++) {
		for (size_t j = 0; j < og_n; j++) {
			A[i * og_n + j] = (*padded_A)[i * n + j];
		}
	}

	free(*padded_A);  // Free padded memory
	*padded_A = A;	  // Point to the depadded version
}

void strassen_matmat(double **A, double **B, double **C, size_t m, size_t n,
		     size_t k) {
	// If matrices are too small, fallback to naive matrix multiplication
	if (m < 512 && n < 512 && k < 512) {
		naive_matmat(*A, *B, *C, m, n, k);
	} else if (m == 1 && n == 1) {
		// Handle base case of single-element multiplication
		for (size_t i = 0; i < m * k; i++) (*C)[i] = (*A)[0] * (*B)[i];
	} else if (n == 1 && k == 1) {
		// Handle special single-row/column edge case
		for (size_t i = 0; i < m * k; i++) (*C)[i] = (*A)[i] * (*B)[0];
	} else {
		// Initialize result matrix to zero
		for (size_t i = 0; i < m * k; i++) (*C)[i] = 0;

		// Pad matrices to ensure they have compatible even dimensions
		size_t og_m = m;
		size_t og_n = n;
		size_t og_k = k;
		pad_matrix(A, og_m, og_n, &m, &n);
		pad_matrix(B, og_n, og_k, &n, &k);
		pad_matrix(C, og_m, og_k, &m, &k);

		// Subdivide matrices into blocks for Strassen's recursive
		// computation
		const size_t start_a = 0;
		const size_t start_x = 0;
		const size_t start_r11 = 0;
		const size_t start_b = n / 2;
		const size_t start_y = k / 2;
		const size_t start_r12 = k / 2;
		const size_t start_c = m / 2 * n;
		const size_t start_z = n / 2 * k;
		const size_t start_r21 = m / 2 * k;
		const size_t start_d = start_c + start_b;
		const size_t start_t = start_y + start_z;
		const size_t start_r22 = start_r12 + start_r21;

		// Extract matrix blocks (allocate contiguous memory dynamically
		// for sub-matrices)
		double *a = create_block(*A, start_a, m, n);
		double *d = create_block(*A, start_d, m, n);
		double *y = create_block(*B, start_y, n, k);
		double *z = create_block(*B, start_z, n, k);

		double *tempA =
		    (double *)malloc(m / 2 * n / 2 * sizeof(double));
		double *tempB =
		    (double *)malloc(n / 2 * k / 2 * sizeof(double));

		// Strassen's recursive multiplications
		double *q1 = (double *)malloc(m * k / 4 * sizeof(double));
		darray_block_add(*B, tempB, start_x, start_z, n, k, 1.0);
		strassen_matmat(&a, &tempB, &q1, m / 2, n / 2, k / 2);
		free(a);

		double *q2 = (double *)malloc(m * k / 4 * sizeof(double));
		darray_block_add(*B, tempB, start_y, start_t, n, k, 1.0);
		strassen_matmat(&d, &tempB, &q2, m / 2, n / 2, k / 2);
		free(d);

		double *q3 = (double *)malloc(m * k / 4 * sizeof(double));
		darray_block_add(*A, tempA, start_d, start_a, m, n, -1.0);
		darray_block_add(*B, tempB, start_z, start_y, n, k, -1.0);
		strassen_matmat(&tempA, &tempB, &q3, m / 2, n / 2, k / 2);

		double *q4 = (double *)malloc(m * k / 4 * sizeof(double));
		darray_block_add(*A, tempA, start_b, start_d, m, n, -1.0);
		darray_block_add(*B, tempB, start_z, start_t, n, k, 1.0);
		strassen_matmat(&tempA, &tempB, &q4, m / 2, n / 2, k / 2);

		double *q5 = (double *)malloc(m * k / 4 * sizeof(double));
		darray_block_add(*A, tempA, start_b, start_a, m, n, -1.0);
		strassen_matmat(&tempA, &z, &q5, m / 2, n / 2, k / 2);
		free(z);

		double *q6 = (double *)malloc(m * k / 4 * sizeof(double));
		darray_block_add(*A, tempA, start_c, start_a, m, n, -1.0);
		darray_block_add(*B, tempB, start_x, start_y, n, k, 1.0);
		strassen_matmat(&tempA, &tempB, &q6, m / 2, n / 2, k / 2);

		double *q7 = (double *)malloc(m * k / 4 * sizeof(double));
		darray_block_add(*A, tempA, start_c, start_d, m, n, -1.0);
		strassen_matmat(&tempA, &y, &q7, m / 2, n / 2, k / 2);
		free(y);

		// Calculate the R blocks
		// R11
		mat_inplace_block_add(*C, q1, q5, start_r11, m, k, 1.0, 1.0);
		// R21
		mat_inplace_block_add(*C, q1, q3, start_r21, m, k, 1.0, 1.0);
		mat_inplace_block_add(*C, q6, q7, start_r21, m, k, 1.0, -1.0);
		// R12
		mat_inplace_block_add(*C, q2, q3, start_r12, m, k, 1.0, 1.0);
		mat_inplace_block_add(*C, q4, q5, start_r12, m, k, 1.0, -1.0);
		// R22
		mat_inplace_block_add(*C, q2, q7, start_r22, m, k, 1.0, 1.0);

		free(q1);
		free(q2);
		free(q3);
		free(q4);
		free(q5);
		free(q6);
		free(q7);

		// Cleanup padding
		depad_matrix(A, m, n, og_m, og_n);
		depad_matrix(B, n, k, og_n, og_k);
		depad_matrix(C, m, k, og_m, og_k);
	}
}


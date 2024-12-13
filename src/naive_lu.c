/*
 * DESC: Module for LU decomposition and inversion.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */
#include "naive_lu.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/IO.h"

// Initialize T with a copy of matrix A
static void init_T(const double *A, double *T, const size_t n) {
	memcpy(T, A, sizeof(double) * n * n);  // Copy the contents of A to T
}

// Perform the elimination step for the i-th row
static void eliminate_step_i(double *T, const size_t n, int i) {
	for (size_t j = i + 1; j < n; j++) {
		T[j * n + i] =
		    T[j * n + i] / T[i * n + i];  // Compute multiplier
		for (size_t k = i + 1; k < n; k++) {
			T[j * n + k] =
			    T[j * n + k] -
			    T[j * n + i] *
				T[i * n + k];  // Update remaining matrix
		}
	}
}

// Perform LU decomposition of matrix A, storing the result in T
void lu_decomposition(const double *const A, double *T, const size_t n) {
	init_T(A, T, n);  // Initialize T with a copy of A
	for (size_t i = 0; i < n - 1; i++) {
		eliminate_step_i(T, n, i);  // Apply elimination for each row
	}
}

// Compute the inverse of matrix A using its LU decomposition
void lu_invert(const double *const A, double *inverse_A, const size_t n) {
	double *T = (double *)malloc(sizeof(double) * n *
				     n);  // Allocate space for LU matrix
	lu_decomposition(A, T, n);	  // Perform LU decomposition

	double *temp = (double *)malloc(
	    sizeof(double) * n);  // Temporary storage for solving equations

	// Solve for each column of the inverse matrix
	for (size_t col = 0; col < n; col++) {
		// Initialize the current column as a canonical basis vector
		for (size_t i = 0; i < n; i++) {
			temp[i] =
			    (i == col) ? 1.0 : 0.0;  // Identity matrix column
		}

		// Forward substitution to solve L * Y = e_col
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < i; j++) {
				temp[i] -= T[i * n + j] *
					   temp[j];  // Subtract known terms
			}
		}

		// Backward substitution to solve U * X = Y
		for (int i = n - 1; i >= 0; i--) {
			for (size_t j = i + 1; j < n; j++) {
				temp[i] -= T[i * n + j] *
					   temp[j];  // Subtract known terms
			}
			temp[i] /= T[i * n + i];  // Divide by diagonal element
		}

		// Write the solution into the appropriate column of the inverse
		// matrix
		for (size_t i = 0; i < n; i++) {
			inverse_A[i * n + col] = temp[i];
		}
	}

	free(T);
	free(temp);
}


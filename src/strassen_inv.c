/*
 * DESC: Module for naive matrix multiplication.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */

#include <stdio.h>
#include <stdlib.h>

#include "../include/IO.h"
#include "../include/block_utilities.h"
#include "../include/naive_matmat.h"
#include "../include/strassen_matmat.h"
// TODO: Make identity padding

static void id_pad_matrix(double **A, size_t n, size_t *new_n) {
	*new_n = (n % 2 == 0) ? n : n + 1;
	if (!(n % 2 == 0)) {
		// Calculate new dimensions
		// Allocate memory for the new padded matrix
		double *padded_A =
		    (double *)calloc((*new_n) * (*new_n), sizeof(double));
		// Copy elements from the original matrix to the padded matrix
		for (size_t i = 0; i < n; i++)
			for (size_t j = 0; j < n; j++)
				padded_A[i * (*new_n) + j] = (*A)[i * n + j];
		for (size_t i = n; i < *new_n; i++)
			padded_A[i * (*new_n) + i] = 1;
		free(*A);
		*A = padded_A;
	}
}

void id_depad_matrix(double **padded_A, const size_t n, const size_t og_n) {
	// If the padded dimensions match the original, no need to depad
	if (n == og_n) {
		return;
	}
	// Allocate memory for the original matrix
	double *A = (double *)malloc(og_n * og_n * sizeof(double));
	// Copy elements from the padded matrix to the original matrix
	for (size_t i = 0; i < og_n; i++) {
		for (size_t j = 0; j < og_n; j++) {
			A[i * og_n + j] = (*padded_A)[i * n + j];
		}
	}
	free(*padded_A);
	*padded_A = A;
}
void strassen_invert_strassen_matmat(double **A, double **inverse_A, size_t n) {
	if (n == 1) {
		*inverse_A[0] = 1 / *A[0];
		return;
	}
	size_t og_n = n;
	printf("Before \n");
	print_mat_double(*A, n, n);
	id_pad_matrix(A, og_n, &n);
	printf("After \n");
	print_mat_double(*A, n, n);
	printf("\n");
	id_pad_matrix(inverse_A, og_n, &n);
	const size_t start_a = 0;
	const size_t start_b = n / 2;
	const size_t start_c = n / 2 * n;
	const size_t start_d = start_c + start_b;
	double *a = create_block(*A, start_a, n, n);
	double *b = create_block(*A, start_b, n, n);
	double *c = create_block(*A, start_c, n, n);
	double *d = create_block(*A, start_d, n, n);
	double *e = (double *)calloc(n * n / 4, sizeof(double));
	double *t = (double *)calloc(n * n / 4, sizeof(double));
	strassen_invert_strassen_matmat(&a, &e, n / 2);
	double *ce = (double *)malloc(n * n / 4 * sizeof(double));
	strassen_matmat(&c, &e, &ce, n / 2, n / 2, n / 2);

	double *temp1 = (double *)malloc(n * n / 4 * sizeof(double));
	strassen_matmat(&ce, &b, &temp1, n / 2, n / 2, n / 2);
	double *Z = darray_add(d, temp1, n * n / 4, -1.0);
	strassen_invert_strassen_matmat(&Z, &t, n / 2);

	double *temp2 = (double *)malloc(n * n / 4 * sizeof(double));
	double *ebt = (double *)malloc(n * n / 4 * sizeof(double));
	strassen_matmat(&e, &b, &temp1, n / 2, n / 2, n / 2);
	strassen_matmat(&temp1, &t, &ebt, n / 2, n / 2, n / 2);
	strassen_matmat(&ebt, &ce, &temp2, n / 2, n / 2, n / 2);
	strassen_matmat(&t, &ce, &temp1, n / 2, n / 2, n / 2);

	// Calculate the inverse A blocks
	// inverse_A11
	mat_inplace_block_add(*inverse_A, e, temp2, start_a, n, n, 1.0, 1.0);
	// inverse_A12
	mat_inplace_block_add(*inverse_A, ebt, NULL, start_b, n, n, -1.0, 0.0);
	// inverse_A21
	mat_inplace_block_add(*inverse_A, temp1, NULL, start_c, n, n, -1.0,
			      0.0);
	// inverse_A22
	mat_inplace_block_add(*inverse_A, t, NULL, start_d, n, n, 1.0, 0.0);

	free(a);
	free(b);
	free(c);
	free(d);
	free(e);
	free(t);
	free(temp1);
	free(temp2);
	free(ce);
	free(ebt);
	free(Z);
	// Depad
	id_depad_matrix(A, n, og_n);
	id_depad_matrix(inverse_A, n, og_n);
}

void strassen_invert_naive_matmat(double *A, double *inverse_A, size_t n) {
	if (n == 1) {
		inverse_A[0] = 1 / A[0];
		return;
	}
	const size_t start_a = 0;
	const size_t start_b = n / 2;
	const size_t start_c = n / 2 * n;
	const size_t start_d = start_c + start_b;
	double *a = create_block(A, start_a, n, n);
	double *b = create_block(A, start_b, n, n);
	double *c = create_block(A, start_c, n, n);
	double *d = create_block(A, start_d, n, n);
	double *e = (double *)calloc(n * n / 4, sizeof(double));
	double *t = (double *)calloc(n * n / 4, sizeof(double));
	strassen_invert_naive_matmat(a, e, n / 2);
	double *ce = (double *)malloc(n * n / 4 * sizeof(double));
	naive_matmat(c, e, ce, n / 2, n / 2, n / 2);
	double *temp1 = (double *)malloc(n * n / 4 * sizeof(double));
	naive_matmat(ce, b, temp1, n / 2, n / 2, n / 2);
	double *Z = darray_add(d, temp1, n * n / 4, -1.0);
	strassen_invert_naive_matmat(Z, t, n / 2);

	double *temp2 = (double *)malloc(n * n / 4 * sizeof(double));
	double *ebt = (double *)malloc(n * n / 4 * sizeof(double));
	naive_matmat(e, b, temp1, n / 2, n / 2, n / 2);
	naive_matmat(temp1, t, ebt, n / 2, n / 2, n / 2);
	naive_matmat(ebt, ce, temp2, n / 2, n / 2, n / 2);
	naive_matmat(t, ce, temp1, n / 2, n / 2, n / 2);

	// Calculate the inverse A blocks
	// inverse_A11
	mat_inplace_block_add(inverse_A, e, temp2, start_a, n, n, 1.0, 1.0);
	// inverse_A12
	mat_inplace_block_add(inverse_A, ebt, NULL, start_b, n, n, -1.0, 0.0);
	// inverse_A21
	mat_inplace_block_add(inverse_A, temp1, NULL, start_c, n, n, -1.0, 0.0);
	// inverse_A22

	mat_inplace_block_add(inverse_A, t, NULL, start_d, n, n, 1.0, 0.0);

	free(a);
	free(b);
	free(c);
	free(d);
	free(e);
	free(t);
	free(temp1);
	free(temp2);
	free(ce);
	free(ebt);
	free(Z);
}

/*void strassen_invert_naive_matmat(double *A, double *inverse_A,
				  const size_t n) {
	if (n == 1) {
		inverse_A[0] = 1 / A[0];
		return;
	}
	const size_t start_a = 0;
	const size_t start_b = n / 2;
	const size_t start_c = n / 2 * n;
	const size_t start_d = start_c + start_b;
	double *a = create_block(A, start_a, n, n);
	double *b = create_block(A, start_b, n, n);
	double *c = create_block(A, start_c, n, n);
	double *d = create_block(A, start_d, n, n);
	double *e = (double *)malloc(n * n / 4 * sizeof(double));
	double *t = (double *)malloc(n * n / 4 * sizeof(double));
	strassen_invert_naive_matmat(a, e, n / 2);
	double *ce = (double *)malloc(n * n / 4 * sizeof(double));
	naive_matmat(c, e, ce, n / 2, n / 2, n / 2);

	double *temp1 = (double *)malloc(n * n / 4 * sizeof(double));
	naive_matmat(ce, b, temp1, n / 2, n / 2, n / 2);
	double *Z = darray_add(d, temp1, n * n / 4, -1.0);
	strassen_invert_naive_matmat(Z, t, n / 2);
	double *temp2 = (double *)malloc(n * n / 4 * sizeof(double));
	double *ebt = (double *)malloc(n * n / 4 * sizeof(double));
	naive_matmat(e, b, temp1, n / 2, n / 2, n / 2);
	naive_matmat(temp1, t, ebt, n / 2, n / 2, n / 2);
	naive_matmat(ebt, ce, temp2, n / 2, n / 2, n / 2);
	naive_matmat(t, ce, temp1, n / 2, n / 2, n / 2);

	// Calculate the inverse A blocks
	// inverse_A11
	mat_inplace_block_add(inverse_A, e, temp2, start_a, n, n, 1.0, 1.0);
	// inverse_A12 TODO: add a codition where 0 is second arg.
	mat_inplace_block_add(inverse_A, ebt, ebt, start_b, n, n, -1.0, 0.0);
	// inverse_A21
	mat_inplace_block_add(inverse_A, temp1, temp1, start_c, n, n, -1.0,
			      0.0);
	// inverse_A22
	mat_inplace_block_add(inverse_A, t, t, start_d, n, n, 1.0, 0.0);

	free(a);
	free(b);
	free(c);
	free(d);
	free(e);
	free(t);
	free(temp1);
	free(temp2);
	free(ce);
	free(ebt);
	free(Z);
} */

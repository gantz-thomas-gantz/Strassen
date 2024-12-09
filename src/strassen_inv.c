/*
 * DESC: Module for naive matrix multiplication.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */

#include <stdlib.h>

#include "../include/block_utilities.h"
#include "../include/naive_matmat.h"
#include "../include/strassen_matmat.h"
// TODO: Make identity padding

void strassen_invert_strassen_matmat(double *A, double *inverse_A,
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
	strassen_invert_strassen_matmat(a, e, n / 2);
	double *ce = (double *)malloc(n * n / 4 * sizeof(double));
	strassen_matmat(&c, &e, &ce, n / 2, n / 2, n / 2);

	double *temp1 = (double *)malloc(n * n / 4 * sizeof(double));
	strassen_matmat(&ce, &b, &temp1, n / 2, n / 2, n / 2);
	double *Z = darray_add(d, temp1, n * n / 4, -1.0);
	strassen_invert_strassen_matmat(Z, t, n / 2);
	double *temp2 = (double *)malloc(n * n / 4 * sizeof(double));
	double *ebt = (double *)malloc(n * n / 4 * sizeof(double));
	strassen_matmat(&e, &b, &temp1, n / 2, n / 2, n / 2);
	strassen_matmat(&temp1, &t, &ebt, n / 2, n / 2, n / 2);
	strassen_matmat(&ebt, &ce, &temp2, n / 2, n / 2, n / 2);
	strassen_matmat(&t, &ce, &temp1, n / 2, n / 2, n / 2);

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
}

void strassen_invert_naive_matmat(double *A, double *inverse_A,
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
	printf("Start %lf and End %lf \n", (double)start, (double)end);

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
}

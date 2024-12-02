/*
 * DESC: Module for naive matrix multiplication.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */
#include <assert.h>
#include <stddef.h>  // for size_t
#include <stdio.h>
#include <stdlib.h>

double *darray_add(const double *const A, const double *const B,
		   const size_t size, const double alpha) {
	double *C = (double *)malloc(size * sizeof(double));
	for (size_t i = 0; i < size; i++) C[i] = A[i] + alpha * B[i];
	return C;
}

double *darray_block_add(const double *const A, const size_t start1,
			 const size_t start2, const size_t m, const size_t n,
			 const double alpha) {
	assert((start1 + (m / 2 - 1) * n + (n / 2 - 1) < m * n) ||
	       (start2 + (m / 2 - 1) * n + (n / 2 - 1) < m * n));
	double *C = (double *)malloc(m / 2 * n / 2 * sizeof(double));
	for (size_t i = 0; i < m / 2; i++)
		for (size_t j = 0; j < n / 2; j++)
			C[i * n / 2 + j] = A[start1 + i * n + j] +
					   alpha * A[start2 + i * n + j];
	return C;
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

void mat_inplace_block_add(double *C, const double *const a, const size_t start,
			   const size_t m, const size_t n, const double alpha) {
	assert((start + (m / 2 - 1) * n + (n / 2 - 1) < m * n));
	for (size_t i = 0; i < m / 2; i++)
		for (size_t j = 0; j < n / 2; j++)
			C[start + i * n + j] += alpha * a[i * n / 2 + j];
}

// TODO: Maybe instead of creating a new matrix, simply realloc the old one and
// memcpy the values correctly within new memory to account for padding?
void pad_matrix(double **A, size_t m, size_t n, size_t *new_m, size_t *new_n) {
	*new_m = (m % 2 == 0) ? m : m + 1;
	*new_n = (n % 2 == 0) ? n : n + 1;
	if (!(m % 2 == 0) || !(n % 2 == 0)) {
		// Calculate new dimensions
		// Allocate memory for the new padded matrix
		double *padded_A =
		    (double *)calloc((*new_m) * (*new_n), sizeof(double));
		// Copy elements from the original matrix to the padded matrix
		for (size_t i = 0; i < m; i++)
			for (size_t j = 0; j < n; j++)
				padded_A[i * (*new_n) + j] = (*A)[i * n + j];
		free(*A);
		*A = padded_A;
	}
}

void depad_matrix(double **padded_A, size_t m, size_t n, size_t og_m,
		  size_t og_n) {
	// If the padded dimensions match the original, no need to depad
	if (m == og_m && n == og_n) {
		return;
	}
	// Allocate memory for the original matrix
	double *A = (double *)malloc(og_m * og_n * sizeof(double));
	// Copy elements from the padded matrix to the original matrix
	for (size_t i = 0; i < og_m; i++) {
		for (size_t j = 0; j < og_n; j++) {
			A[i * og_n + j] = (*padded_A)[i * n + j];
		}
	}
	free(*padded_A);
	*padded_A = A;
}

// maybe add recursive func & non-recursive func
// A, B & C: data doesn't change but location changes
// m, n, k change throughout but are the same at the end of recursion
// this is reason for not having constant args ()
void strassen_matmat(double **A, double **B, double **C, size_t m, size_t n,
		     size_t k) {
	if (m == 1 && n == 1) {
		for (size_t i = 0; i < m * k; i++) (*C)[i] = (*A)[0] * (*B)[i];
	} else if (n == 1 && k == 1) {
		for (size_t i = 0; i < m * k; i++) (*C)[i] = (*A)[i] * (*B)[0];
	} else {
		for (size_t i = 0; i < m * k; i++) (*C)[i] = 0;

		// Pad matrices if needed
		size_t og_m = m;
		size_t og_n = n;
		size_t og_k = k;
		pad_matrix(A, og_m, og_n, &m, &n);
		pad_matrix(B, og_n, og_k, &n, &k);
		pad_matrix(C, og_m, og_k, &m, &k);
		// Define block start indices
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
		// Deep copy block matrices into contiguous memory
		// TODO: Instead of creating blocks for adding their values
		double *a = create_block(*A, start_a, m, n);
		double *d = create_block(*A, start_d, m, n);
		double *y = create_block(*B, start_y, n, k);
		double *z = create_block(*B, start_z, n, k);

		// Compute q coefficients
		// TODO: We keep instantiating new temps, although we just need
		// 2. Also, instead of having darray_block_add return an array,
		// pass in the temporary memory location to store the
		// intermediate results.
		double *q1 = (double *)malloc(m * k / 4 * sizeof(double));
		double *temp1 =
		    darray_block_add(*B, start_x, start_z, n, k, 1.0);
		strassen_matmat(&a, &temp1, &q1, m / 2, n / 2, k / 2);
		free(temp1);  // Free after use
		free(a);

		double *q2 = (double *)malloc(m * k / 4 * sizeof(double));
		double *temp2 =
		    darray_block_add(*B, start_y, start_t, n, k, 1.0);
		strassen_matmat(&d, &temp2, &q2, m / 2, n / 2, k / 2);
		free(temp2);
		free(d);

		double *q3 = (double *)malloc(m * k / 4 * sizeof(double));
		double *temp3_1 =
		    darray_block_add(*A, start_d, start_a, m, n, -1.0);
		double *temp3_2 =
		    darray_block_add(*B, start_z, start_y, n, k, -1.0);
		strassen_matmat(&temp3_1, &temp3_2, &q3, m / 2, n / 2, k / 2);
		free(temp3_1);
		free(temp3_2);

		double *q4 = (double *)malloc(m * k / 4 * sizeof(double));
		double *temp4_1 =
		    darray_block_add(*A, start_b, start_d, m, n, -1.0);
		double *temp4_2 =
		    darray_block_add(*B, start_z, start_t, n, k, 1.0);
		strassen_matmat(&temp4_1, &temp4_2, &q4, m / 2, n / 2, k / 2);
		free(temp4_1);
		free(temp4_2);

		double *q5 = (double *)malloc(m * k / 4 * sizeof(double));
		double *temp5 =
		    darray_block_add(*A, start_b, start_a, m, n, -1.0);
		strassen_matmat(&temp5, &z, &q5, m / 2, n / 2, k / 2);
		free(temp5);
		free(z);

		double *q6 = (double *)malloc(m * k / 4 * sizeof(double));
		double *temp6_1 =
		    darray_block_add(*A, start_c, start_a, m, n, -1.0);
		double *temp6_2 =
		    darray_block_add(*B, start_x, start_y, n, k, 1.0);
		strassen_matmat(&temp6_1, &temp6_2, &q6, m / 2, n / 2, k / 2);
		free(temp6_1);
		free(temp6_2);

		double *q7 = (double *)malloc(m * k / 4 * sizeof(double));
		double *temp7 =
		    darray_block_add(*A, start_c, start_d, m, n, -1.0);
		strassen_matmat(&temp7, &y, &q7, m / 2, n / 2, k / 2);
		free(temp7);
		free(y);
		// Calculate the R blocks
		// TODO: Instead of having darray_block_add return an array,
		// pass in the temporary memory location to store the
		// intermediate result. This will result in less additions here.
		// R11
		mat_inplace_block_add(*C, q1, start_r11, m, k, 1.0);
		mat_inplace_block_add(*C, q5, start_r11, m, k, 1.0);
		// R21
		mat_inplace_block_add(*C, q1, start_r21, m, k, 1.0);
		mat_inplace_block_add(*C, q3, start_r21, m, k, 1.0);
		mat_inplace_block_add(*C, q6, start_r21, m, k, 1.0);
		mat_inplace_block_add(*C, q7, start_r21, m, k, -1.0);
		// R12
		mat_inplace_block_add(*C, q2, start_r12, m, k, 1.0);
		mat_inplace_block_add(*C, q3, start_r12, m, k, 1.0);
		mat_inplace_block_add(*C, q4, start_r12, m, k, 1.0);
		mat_inplace_block_add(*C, q5, start_r12, m, k, -1.0);
		// R22
		mat_inplace_block_add(*C, q2, start_r22, m, k, 1.0);
		mat_inplace_block_add(*C, q7, start_r22, m, k, 1.0);
		// Write back to C
		free(q1);
		free(q2);
		free(q3);
		free(q4);
		free(q5);
		free(q6);
		free(q7);

		// DEPAD HERE
		depad_matrix(A, m, n, og_m, og_n);
		depad_matrix(B, n, k, og_n, og_k);
		depad_matrix(C, m, k, og_m, og_k);
	}
}


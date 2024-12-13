/*
 * DESC: Header of module for naive matrix multiplication.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */

#include <stddef.h>

/*
 * Description:
 * Multiply A (size mxn) with B (size nxk) naively, store the result in C (size
 * mxk).
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
void naive_matmat(double *A, double *B, double *C, const size_t m,
		  const size_t n, const size_t k);


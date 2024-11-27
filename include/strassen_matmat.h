/*
 * DESC: Header of module for strassen matrix multiplication.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */

#include <stddef.h>  // for size_t

/*
 * Description:
 * Multiply A (size mxn) with B (size nxk) using Strassen's multiplication
 * algorithm, store the result in C (size mxk).
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
void strassen_matmat(double **A, double **B, double **C, size_t m, size_t n,
		     size_t k);


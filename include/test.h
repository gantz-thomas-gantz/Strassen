/*
 * DESC: Header of module for testing all implementations.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */

#include <stddef.h>  // for size_t

/*
 * Description:
 * Call naive_matmat implementation and time, also compare to CBLAS to assert
 * correctness of result.
 *
 * Return:
 * time in seconds. If -1, wrong result.
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
double test_naive_matmat(double *A, double *B, double *C, const size_t m,
			 const size_t n, const size_t k);


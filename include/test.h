/*
 * DESC: Header of module for testing all implementations.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */

#include <stddef.h>  // for size_t

/*
 * Description:
 * Compare two double matrices.
 *
 * Return:
 * If 0, wrong result, if 1 correct.
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
int compare_mat(const double *const A, const double *const B, const size_t m,
		const size_t n);

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
double test_naive_matmat(const size_t N, const double eps);


/*
 * DESC: Header of module for testing all implementations.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */

#include <stddef.h>  // for size_t

void flush_cache();

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

void gen_rand_matrix(double *A, const size_t m, const size_t n);
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

/*
 * Description:
 * Call strassen_matmat implementation and time, also compare to CBLAS to assert
 * correctness of result.
 *
 * Return:
 * time in seconds. If -1, wrong result.
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
double test_strassen_matmat(const size_t N, const double eps);

int is_invertible(double *A, int n);

// note: reason for strassen_inverse_strassen_matmat is a bit slower for small
// matrices is that strassen_matmat calls on naive_matmat function. overhead !
double test_strassen_invert_strassen_matmat(double **A, const size_t n,
					    const double eps);

double test_strassen_invert_naive_matmat(double **A, const size_t n,
					 const double eps);

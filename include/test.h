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

/*
 * Description:
 * Generate a random double matrix.
 *
 * Arguments:
 * - `A`: Pointer to the array where the generated matrix will be stored.
 * - `m`: Number of rows in the matrix.
 * - `n`: Number of columns in the matrix.
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
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
double test_naive_matmat(double **A, double **B, const size_t m, const size_t n, const size_t k, const double eps);

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
double test_strassen_matmat(double **A, double **B, const size_t m, const size_t n, const size_t k, const double eps);

/*
 * Description:
 * Check if a square matrix is invertible using LAPACK's LU inversion function.
 *
 * Arguments:
 * - `A`: Pointer to the matrix.
 * - `n`: Dimension of the square matrix.
 *
 * Return:
 * 1 if the matrix is invertible, 0 otherwise.
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
int is_invertible(double *A, int n);

/*
 * Description:
 * Test the Strassen inversion algorithm implementation.
 * Compares the result to LAPACK's expected output to validate
 * correctness.
 *
 * Arguments:
 * - `A`: Pointer to the matrix.
 * - `n`: Dimension of the square matrix.
 * - `eps`: Tolerance for comparison.
 *
 * Return:
 * Time in seconds. If -1, wrong result.
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
double test_strassen_invert_strassen_matmat(double **A, const size_t n,
					    const double eps);

/*
 * Description:
 * Test the naive block inversion algorithm implementation.
 * Compares the result to LAPACK's output to validate
 * correctness.
 *
 * Arguments:
 * - `A`: Pointer to the matrix.
 * - `n`: Dimension of the square matrix.
 * - `eps`: Tolerance for comparison.
 *
 * Return:
 * Time in seconds. If -1, wrong result.
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
double test_strassen_invert_naive_matmat(double **A, const size_t n,
					 const double eps);

/*
 * Description:
 * Test the implementation of LU decomposition inversion of a matrix.
 * Compares the result to LAPACK's output to validate correctness.
 *
 * Arguments:
 * - `A`: Pointer to the matrix to be inverted.
 * - `n`: Dimension of the square matrix.
 * - `eps`: Tolerance for comparison.
 *
 * Return:
 * Time in seconds. If -1, wrong result.
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
double test_lu_invert(const double *const A, const size_t n, const double eps);


#include <stddef.h>  // for size_t

/*
 * Description:
 * Perform element-wise addition of two arrays, scaling the second array
 * by a scalar factor, and store the result in a newly allocated array.
 *
 * Arguments:
 * - `A`: Pointer to the first input array.
 * - `B`: Pointer to the second input array.
 * - `size`: Number of elements in each array.
 * - `alpha`: Scalar to multiply elements of array `B`.
 *
 * Return:
 * Pointer to a newly allocated array containing the result of A + alpha * B.
 */
double *darray_add(const double *const A, const double *const B,
		   const size_t size, const double alpha);

/*
 * Description:
 * Perform element-wise addition of submatrices (blocks) of a 2D mxn matrix,
 * scaling block starting at start2 by a scalar factor, and store the result in
 * a specified block of the output matrix.
 *
 * Arguments:
 * - `A`: Pointer to the input matrix.
 * - `C`: Pointer to the output matrix.
 * - `start1`: Starting index of the first block in matrix `A`.
 * - `start2`: Starting index of the second block in matrix `A`.
 * - `m`: Number of rows in the original matrix.
 * - `n`: Number of columns in the original matrix.
 * - `alpha`: Scalar to multiply the elements of the second block.
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
void darray_block_add(const double *const A, double *C, const size_t start1,
		      const size_t start2, const size_t m, const size_t n,
		      const double alpha);

/*
 * Description:
 * Extract a submatrix (block) from A starting at a specified
 * position (start). The block size is half the number of rows and columns of
 * the original matrix.
 *
 * Arguments:
 * - `A`: Pointer to the input matrix.
 * - `start`: Starting index of the block in matrix `A`.
 * - `m`: Number of rows in A.
 * - `n`: Number of columns in A.
 *
 * Return:
 * Pointer to a newly allocated array containing the extracted block.
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
double *create_block(const double *const A, const size_t start, const size_t m,
		     const size_t n);

/*
 * Description:
 * Perform in-place addition of submatrices (blocks) into a specified block
 * of a matrix. Optionally, scale the input blocks by scalar factors before
 * adding.
 *
 * Arguments:
 * - `C`: Pointer to the output matrix where the result will be stored.
 * - `a`: Pointer to the first input block.
 * - `b`: Pointer to the second input block (can be `NULL` in case none to be
 *   added).
 * - `start`: Starting index of the block in matrix `C` to update.
 * - `m`: Number of rows in the original matrix.
 * - `n`: Number of columns in the original matrix.
 * - `alpha`: Scalar to multiply the elements of block `a`.
 * - `beta`: Scalar to multiply the elements of block `b` (if `b` is not
 * `NULL`).
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
void mat_inplace_block_add(double *C, const double *const a,
			   const double *const b, const size_t start,
			   const size_t m, const size_t n, const double alpha,
			   const double beta);


/*
 * DESC: Module for Input/Output functions.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */
#include <stdbool.h>
#include <stddef.h>  // for size_t

/*
 * Description:
 * Print the elements of a matrix of type `double` to the standard output in
 * a readable row-major format.
 *
 * Arguments:
 * - `M`: Pointer to the matrix to be printed.
 * - `m`: Number of rows in the matrix.
 * - `n`: Number of columns in the matrix.
 *
 * Matrix format:
 * The matrix should be a flattened array in row-major format.
 */
void print_mat_double(double *M, const size_t m, const size_t n);


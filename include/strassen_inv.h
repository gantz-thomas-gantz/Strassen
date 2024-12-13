/*
 * DESC: Header of module for Strassen invert.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */

#include <stddef.h>

/*
 * Description:
 * Invert A (size nxn) using recursive block inversion and strassen
 * multiplication.
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
void strassen_invert_naive_matmat(double **A, double **inverse_A, size_t n);

/*
 * Description:
 * Invert A (size nxn) using recursive block inversion and naive
 * multiplication.
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
void strassen_invert_strassen_matmat(double **A, double **inverse_A, size_t n);

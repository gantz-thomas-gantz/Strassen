/*
 * DESC: Header of module for Strassen invert.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */

#include <stddef.h>  // for size_t

/*
 * Description:
 * Invert A (size nxn), by ...
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
void strassen_invert_naive_matmat(double **A, double **inverse_A, size_t n);

/*
 * Description:
 * Invert A (size nxn), by ...
 *
 * Matrix format:
 * Matrices should be flattened arrays in row-major format.
 */
void strassen_invert_strassen_matmat(double **A, double **inverse_A, size_t n);

/*
 * DESC: Input/Output utilities module.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */
#include "IO.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void print_mat_double(double *M, const size_t m, const size_t n) {
	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < n; j++) {
			printf("%lf\t", M[i * n + j]);
		}
		printf("\n");
	}
}


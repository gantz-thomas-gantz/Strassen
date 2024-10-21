#include <stdbool.h>
#include <stddef.h>  // for size_t

int read_matrix(const char *filename, double **A, size_t *m, size_t *n);

void print_mat_double(double *M, const size_t m, const size_t n);

#include "IO.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Function to read matrix from file and store it in a flattened array
int read_matrix(const char *filename, double **A, size_t *m, size_t *n) {
	FILE *file = fopen(filename, "r");
	// Read the size of the matrix
	fscanf(file, "%zu %zu", m, n);
	fgetc(file);
	// Allocate memory for the flattened matrix
	*A = (double *)malloc((*m) * (*n) * sizeof(double));
	// Read the matrix elements row by row
	char buffer[1024];
	for (int i = 0; i < *m; i++) {
		fgets(buffer, sizeof(buffer), file);
		// Split the line into columns using semicolon as delimiter
		char *token = strtok(buffer, " ");
		int j = 0;
		while (token != NULL) {
			(*A)[i * (*n) + j++] = atof(token);
			token = strtok(NULL, " ");
		}
	}
	fclose(file);
	return 0;
}
// Function to print a flattened matrix for 'double' type
void print_mat_double(double *M, const size_t m, const size_t n) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			printf("%lf\t", M[i * n + j]);
		}
		printf("\n");
	}
}


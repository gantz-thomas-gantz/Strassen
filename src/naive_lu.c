/*
 * DESC: Module for LU decomposition.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */

// Function to initialize T to A
static void init_T(const double *A, double *T, const size_t n) {
	memcpy(T, A, sizeof(double) * n * n);
}

static void eliminate_step_i(double *T, const int n, const int i) {
	for (int j = i + 1; j < n; j++) {
		T[j * n + i] = T[j * n + i] / T[i * n + i];
		for (int k = i + 1; k < n; k++) {
			T[j * n + k] =
			    T[j * n + k] - T[j * n + i] * T[i * n + k];
		}
	}
}

void lu_decomposition(const double *const A, double *T, size_t n) {
	// L and U
	init_T(A, T, n);
	// Comput L and U
	printf("Printing A \n");
	print_mat_double(A, n);
	for (int i = 0; i < n - 1; i++) {
		printf("Step %d:\n", i);
		printf("Printing T before elimination\n");
		print_mat_double(T, n);
		eliminate_step_i(T, n, i);
		printf("Printing T after elimination\n");
		print_mat_double(T, n);
		printf("\n");
	}
}

void lu_invert(const double *const A, double *inverse_A, size_t n) {
	// Allocate space for the LU matrix
	double *T = (double *)malloc(sizeof(double) * n * n);
	// Perform LU decomposition
	lu_decomposition(A, T, n);
	// Temporary vector for storing canonical vectors
	double *temp = (double *)malloc(sizeof(double) * n);
	// Solve for each column of the inverse matrix
	for (size_t col = 0; col < n; col++) {
		// Create the identity column (current column of identity
		// matrix)
		for (size_t i = 0; i < n; i++) {
			temp[i] = (i == col) ? 1.0 : 0.0;
		}
		// Solve L * Y = ei
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < i; j++) {
				temp[i] -= T[i * n + j] * temp[j];
			}
		}
		// Solve U * X = Y
		for (int i = n - 1; i >= 0; i--) {
			for (size_t j = i + 1; j < n; j++) {
				temp[i] -= T[i * n + j] * temp[j];
			}
			temp[i] /= T[i * n + i];  // Diagonal element of U
		}
		// Write the solved column into the inverse matrix
		for (size_t i = 0; i < n; i++) {
			inverse_A[i * n + col] = temp[i];
		}
	}
	free(T);
	free(temp);
}


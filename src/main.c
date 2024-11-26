/*
 * DESC: Module for naive matrix multiplication.
 * AUTHORS: Thomas Gantz, Laura Paxton, Jan Marxen
 */
#include <math.h>  // for pow()
#include <stdio.h>
#include <stdlib.h>

#include "../include/IO.h"
#include "../include/test.h"

int main(int argc, char *argv[]) {
	const size_t N = 10;  // max power dimension of matrix
	for (size_t i = 1; i < N; i++) {
		printf("Size of test: %zu \n", i);
		printf("test_naive_matmat time: %lf \n",
		       test_naive_matmat(i, 0.001));
		printf("test_strassen_matmat time: %lf \n",
		       test_strassen_matmat(i, 0.001));
	}
}

#include <stddef.h>  // for size_t
// TODO: Create function descriptions

double *darray_add(const double *const A, const double *const B,
		   const size_t size, const double alpha);

void darray_block_add(const double *const A, double *C, const size_t start1,
		      const size_t start2, const size_t m, const size_t n,
		      const double alpha);

double *create_block(const double *const A, const size_t start, const size_t m,
		     const size_t n);

void mat_inplace_block_add(double *C, const double *const a,
			   const double *const b, const size_t start,
			   const size_t m, const size_t n, const double alpha,
			   const double beta);


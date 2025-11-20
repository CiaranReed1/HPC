#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <limits.h>

void init(double *a, double *b, size_t N){
	/* Intialise with some values */
	// Parallelize this loop!
	for (size_t i = 0; i < N; i++) {
		a[i] = i+1;
		b[i] = (-1)*(i%2) * i;
	}
}

double dot(double *a, double *b, size_t N){
	double dotabparallel = 0;
	double start = omp_get_wtime();
	{
	   /* Implement a parallel dot product using
	   *
	   * 1. The manual approach,
	   * 2. critical sections to protect the shared updates,
	   * 3. atomics to protect the shared updates,
	   * 4. the reduction clause.
	   */
	   for (size_t i = 0; i < N; i++) {
	   }
	}
  	return dotabparallel;
}

int main(int argc, char **argv)
{
  const size_t N = argc > 1 ? strtoul(argv[1], NULL, 10) : 1024u;
  double *a = malloc(N * sizeof(*a));
  double *b = malloc(N * sizeof(*b));
  double adotb;
  init(a, b, N);
  double start = omp_get_wtime();
  adotb = dot(a, b, N);
  printf("%d\t%e\t%lu\t%.4g\n", omp_get_max_threads(), omp_get_wtime()-start, N, adotb);
  free(a);
  free(b);
  return 0;
}

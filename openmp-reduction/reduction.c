#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <limits.h>


void init(double *a, double *b, size_t N){
	/* Intialise with some values */
	// Parallelize this loop!
	#pragma omp parallel for
	for (size_t i = 0; i < N; i++) {
		a[i] = i+1;
		b[i] = (-1)*(i%2) * i;
	}
}

double dot_manual(double *a, double *b, size_t N){
	double dotabparallel = 0;
	int partial_sums[omp_get_num_threads()];
	#pragma omp parallel for
	for (size_t i = 0; i < N; i++) {
		partial_sums[omp_get_thread_num()] += a[i] * b[i];
	}
	for (int i = 0; i < omp_get_num_threads(); i++) {
		dotabparallel += partial_sums[i];
	}
  	return dotabparallel;
}

double dot_critical(double *a, double *bm size_t N){
	double dotabparallel = 0;
	#pragma omp parallel for
	for (size_t i = 0; i < N; i++) {
		#pragma omp critical
		{
			dotabparallel += a[i] * b[i];
		}
	}
  	return dotabparallel;
}

double dot_atomic(double *a, double *b, size_t N){
	double dotabparallel[N];
	double dottotal = 0;
	#pragma omp parallel for
	for (size_t i = 0; i < N; i++) {
		dotabparallel[i] = a[i] * b[i];
		#pragma omp atomic
		dottotal += dotabparallel[i];
	}
  	return dottotal;
}

int main(int argc, char **argv)
{
  const size_t N = argc > 1 ? strtoul(argv[1], NULL, 10) : 1024u;
  double *a = malloc(N * sizeof(*a));
  double *b = malloc(N * sizeof(*b));
  double adotb;
  init(a, b, N);
  omp_set_num_threads(8);
  double start = omp_get_wtime();
  adotb = dot(a, b, N);
  printf("%d\t%e\t%lu\t%.4g\n", omp_get_max_threads(), omp_get_wtime()-start, N, adotb);
  free(a);
  free(b);
  return 0;
}

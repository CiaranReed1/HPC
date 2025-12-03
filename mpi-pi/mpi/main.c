
// This file is part of the HPC workshop 2019 at Durham University
// Author: Christian Arnold

/*
 * This program calculates the numerical value of the constant PI.
 */

#include "proto.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#ifndef M_PI
#define M_PI 3.141592653589793
#endif
/*
 * The program takes one argument: the number of samples to use.
 */
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("Hello from process %d of %d\n", rank, size);

  double start, end;

  int N; /* the total count of random coordinates vectors */
  if (argc != 2) {
    fprintf(stderr, "Usage: %s N\n", argv[0]);
    fprintf(stderr, "\nEstimate PI by Monte-Carlo using N samples\n");
    return 1;
  }
  N = atoi(argv[1]);

  int K = 1000; // number of independent runs
  double error_sum = 0.0;
  double runtime_sum = 0.0;
  double pi_est_sum = 0.0;

  for(int k = 0; k < K; k++) {
      start = MPI_Wtime();
      double pi_est = calculate_pi(N, k+1*(rank+1)); // varying seed
      end = MPI_Wtime();
      pi_est_sum += pi_est;
      runtime_sum += (end - start);
      double err = pi_est - M_PI;
      error_sum += err*err;
  }
  double my_pi = pi_est_sum / K;
  double rms_error = sqrt(error_sum / K);
  double avg_runtime = runtime_sum / K;
  printf("%d,%f,%f,%.6e\n", N, my_pi, rms_error, avg_runtime);
  MPI_Finalize();
  return 0;
}


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
#ifndef M_PI
#define M_PI 3.141592653589793
#endif
/*
 * The program takes one argument: the number of samples to use.
 */
int main(int argc, char *argv[]) {
  clock_t start, end;

  int N; /* the total count of random coordinates vectors */
  if (argc != 2) {
    fprintf(stderr, "Usage: %s N\n", argv[0]);
    fprintf(stderr, "\nEstimate PI by Monte-Carlo using N samples\n");
    return 1;
  }
  N = atoi(argv[1]);

  int nIter = 1000;
  double my_pi_total = 0.0;
  double time_total = 0.0;
  
for (int i = 0; i < nIter; i++){
  start = clock();
  my_pi_total += calculate_pi(N,i);
  end = clock();
  time_total += ((double)end - start) / CLOCKS_PER_SEC;
}


  double my_pi = my_pi_total / nIter;
  double time = time_total / nIter;
  double error = M_PI - my_pi;
  printf("%d, %.20f,%.20f,%g\n",N, my_pi, error,time);
 

  return 0;
}

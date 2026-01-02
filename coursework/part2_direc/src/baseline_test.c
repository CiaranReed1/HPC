#include "params.h" // model & simulation parameters
#include <math.h>   // needed for tanh, used in init function
#include <stdio.h>  // needed for printing
#include <stdlib.h> // needed for malloc and free
#include <mpi.h>  // MPI header

int M = 1024;  // x domain size
int N = 512;   // y domain size

void init(double *u, double *v) {
  int idx;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      idx = i * N + j;
      u[idx] = ulo + (uhi - ulo) * 0.5 * (1.0 + tanh((i - M / 2) / 16.0));
      v[idx] = vlo + (vhi - vlo) * 0.5 * (1.0 + tanh((j - N / 2) / 16.0));
    }
  }
}

void dxdt(double *du, double *dv, const double *u, const double *v) {
  double lapu, lapv, _u, _v;
  int up, down, left, right, idx;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      idx = i * N + j;
      if (i == 0) {
        down = i * N + j;
      } else {
        down = (i - 1) * N + j;
      }
      if (i == M - 1) {
        up = i * N + j;
      } else {
        up = (i + 1) * N + j;
      }
      if (j == 0) {
        left = i * N + j;
      } else {
        left = i * N + j - 1;
      }
      if (j == N - 1) {
        right = i * N + j;
      } else {
        right = i * N + j + 1;
      }
      _u = u[idx];
      _v = v[idx];
      lapu = u[up] + u[down] + u[left] + u[right] + -4.0 * _u;
      lapv = v[up] + v[down] + v[left] + v[right] + -4.0 * _v;
      du[idx] = DD * lapu + f(_u, _v) + R * stim(i, j);
      dv[idx] = d * DD * lapv + g(_u, _v);
    }
  }
}

void step(const double *du, const double *dv, double *u, double *v) {
  int idx;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      idx = i * N + j;
      u[idx] += dt * du[idx];
      v[idx] += dt * dv[idx];
    }
  }
}

double norm(const double *x) {
  double nrmx = 0.0;
  int idx;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      idx = i * N + j;
      nrmx += x[idx] * x[idx];
    }
  }
  return nrmx;
}

int main(int argc, char **argv) {

  int problem_size = argc > 1 ? atoi(argv[1]) : 1;
  int nruns = argc > 2 ? atoi(argv[2]) : 1;

  double scale = sqrt(problem_size);
  M = (int)(M * scale + 0.5);  // round to nearest int
  N = (int)(N * scale + 0.5);
  M += M % 2;  // add 1 if odd
  N += N % 2;  // add 1 if odd

  MPI_Init(&argc, &argv);
  double t = 0.0, nrmu, nrmv;
  int writeInd = 0;
  double stats[T / m][3];


  // Allocate memory for 1D arrays representing 2D grids
  double *u = (double *)malloc(M * N * sizeof(double));
  double *v = (double *)malloc(M * N * sizeof(double));
  double *du = (double *)malloc(M * N * sizeof(double));
  double *dv = (double *)malloc(M * N * sizeof(double));

  double init_time = 0.0, step_time = 0.0, norm_time = 0.0, dxdt_time = 0.0;
  double start_time, end_time, temp_start, temp_end;


 

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // Check for allocation success
  if (!u || !v || !du || !dv) {
    fprintf(stderr, "Error: Failed to allocate memory\n");
    return 1;
  }
  // initialize the state
  for (int run = 0; run < nruns; run++) {
    init_time = 0.0;
    step_time = 0.0;
    norm_time = 0.0;
    dxdt_time = 0.0;
    start_time = MPI_Wtime();
    temp_start = MPI_Wtime();
    init(u, v);
    temp_end = MPI_Wtime();
    init_time += (temp_end - temp_start);
    // time-loop
    for (int k = 0; k < T; k++) {
      // track the time
      t = dt * k;
      // evaluate the PDE
      temp_start = MPI_Wtime();
      dxdt(du, dv, u, v);
      temp_end = MPI_Wtime();
      dxdt_time += (temp_end - temp_start);
      // update the state variables u,v
      temp_start = MPI_Wtime();
      step(du, dv, u, v);
      temp_end = MPI_Wtime();
      step_time += (temp_end - temp_start);
      if (k % m == 0) {
        writeInd = k / m;
        // calculate the norms
        temp_start = MPI_Wtime();
        nrmu = norm(u);
        nrmv = norm(v);
        temp_end = MPI_Wtime();
        norm_time += (temp_end - temp_start);
        stats[writeInd][0] = t;
        stats[writeInd][1] = nrmu;
        stats[writeInd][2] = nrmv;
      }
    }
    // write norms output
    if (rank == 0){
      char filename[50];
      sprintf(filename, "programData/test_problemsize-%d_%d-ranks_run-%d.dat", problem_size, size, run);
      FILE *fptr = fopen(filename, "w");
      fprintf(fptr, "#t\t\tnrmu\t\tnrmv\n");
      for (int k = 0; k < (T / m); k++) {
        fprintf(fptr, "%02.5f\t%02.5f\t%02.5f\n", stats[k][0], stats[k][1],
                stats[k][2]);
      }
      fclose(fptr);
    }

    end_time = MPI_Wtime();
    if (rank == 0){
      printf("%d,%d,%d,%f,%f,%f,%f,%f\n",problem_size, size, run, init_time, step_time, dxdt_time, norm_time,
            end_time - start_time);
      }
    }
  // Free allocated memory
  free(u);
  free(v);
  free(du);
  free(dv);
 
  MPI_Finalize();
  return 0;
}

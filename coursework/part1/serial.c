#include "params.h" // model & simulation parameters
#include <math.h>   // needed for tanh, used in init function
#include <stdio.h>  // needed for printing
#include <stdlib.h> // needed for malloc and free
#include <time.h>

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

  clock_t start,end;
  start=clock();
  double t = 0.0, nrmu, nrmv;
  int writeInd = 0;
  double stats[T / m][3];

  // Allocate memory for 1D arrays representing 2D grids
  double *u = (double *)malloc(M * N * sizeof(double));
  double *v = (double *)malloc(M * N * sizeof(double));
  double *du = (double *)malloc(M * N * sizeof(double));
  double *dv = (double *)malloc(M * N * sizeof(double));

  // Check for allocation success
  if (!u || !v || !du || !dv) {
    fprintf(stderr, "Error: Failed to allocate memory\n");
    return 1;
  }

  // initialize the state
  init(u, v);
  // time-loop
  for (int k = 0; k < T; k++) {
    // track the time
    t = dt * k;
    // evaluate the PDE
    dxdt(du, dv, u, v);
    // update the state variables u,v
    step(du, dv, u, v);
    if (k % m == 0) {
      writeInd = k / m;
      // calculate the norms
      nrmu = norm(u);
      nrmv = norm(v);
      stats[writeInd][0] = t;
      stats[writeInd][1] = nrmu;
      stats[writeInd][2] = nrmv;
    }
  }
  // write norms output
  FILE *fptr = fopen("serial.dat", "w");
  fprintf(fptr, "#t\t\tnrmu\t\tnrmv\n");
  for (int k = 0; k < (T / m); k++) {
    fprintf(fptr, "%02.5f\t%02.5f\t%02.5f\n", stats[k][0], stats[k][1],
            stats[k][2]);
  }
  fclose(fptr);

  // Free allocated memory
  free(u);
  free(v);
  free(du);
  free(dv);
  end = clock();
  double time_taken = ((double)(end - start))/CLOCKS_PER_SEC;
  printf("Successful completion of serial execution\n");
  printf("Time taken for serial execution: %f seconds\n", time_taken);

  return 0;
}

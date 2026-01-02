#define inline static inline
#include "params.h" // model & simulation parameters
#undef inline
#include <math.h>   // needed for tanh, used in init function
#include <stdio.h>  // needed for printing
#include <stdlib.h> // needed for malloc and free
#include <omp.h>


void init(double *u, double *v) {
  /* this has no data dependencies and is such trivially parallellised, 
testing reveals that collapsing both loops yields minor performance gains*/
  int idx,i,j;
  #pragma omp parallel for default(none) shared(u,v,N,M,uhi,ulo,vhi,vlo) private(idx,i,j) collapse(2)
  for (i = 0;i<M;i++){
    for (j = 0;j<N;j++){
      idx = i * N + j;
      u[idx] = ulo + (uhi - ulo) * 0.5 * (1.0 + tanh((i - M / 2) / 16.0)); //smooth gradient in vertical direction
      v[idx] = vlo + (vhi - vlo) * 0.5 * (1.0 + tanh((j - N / 2) / 16.0)); //smooth gradient in horizontal direction
    }
  }
}
void dxdt(double *du, double *dv, const double *u, const double *v) {
  /*Calculation of boundary points can be done independently of the interior points, 
  therefore by splitting the evaluation of the function into interior and boundary sections, 
  one can utilise a nowait clause to overlap calculation, 
  the implicit barrier at the end of interior points forces syncrhonisation before continuing*/
  double lapu, lapv, _u, _v;
  int up, down, left, right, idx,i,j; // indicies of neighboring grid positions
  #pragma omp parallel default(none) shared(u,v,du,dv,M,N,DD,R,d) private(idx,i,j,up,down,left,right,lapu,lapv,_u,_v)
  {
    #pragma omp for nowait 
    for (i = 1; i < M-1; i++) { //left and right edges (excluding corners)
      //start with left edge : j = 0
      up = (i + 1) * N;
      down = (i - 1) * N;
      left = i * N;
      right = i * N + 1;
      idx = i * N;
      _u = u[idx]; 
      _v = v[idx];
      lapu = u[up] + u[down] + u[left] + u[right] + -4.0 * _u; 
      lapv = v[up] + v[down] + v[left] + v[right] + -4.0 * _v;
      du[idx] = DD * lapu + f(_u, _v) + R * stim(i, 0); 
      dv[idx] = d * DD * lapv + g(_u, _v);
      //now right edge: j = N-1
      j = N-1;
      up = (i + 1) * N + (j);
      down = (i - 1) * N + (j);
      left = i * N + (N - 2);
      right = i * N + (j);
      idx = i * N + (j);
      _u = u[idx];
      _v = v[idx];
      lapu = u[up] + u[down] + u[left] + u[right] + -4.0 * _u; 
      lapv = v[up] + v[down] + v[left] + v[right] + -4.0 * _v;
      du[idx] = DD * lapu + f(_u, _v) + R * stim(i, j); 
      dv[idx] = d * DD * lapv + g(_u, _v);
   }
    #pragma omp for nowait
    for (j = 0; j < N; j++) { //bottom and top edges (incuding corners)
      //start with bottom edge : i = 0
      up = N + j;
      down = j;
      left = j == 0 ? j : j - 1;
      right = j == N - 1 ? j : j + 1;
      //here idx = j
      _u = u[j]; 
      _v = v[j];
      lapu = u[up] + u[down] + u[left] + u[right] + -4.0 * _u; 
      lapv = v[up] + v[down] + v[left] + v[right] + -4.0 * _v;
      du[j] = DD * lapu + f(_u, _v) + R * stim(0, j); 
      dv[j] = d * DD * lapv + g(_u, _v);
      //now top edge : i = M-1
      i = M-1;
      up = (i) * N + j;
      down = (M - 2) * N + j;
      left = (i) * N + (j == 0 ? j : j - 1);
      right = (i) * N + (j == N - 1 ? j : j + 1);
      idx = (i) * N + j;
      _u = u[idx]; 
      _v = v[idx];
      lapu = u[up] + u[down] + u[left] + u[right] + -4.0 * _u; 
      lapv = v[up] + v[down] + v[left] + v[right] + -4.0 * _v;
      du[idx] = DD * lapu + f(_u, _v) + R * stim(i, j); 
      dv[idx] = d * DD * lapv + g(_u, _v);
    }
    #pragma omp for  //interior points, testing revealed that collapsing loops here had negative performance impacts
    for (i = 1; i < M-1; i++) {
      for (j = 1; j < N-1; j++) {
        idx = i * N + j;
        down = (i - 1) * N + j;
        up = (i + 1) * N + j;
        left = i * N + j - 1;
        right = i * N + j + 1;
        _u = u[idx]; //create local copies
        _v = v[idx];
        lapu = u[up] + u[down] + u[left] + u[right] + -4.0 * _u; // discrete Laplacian
        lapv = v[up] + v[down] + v[left] + v[right] + -4.0 * _v;
        du[idx] = DD * lapu + f(_u, _v) + R * stim(i, j); // uses the differential equations to update du and dv
        dv[idx] = d * DD * lapv + g(_u, _v);
      }
    }
  }
  }


void step(const double *du, const double *dv, double *u, double *v) {
  /*Testing of several different parallelisation strategies all had negative performance impacts
  thus the final program is just executed serially, justification is available within the report*/
  int idx,i,j;
  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      idx = i * N + j;
     //for every grid point update u and v
      u[idx] += dt * du[idx];
      v[idx] += dt * dv[idx];
    }
  }
}

double norm(const double *x) {
  /*this is a classis example of a reduction.
   Testing revealed that collapsing the inner loop had worse performance compared to a trivial parallelisation */
  double nrmx = 0.0;
  int idx,i,j;
  #pragma omp parallel for default(none) shared(x,M,N) private(idx,i,j) reduction(+:nrmx)
  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      idx = i * N + j;
      nrmx += x[idx] * x[idx];
    }
  }
  return nrmx;
}


int main(int argc, char **argv) {
  double start,end;
  start = omp_get_wtime();
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
    if (k % m == 0) { // every m time steps, store norms
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
  FILE *fptr = fopen("part1.dat", "w");
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
  end = omp_get_wtime();

  //this is how i recorded the wall time for part 1
  //printf("%d,%.6f\n",omp_get_max_threads(), end - start);
  return 0;
}

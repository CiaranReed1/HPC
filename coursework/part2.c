#define inline static inline
#include "params.h" // model & simulation parameters
#undef inline
#include <math.h>   // needed for tanh, used in init function
#include <stdio.h>  // needed for printing
#include <stdlib.h> // needed for malloc and free
#include <mpi.h>  // MPI header


/*These have been defined within global scope as the function signatures cannot be changed*/
int sub_M,sub_N; //subdomain sizes for each rank
int i_first, i_last; //local indices for first and last rows in each rank
int global_i_first; //global index for first row in each rank
int rank, size;

void exchange_ghost_cells(double *u, double *v){
  //function to exchange ghost cells between ranks
  //using Sendrecv to avoid deadlocks and improve performance compared to separate blocking sends and recvs
  MPI_Status status;
  //sending with tag 0 and receiving with tag 1 for lower neighbor
  //sending with tag 1 and receiving with tag 0 for upper neighbor
  if (rank > 0){ //exchange from rank above  (rank = 0 represnts "highest" rows in global domain)
    MPI_Sendrecv(&u[i_first * sub_N],sub_N,MPI_DOUBLE,rank-1,0,&u[0],sub_N,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD,&status);
    MPI_Sendrecv(&v[i_first * sub_N],sub_N,MPI_DOUBLE,rank-1,0,&v[0],sub_N,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD,&status);
  }
  if (rank < size -1)//exchange from rank below
  {
    MPI_Sendrecv(&u[(i_last -1) * sub_N],sub_N,MPI_DOUBLE,rank+1,1,&u[(sub_M -1)*sub_N],sub_N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,&status);
    MPI_Sendrecv(&v[(i_last -1) * sub_N],sub_N,MPI_DOUBLE,rank+1,1,&v[(sub_M -1)*sub_N],sub_N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,&status);
  }
}

void init(double *u, double *v) {
  /*No dependencies within init, so each rank performs its own calculation */
  int idx;
  for (int i = i_first; i < i_last; i++) { //this iterates over local interior rows only, not ghost cells
    for (int j = 0; j < sub_N; j++) {
      idx = i * sub_N + j;
      //for calculatoin of u, the i needs to be scaled to the global index
      int i_scaled = global_i_first + (i - i_first);
      u[idx] = ulo + (uhi - ulo) * 0.5 * (1.0 + tanh((i_scaled - 0.5 * M) / 16.0));
      v[idx] = vlo + (vhi - vlo) * 0.5 * (1.0 + tanh((j - 0.5 * N) / 16.0));
    }
  }
   // top global boundary ghost cells 
  if (rank == 0) {
      for (int j = 0; j < sub_N; j++) {
          u[0 * sub_N + j] = u[i_first * sub_N + j];
          v[0 * sub_N + j] = v[i_first * sub_N + j];
      }
  }

  // bottom global boundary ghost cells
  if (rank == size - 1) {
      for (int j = 0; j < sub_N; j++) {
          u[(sub_M - 1) * sub_N + j] = u[(i_last - 1) * sub_N + j];
          v[(sub_M - 1) * sub_N + j] = v[(i_last - 1) * sub_N + j];
      }
  }
  // now do a halo exchange to fill interior ghost cells
  exchange_ghost_cells(u,v); 
}

void dxdt(double *du, double *dv, const double *u, const double *v) {
  double lapu, lapv, _u, _v;
  int up, down, left, right, idx;
  //this defines the global boundary conditions 
  for (int i = i_first; i < i_last; i++) {
    for (int j = 0; j < sub_N; j++) {
      idx = i * sub_N + j;
      int i_scaled = global_i_first + (i - i_first);
      if (i == i_first && rank == 0) { //global boundary condition at top edge only on top rank
        up = idx;
      } else {
        up = (i - 1) * sub_N + j;
      }
      if (i == i_last -1 && rank == size -1) { //global boundary condition at bottom edge only on bottom rank
        down = idx;
      } else {
        down = (i + 1) * sub_N + j;
      }
      if (j == 0) { //global boundary condition at left edge (all ranks)
        left = idx;
      } else {
        left = idx - 1;
      }
      if (j == sub_N - 1) { //global boundary condition at right edge (all ranks)
        right = idx;
      } else {
        right = idx + 1;
      }
      _u = u[idx];
      _v = v[idx];
      lapu = u[up] + u[down] + u[left] + u[right] + -4.0 * _u;
      lapv = v[up] + v[down] + v[left] + v[right] + -4.0 * _v;
      du[idx] = DD * lapu + f(_u, _v) + R * stim(i_scaled, j);
      dv[idx] = d * DD * lapv + g(_u, _v);
    }
  }
  //ghost cells for du and dv are not needed, so no exchange is required here
}

void step(const double *du, const double *dv, double *u, double *v) {
  int idx;
  for (int i = i_first; i < i_last; i++) {
    for (int j = 0; j < sub_N; j++) {
      idx = i * sub_N + j;
      u[idx] += dt * du[idx];
      v[idx] += dt * dv[idx];
    }
  }
  // now do a halo exchange to update ghost cells
  exchange_ghost_cells(u,v);
}

double norm(const double *x) {
  /*Initially I had solved this using an MPI_Allreduce, however this causes floating point errors and the result differed from the serial result by a tiny margin.
  Given that the assignment brief tells us to ensure the resulting data exactly matches the serial data, I have attempted to fix this.
  My solution gathers all the local arrays onto one rank and calculates the norm serially to ensure coherant order of operations with a serial execution*/

  int local_count = (i_last - i_first) * sub_N; //number of elements in local array excluding ghost cells
  int *counts = NULL, *displs = NULL;;
  double *global_x = NULL;
  if (rank ==0){
      counts = malloc(size * sizeof(int));
      displs = malloc(size * sizeof(int));
  }
  MPI_Gather(&local_count,1,MPI_INT,counts,1,MPI_INT,0,MPI_COMM_WORLD);  //gathers the counts of elements from each rank
  //now gather all local counts to rank 0, work out the displacements for each rank
  if (rank ==0){
      int total_count =0;
      displs[0]=0;
      for (int r =0; r < size; r++){
          total_count += counts[r];
          if (r > 0){
            //works out displacements, to know where to place the data from each rank
            displs[r] = displs[r-1] + counts[r-1]; 
          }
      }
      global_x = malloc(total_count * sizeof(double));
  }
  //gathers all local arrays to rank 0 to "global_x", using counts and displs to know how much data to expect from each rank and where to place it
  MPI_Gatherv(&x[i_first * sub_N], local_count, MPI_DOUBLE,
          global_x, counts, displs, MPI_DOUBLE,
          0, MPI_COMM_WORLD); 
  /*Calculates the global norm from the gathered global array on rank 0 */
  double global_nrmx = 0.0;
  if (rank ==0){
    free(counts);
    free(displs);
    global_nrmx =0.0;
    for (int i =0; i < (M * N); i++){
        global_nrmx += global_x[i] * global_x[i];
    }
    free(global_x);
  }
  /*I could broadcast the global norm from rank 0 to all other ranks if needed, but this is unnecessary here and returning 0.0 for other ranks is fine*/
  return global_nrmx; 
}

int main(int argc, char **argv) {

  /*
  This is how the problem size was scaled to perform weak scaling tests, 
  note each dmimenson was scaled by sqrt(problem_size) to keep the total number of grid points proportional to problem_size 

  int problem_size = argc > 1 ? atoi(argv[1]) : 1;
  double scale = sqrt(problem_size);
  M = (int)(M * scale + 0.5);  // round to nearest int
  N = (int)(N * scale + 0.5);
  M += M % 2;  // add 1 if odd
  N += N % 2;  // add 1 if odd
  */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double t = 0.0, nrmu, nrmv;
  int writeInd = 0;
  double stats[T / m][3];

  //Split into a 4x1 Grid add remainder to last rank (splitting along slower axis leaving faster computation for each rank)
  int local_rows = M / size;
  if (rank == size - 1) {
      local_rows += M % size;
  } 

  //set subdomain sizes
  sub_M = local_rows + 2; 
  sub_N = N; 

  //set local and global indices
  i_first = 1; 
  i_last = sub_M -1; 
  global_i_first = rank * (M / size);


  // Allocate memory for 1D arrays representing 2D grids using the subdomain sizes
  double *u = (double *)malloc(sub_M * sub_N * sizeof(double));
  double *v = (double *)malloc(sub_M * sub_N * sizeof(double));
  double *du = (double *)malloc(sub_M * sub_N * sizeof(double));
  double *dv = (double *)malloc(sub_M * sub_N * sizeof(double));

  double init_time = 0.0, step_time = 0.0, norm_time = 0.0, dxdt_time = 0.0;
  double start_time, end_time, temp_start, temp_end;

  // Check for allocation success
  if (!u || !v || !du || !dv) {
    fprintf(stderr, "Error: Failed to allocate memory\n");
    return 1;
  }
  // initialize the state

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
      //note only rank 0 will have the correct norm values, other ranks will have 0.0
      stats[writeInd][0] = t;
      stats[writeInd][1] = nrmu;
      stats[writeInd][2] = nrmv;
    }
  }
  // write norms output
  if (rank == 0){
    FILE *fptr = fopen("part2.dat", "w");
    fprintf(fptr, "#t\t\tnrmu\t\tnrmv\n");
    for (int k = 0; k < (T / m); k++) {
      fprintf(fptr, "%02.5f\t%02.5f\t%02.5f\n", stats[k][0], stats[k][1],
              stats[k][2]);
    }
    fclose(fptr);
  }

  end_time = MPI_Wtime();
  /*
  the following was how i measured the wall time for weak scaing tests

  if (rank == 0){
    printf("%d,%d,%f,%f,%f,%f,%f\n",problem_size, size, init_time, step_time, dxdt_time, norm_time,
          end_time - start_time);
    }
  */

  // Free allocated memory
  free(u);
  free(v);
  free(du);
  free(dv);
 
  MPI_Finalize();
  return 0;
}

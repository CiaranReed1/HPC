#include "params.h" // model & simulation parameters
#include <math.h>   // needed for tanh, used in init function
#include <stdio.h>  // needed for printing
#include <stdlib.h> // needed for malloc and free
#include <mpi.h>  // MPI header

int M = 1024;  // x domain size
int N = 512;   // y domain size
int sub_M,sub_N; //subdomain sizes for each rank
int i_first, i_last; //local indices for first and last rows in each rank
int global_i_first; //global index for first row in each rank
int rank, size;

void exchange_ghost_cells(double *u, double *v){
  //function to exchange ghost cells between ranks
  //using Sendrecv to avoid deadlocks and improve performance compared to separate sends and recvs
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
  for (int i = i_first; i < i_last; i++) { //this iterated over local rows only, not ghost cells
    for (int j = 0; j < sub_N; j++) {
      idx = i * sub_N + j;
      //for calculatoin of u, the i needs to be scaled to the global index
      int i_scaled = global_i_first + (i - i_first);
      
      if (i == i_first && j < 5) {
          printf("Rank %d: i=%d, j=%d, i_scaled=%d, u=%f, v=%f\n",
                rank, i, j, i_scaled, u[idx], v[idx]);
          fflush(stdout);
      }
     
      u[idx] = ulo + (uhi - ulo) * 0.5 * (1.0 + tanh((i_scaled - 0.5 * M) / 16.0));
      v[idx] = vlo + (vhi - vlo) * 0.5 * (1.0 + tanh((j - 0.5 * N) / 16.0));
    }
  }
  // now do a halo exchange to fill ghost cells
  //exchange_ghost_cells(u,v);
}

void dxdt(double *du, double *dv, const double *u, const double *v) {
  double lapu, lapv, _u, _v;
  int up, down, left, right, idx;
  //this defines the global boundary conditions 
  for (int i = i_first; i < i_last; i++) {
    for (int j = 0; j < sub_N; j++) {
      idx = i * sub_N + j;
      int i_scaled = global_i_first + (i - i_first);
      if (i == i_first && rank == 0) { //global boundary condition at top edge
        up = idx;
      } else {
        up = (i - 1) * sub_N + j;
      }
      if (i == i_last -1 && rank == size -1) {
        down = idx;
      } else {
        down = (i + 1) * sub_N + j;
      }
      if (j == 0) {
        left = idx;
      } else {
        left = idx - 1;
      }
      if (j == sub_N - 1) {
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
  double nrmx = 0.0;
  int idx;
  for (int i = i_first; i < i_last; i++) {
    for (int j = 0; j < sub_N; j++) {
      idx = i * sub_N + j;
      nrmx += x[idx] * x[idx];
    }
  }
  MPI_Allreduce(&nrmx,&nrmx,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return nrmx;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    M = 8;
    N = 8;
    int local_rows = M / size;
    if (rank == size -1) local_rows += M % size;

    sub_M = local_rows + 2; // +2 for ghost cells
    sub_N = N;

    printf("Rank %d of %d: local_rows=%d, sub_M=%d, sub_N=%d\n", rank, size, local_rows, sub_M, sub_N);
    i_first = 1;
    i_last = sub_M -1;
    if (rank < size-1)
        global_i_first = rank * (M/size);
    else
        global_i_first = M - local_rows; // last rank gets remainder


    double *u = (double*)malloc(sub_M * sub_N * sizeof(double));
    double *v = (double*)malloc(sub_M * sub_N * sizeof(double));
    double *du = (double*)malloc(sub_M * sub_N * sizeof(double));
    double *dv = (double*)malloc(sub_M * sub_N * sizeof(double));

    if (!u || !v || !du || !dv) {
        fprintf(stderr, "Memory allocation failed\n");
        MPI_Finalize();
        return 1;
    }

    // Initialize arrays
    init(u, v);

    // Print initial values (first rank only)
    if (rank == 0) {
        printf("Initial u:\n");
        for (int i = i_first; i < i_last; i++) {
            for (int j = 0; j < sub_N; j++)
                printf("%f ", u[i*sub_N + j]);
            printf("\n");
        }
        printf("Initial v:\n");
        for (int i = i_first; i < i_last; i++) {
            for (int j = 0; j < sub_N; j++)
                printf("%f ", v[i*sub_N + j]);
            printf("\n");
        }
    }

    

    free(u); free(v); free(du); free(dv);
    MPI_Finalize();
    return 0;
}



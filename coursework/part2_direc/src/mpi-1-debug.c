#include "params.h" // model & simulation parameters
#include <math.h>   // needed for tanh, used in init function
#include <stdio.h>  // needed for printing
#include <stdlib.h> // needed for malloc and free
#include <mpi.h>  // MPI header

int M = 9;  // x domain size
int N = 8;   // y domain size
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

int main(int argc, char **argv)
{
    
  int problem_size = argc > 1 ? atoi(argv[1]) : 1;
  int nruns = argc > 2 ? atoi(argv[2]) : 1;

  double scale = sqrt(problem_size);
  M = (int)(M * scale + 0.5);  // round to nearest int
  N = (int)(N * scale + 0.5);
  M += M % 2;  // add 1 if odd
  N += N % 2;  // add 1 if odd

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double t = 0.0, nrmu, nrmv;
  int writeInd = 0;
  double stats[T / m][3];

  //Split into a 4x1 Grid (splitting along slower axis leaving faster computation for each rank)
  int local_rows = M / size;
  if (rank == size - 1) {
      local_rows += M % size;
  } // add remainder to last rank


  sub_M = local_rows + 2; // +2 for ghost cells
  sub_N = N; // full N for each rank
  
  i_first = 1; //0 is ghost cell
  i_last = sub_M -1; //last is ghost cell 
  
  global_i_first = rank * (M / size);

    /* --- allocate arrays --- */
    double *u  = malloc(sub_M * sub_N * sizeof(double));
    double *v  = malloc(sub_M * sub_N * sizeof(double));

    if (!u || !v) {
        fprintf(stderr, "Rank %d: allocation failed\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* --- initialize --- */
    init(u, v);

    /* --- ordered printing --- */
    MPI_Barrier(MPI_COMM_WORLD);

    for (int r = 0; r < size; r++) {
        if (rank == r) {
            printf("\n==============================\n");
            printf("Rank %d\n", rank);
            printf("global_i_first = %d\n", global_i_first);
            printf("local rows (physical) = %d\n", i_last - i_first);
            printf("sub_M = %d, sub_N = %d\n", sub_M, sub_N);

            for (int i = 0; i < sub_M; i++) {
                int gi = (i >= i_first && i < i_last)
                           ? global_i_first + (i - i_first)
                           : -1;  // ghost row

                printf("i = %2d (global %2d): ", i, gi);
                for (int j = 0; j < sub_N; j++) {
                    int idx = i * sub_N + j;
                    printf("%+.4e ", u[idx]);
                }
                printf(" | ");
                for (int j = 0; j < sub_N; j++) {
                    int idx = i * sub_N + j;
                    printf("%+.4e ", v[idx]);
                }
                printf("\n");
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    free(u);
    free(v);

    MPI_Finalize();
    return 0;
}

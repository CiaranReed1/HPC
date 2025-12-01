#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void gather_nonblocking(const int *send, int *recvbuf, MPI_Comm comm)
{
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  if (rank == 0) {
   //??
  } else {
   //??
  }
}

void gather_blocking(const int *send, int *recvbuf, MPI_Comm comm)
{
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  if (rank == 0) {
  //??
  } else {
  //??
  }
}
/*
	* Process 0: A
	* Process 1: B
	* Process 2: C
	* Process 3: D
	*
	* gather_...(...) ->
	* Process 0: A, B, C, D
*/
int main(int argc, char **argv){
	MPI_Init(&argc, &argv);

	MPI_Comm comm = MPI_COMM_WORLD;
	int rank, size;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	int local;
	local = rank * rank;

	int *blocking = NULL;
	int *nonblocking = NULL;
	if (rank == 0) {
		/* Allocate space for output arrays -- 1 int per process. */
		blocking = malloc(size*sizeof(*blocking));
		nonblocking = malloc(size*sizeof(*nonblocking));
	}
	double start, end;
	start = MPI_Wtime();
	for (int i = 0; i < 100; i++) {
		gather_blocking(&local, blocking, comm);
	}
	end = MPI_Wtime();
	if (rank == 0) {
		printf("Blocking gather takes %.3g s\n", (end - start)/100);
	}
	start = MPI_Wtime();
	for (int i = 0; i < 100; i++) {
		gather_nonblocking(&local, nonblocking, comm);
	}
	end = MPI_Wtime();
	if (rank == 0) {
		printf("Non-blocking gather takes %.3g s\n", (end - start)/100);
	}
	free(blocking);
	free(nonblocking);
	MPI_Finalize();
	return 0;
}

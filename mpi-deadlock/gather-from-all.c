#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void gather_nonblocking(const int *send, int *recvbuf, MPI_Comm comm,int nEntries)
{
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  

	if (rank == 0) {
		MPI_Request requests[size-1];
		for (int i =1;i < size;i++){
			MPI_Irecv(recvbuf+i*nEntries,nEntries,MPI_INT,MPI_ANY_SOURCE,0,comm,&requests[i-1]);
		}
		MPI_Waitall(size-1, requests, MPI_STATUS_IGNORE);
	} else {
		MPI_Send(send,nEntries,MPI_INT,0,0,comm);
	}
}

void gather_blocking(const int *send, int *recvbuf, MPI_Comm comm,int nEntries)
{
	int size, rank;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	if (rank == 0) {
		for (int i =1;i<size;i++){
			MPI_Recv(recvbuf+i*nEntries,nEntries,MPI_INT,MPI_ANY_SOURCE,0,comm,MPI_STATUS_IGNORE);
		}

	} else {
		MPI_Send(send,nEntries,MPI_INT,0,0,comm);
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
	//the messages being passed are the square of the rank
*/
int main(int argc, char **argv){
	MPI_Init(&argc, &argv);

	MPI_Comm comm = MPI_COMM_WORLD;
	int rank, size;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	int *blocking = NULL;
	int *nonblocking = NULL;

	double blocking_time = 0.0;
	double nonblocking_time = 0.0;

	int nMessages = argc > 1 ? atoi(argv[1]) : 100;
	int nEntries = argc > 2 ? atoi(argv[2]) : 1;

	if (rank == 0) {
		/* Allocate space for output arrays -- 1 int per process. */
		blocking = malloc(size*sizeof(*blocking)*nEntries);
		nonblocking = malloc(size*sizeof(*nonblocking)*nEntries);
	}

	int *local = NULL;
	local = malloc(nEntries * sizeof(int));
	for (int i = 0; i < nEntries; i++) {
		local[i] = rank * rank;
	}

	double start, end;
	
	start = MPI_Wtime();
	for (int i = 0; i < nMessages; i++) {
		gather_blocking(local, blocking, comm,nEntries);
	}
	end = MPI_Wtime();
	blocking_time += (end - start);

	if (rank == 0) {
		printf("With %d ranks, Average blocking gather over %d messages with %d entries takes %.3g s\n", size, nMessages, nEntries, blocking_time/(nMessages));
		}
	
	
	start = MPI_Wtime();
	for (int i = 0; i < nMessages; i++) {
		gather_nonblocking(local, nonblocking, comm,nEntries);
	}
	end = MPI_Wtime();
	nonblocking_time += (end - start);
	
	if (rank == 0) {
		printf("With %d ranks, Average non-blocking gather over %d messages with %d entries takes %.3g s\n", size, nMessages, nEntries, nonblocking_time/(nMessages));
		}
	free(blocking);
	free(nonblocking);
	free(local);
	MPI_Finalize();
	return 0;
}

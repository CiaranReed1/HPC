#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

static void ping_pong(void *buffer, int count, MPI_Datatype dtype, MPI_Comm comm, int rank){
	/* Implement a ping pong.
	*
	* rank 0 should send count bytes from buffer to rank 1
	* rank 1 should then send the received data back to rank 0
	*
	*/
	MPI_Status status;
	if (rank == 0){
		MPI_Send(buffer,count,dtype,1,0,comm);
		MPI_Recv(buffer,count,dtype,1,0,comm,&status);
	}
	if (rank == 1){
		MPI_Recv(buffer,count,dtype,0,0,comm,&status);
		MPI_Send(buffer,count,dtype,0,0,comm);
	}
	else{
		/* nothing to do */
	}
}

int main(int argc, char **argv){
	MPI_Init(&argc, &argv);

	int nbytes, rank;
	double start,end;
	char *buffer;
	MPI_Comm comm;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	comm = MPI_COMM_WORLD;
	nbytes = argc > 1 ? atoi(argv[1]) : 1;

	buffer = calloc(nbytes, sizeof(*buffer));

	start = MPI_Wtime();
	ping_pong(buffer, nbytes, MPI_CHAR, comm,rank);
	end = MPI_Wtime();
	printf("Time for ping pong of %d bytes on rank %d: %f seconds\n", nbytes, rank, end - start);
	free(buffer);

	MPI_Finalize();
	printf("Ping pong completed\n");
	return 0;
}

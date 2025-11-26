#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

static void ping_pong(void *buffer, int count, MPI_Datatype dtype, MPI_Comm comm){
	/* Implement a ping pong.
	*
	* rank 0 should send count bytes from buffer to rank 1
	* rank 1 should then send the received data back to rank 0
	*
	*/
}

int main(int argc, char **argv){
	MPI_Init(&argc, &argv);

	int nbytes, rank;
	char *buffer;
	double start, end;
	MPI_Comm comm;

	comm = MPI_COMM_WORLD;
	nbytes = argc > 1 ? atoi(argv[1]) : 1;

	buffer = calloc(nbytes, sizeof(*buffer));

	ping_pong(buffer, nbytes, MPI_CHAR, comm);

	free(buffer);

	MPI_Finalize();
	return 0;
}

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
	int iterations = 10000;
	int byte_size[14] = {1,4,16,64,256,1024,4096,16384,65536,262144,1048576,4194304,16777216,67108864}; 
	double start,end;
	char *buffer;
	MPI_Comm comm;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	comm = MPI_COMM_WORLD;

	for (int j = 0; j < 14;j++){
		nbytes = byte_size[j];
		buffer = calloc(nbytes, sizeof(*buffer));
		start = MPI_Wtime();
		for (int i = 0;i < iterations;i++){
			ping_pong(buffer, nbytes, MPI_CHAR, comm,rank);
		}
		end = MPI_Wtime();
		free(buffer);
		if (rank == 0){
			printf("%d,%d,%f\n",iterations,nbytes,end-start);
		}
			
	}

	MPI_Finalize();
	return 0;
}

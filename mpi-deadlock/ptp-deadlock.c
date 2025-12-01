#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv)
{
	int rank, size;

	int *sendbuf;
	int *recvbuf;
	int buffsize;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int nentries = argc > 1 ? atoi(argv[1]) : 1;

	buffsize = nentries * sizeof(*sendbuf);
	sendbuf = calloc(nentries, buffsize);  // Use calloc to initialize to zero, allocating the size of 1 integer for each entry
	recvbuf = malloc(nentries * buffsize); 

	MPI_Buffer_attach(sendbuf, buffsize);

	for (int i = 0; i < nentries; i++) {
		sendbuf[i] = rank;
	}


	if (rank == 0) {
		printf("Sending %d entries per rank\n", nentries);
		MPI_Bsend(sendbuf, nentries, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(recvbuf, nentries, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	} else if (rank == 1) {
		MPI_Bsend(sendbuf, nentries, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Recv(recvbuf, nentries, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	MPI_Buffer_detach(&sendbuf, &buffsize);
	printf("Rank : [%d] Success! First entry of sendbuf is %d; first of recvbuf is %d\n", rank, sendbuf[0], recvbuf[0]);
	MPI_Finalize();
	return 0;
}
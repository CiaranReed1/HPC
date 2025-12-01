#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv)
{
	int rank, size;

	int *sendbuf;
	int *sendarr;
	int *recvarr;
	int buffsize;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int nentries = argc > 1 ? atoi(argv[1]) : 1;

	buffsize = nentries * sizeof(*sendarr) + MPI_BSEND_OVERHEAD;

	sendbuf = malloc(buffsize);  //initialize buffer

	recvarr = malloc(nentries * sizeof(*recvarr)); //intialise send and recieve arrarrys
	sendarr = calloc(nentries, sizeof(*sendarr));

	MPI_Buffer_attach(sendbuf, buffsize);

	for (int i = 0; i < nentries; i++) {
		sendarr[i] = rank;
	}


	if (rank == 0) {
		printf("Sending %d entries per rank\n", nentries);
		MPI_Bsend(sendarr, nentries, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(recvarr, nentries, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	} else if (rank == 1) {
		MPI_Bsend(sendarr, nentries, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Recv(recvarr, nentries, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	MPI_Buffer_detach(&sendbuf, &buffsize);
	printf("Rank : [%d] Success! First entry of sendbuf is %d; first of recvbuf is %d\n", rank, sendbuf[0], recvarr[0]);
	MPI_Finalize();
	return 0;
}
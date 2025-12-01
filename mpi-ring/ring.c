#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv){
	MPI_Init(&argc, &argv);

	MPI_Comm comm;
	comm = MPI_COMM_WORLD;

	int rank, size,left,right;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	left = (rank - 1 + size) % size;
	right = (rank + 1) % size;

	int local_value = rank;
	int summed_value = 0;	// initialized sum: 0

	printf("[%d] Before reduction: local value is %d; summed value is %d\n",
	      rank, local_value, summed_value);
	double start = MPI_Wtime();
	summed_value = local_value;
    for(int i=0; i <size -1;i++){
		MPI_Request request;
		MPI_Isend(&local_value,1,MPI_INT,right,0,comm,&request);
		MPI_Recv(&local_value,1,MPI_INT,left,0,comm,MPI_STATUS_IGNORE);
		summed_value+=local_value;
		MPI_Wait(&request,MPI_STATUS_IGNORE);
	}
	double end = MPI_Wtime();
	printf("[%d] After reduction: local value is %d; summed value is %d (should be %d)\n",
	      rank, local_value, summed_value, size*(size-1)/2);
	if (rank == 0){
		printf("\tReduction on %d processes took %e seconds.\n", size, end-start);
	}
	MPI_Finalize();
	return 0;
}

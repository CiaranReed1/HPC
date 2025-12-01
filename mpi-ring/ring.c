#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv){
	MPI_Init(&argc, &argv);

	MPI_Comm comm;
	comm = MPI_COMM_WORLD;

	int rank, size,left,right;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	left = (rank - 1 + size) % size;
	right = (rank + 1) % size;

	int nIterations = argc > 1 ? atoi(argv[1]) : 1000;

	int local_value,received_local;
	int summed_value;	
	int desired = (size * (size - 1)) / 2;

	double start, end, mean_time = 0.0,score = 0.0;

	for (int iter = 0; iter < nIterations; iter++){
		summed_value = 0;
		local_value = rank;
		received_local = -1;
		start = MPI_Wtime();
		summed_value = local_value;
		for(int i=0; i <size -1;i++){
			MPI_Request request;
			MPI_Isend(&local_value,1,MPI_INT,right,0,comm,&request);
			MPI_Recv(&received_local,1,MPI_INT,left,0,comm,MPI_STATUS_IGNORE);
			summed_value+=received_local;
			MPI_Wait(&request,MPI_STATUS_IGNORE);
			local_value = received_local;
		}
		end = MPI_Wtime();
		mean_time += (end - start);
		if (summed_value == desired){
			score += 1.0;
		}
	};
	score /= nIterations;
	mean_time /= nIterations;
	if (rank == 0){
		printf("%d,%d,%e,%f \n", nIterations,size, mean_time,score);
	}
	MPI_Finalize();
	return 0;
}

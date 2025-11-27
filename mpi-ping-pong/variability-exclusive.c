#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static void ping_pong(void *buffer, int count, MPI_Datatype dtype, MPI_Comm comm, int rank,int size){
	/* Implement a ping pong.
	*
	* rank 0 should send count bytes from buffer to rank 1
	* rank 1 should then send the received data back to rank 0
	*
	*/
	MPI_Status status;
	if (rank < (size/2)){
		MPI_Send(buffer,count,dtype,rank + (size/2),0,comm);
		MPI_Recv(buffer,count,dtype,rank + (size/2),0,comm,&status);
	}
	else{
		MPI_Recv(buffer,count,dtype,rank - (size/2),0,comm,&status);
		MPI_Send(buffer,count,dtype,rank - (size/2),0,comm);
	}

}

int main(int argc, char **argv){
	MPI_Init(&argc, &argv);

	int nbytes, rank, size; 
	int iterations = 10000;
	int byte_size[14] = {1,4,16,64,256,1024,4096,16384,65536,262144,1048576,4194304,16777216,67108864}; 
	double start,end,mean,stdev,max,min,total;
	char *buffer;
	double *times;
	MPI_Comm comm;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (rank == 0){
		printf("iterations,nbytes,mean,stdev,min,max,total_time\n");
	}
	comm = MPI_COMM_WORLD;

	for (int j = 0; j < 14;j++){
		nbytes = byte_size[j];
		buffer = calloc(nbytes, sizeof(*buffer));
		times = calloc(iterations, sizeof(*times));
		min = 0.0;
	    max = 0.0;
		mean = 0.0;
		stdev = 0.0;
		total = 0.0;
		for (int i = 0;i < iterations;i++){
			start = MPI_Wtime();
			ping_pong(buffer, nbytes, MPI_CHAR, comm,rank,size);
			end = MPI_Wtime();
			times[i] = end - start;
			total += times[i];
			if (i == 0){
				min = times[i];
				max = times[i];
			}
			else{
				times[i] < min ? (min = times[i]) : 0;
				times[i] > max ? (max = times[i]) : 0;
			}
		}
		mean = total / iterations;
		for (int i = 0;i < iterations;i++){
			stdev += (times[i] - mean) * (times[i] - mean);
		}
		stdev = stdev / iterations;
		stdev = sqrt(stdev);
		free(buffer);
		free(times);
		if (rank == 0){
			printf("%d,%d,%f,%f,%f,%f,%f\n",iterations,nbytes,mean,stdev,min,max,total);
		}
			
	}

	MPI_Finalize();
	return 0;
}

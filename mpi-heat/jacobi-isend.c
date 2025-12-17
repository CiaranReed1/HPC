/*
Adapted from Argonne National Lab MPI tutorial
https://www.mcs.anl.gov/research/projects/mpi/tutorial/mpiexmpl/contents.html
*/
#include <stdio.h>
#include <math.h>
#include "mpi.h"

/* This example handles a maxn x maxn mesh, on 4 processors only. */
#define maxn 24

int main( int argc, char **argv ){

    int        rank, value, size, errcnt, toterr, i, j, iter;
    int        i_first, i_last;
    double     diffnorm, gdiffnorm;
    double     xloc[(maxn/4)+2][maxn];	// Note: this is technically oversized for Rank 0 & 3
    double     xnew[(maxn/4)+2][maxn];
	MPI_Status status;
	MPI_Request requests[4] =  {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
	
    MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    if (size != 4){	// Hardcoding a four-process decomposition
    	MPI_Abort( MPI_COMM_WORLD, 1 );
	}
	
    /* xloc[][0] is lower ghostpoints, xloc[][maxn+2] is upper */

	/*
    
    () : Physical Boundary Point
    {} : Ghost Boundary Point
    [] : Interior Point
    == : Rank boundary
    
    {}{}{}{}{}{}{}{}{}{}{}{}
    ()()()()()()()()()()()()		.--> second index
    ()[][][][][][][][][][]()		|
    ()[][][][][][][][][][]()		v
    {}{}{}{}{}{}{}{}{}{}{}{}		first index
   ==========================		
    {}{}{}{}{}{}{}{}{}{}{}{}
    ()[][][][][][][][][][]()
    ()[][][][][][][][][][]()
    ()[][][][][][][][][][]()
    {}{}{}{}{}{}{}{}{}{}{}{}
   ==========================
    {}{}{}{}{}{}{}{}{}{}{}{}
    ()[][][][][][][][][][]()
    ()[][][][][][][][][][]()
    ()[][][][][][][][][][]()
    {}{}{}{}{}{}{}{}{}{}{}{}
   ==========================
    {}{}{}{}{}{}{}{}{}{}{}{}
    ()[][][][][][][][][][]()
    ()[][][][][][][][][][]()
    ()()()()()()()()()()()()
    {}{}{}{}{}{}{}{}{}{}{}{}
    
    Top:				 	T=100
    Left:					T=20
    Right:					T=80;
    Bottom:					T=0
    Interior init:			T=rank
    
    */
	double start_time = MPI_Wtime();

    /* Note that top and bottom processes have one less row of interior points */
    i_first = 1;
    i_last  = maxn/size;
    if (rank == 0){
    	i_first++;
    }
    if (rank == size - 1){
    	i_last--;
    }

	/* Fill the interior with data */ 
    for (i=i_first; i<=i_last; i++){ 
		for (j=0; j<maxn; j++){
		    xloc[i][j] = (double)rank;
		}
	}
	/* Fill Boundaries */
	if (rank==0){
		for (j=0; j<maxn; j++){
			xloc[i_first-1][j] 	= 100.0;	// top
		}
	}
	if (rank==3){
		for (j=0; j<maxn; j++){
			xloc[i_last+1][j] 	= 0.0;		// bottom
		}
	}
	for (i=i_first; i<=i_last; i++){
		xloc[i][0] 				= 20.0;		// left
		xloc[i][maxn-1] 		= 80.0;		// right
	}
	
	/* All the Ghost points / places being recieved into are uninitialized */


	/* Print Initial Data */
	for (int k=0; k<size; k++){
		MPI_Barrier(MPI_COMM_WORLD);
		if (k==rank){
			printf("Rank: %d\n",rank);
			for (i=i_first-1; i<=i_last+1; i++){
				for (j=0; j<maxn; j++){
					printf("%2.1f ",xloc[i][j]);
				}
				printf("\n");
			}
		}
	}
	

    iter = 0; 
    do {
		/* Send up unless I'm at the top, then receive from below */
		/* Note the use of xloc[i] for &xloc[i][0] */
		if (rank > 0){
			printf("Recieving row %d from rank %d to rank %d\n", 0, rank-1, rank);
		    MPI_Irecv( xloc[0], maxn, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, requests );
		}
		if (rank < size - 1){
			printf("Sending row %d from rank %d to rank %d\n", maxn/size, rank, rank+1);
		    MPI_Isend( xloc[maxn/size], maxn, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, requests + 1);
		}
		
		/* this does not deadlock as the top one will receive whilst all the others send, and it will cascade down the sim*/

		/* Send down unless I'm at the bottom */
		if (rank < size - 1){
			printf("Recieving row %d from rank %d to rank %d\n", maxn/size+1, rank+1, rank);
		    MPI_Irecv( xloc[maxn/size+1], maxn, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, requests + 2 );
		}
		if (rank > 0){ 
			printf("Sending row %d from rank %d to rank %d\n", 1, rank, rank-1);
		    MPI_Isend( xloc[1], maxn, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, requests + 3 );
		}
		MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
		
		/* Compute new values (but not on boundary) */
		iter++;
		diffnorm = 0.0;
		for (i=i_first; i<=i_last; i++) {
		    for (j=1; j<maxn-1; j++) {
				xnew[i][j] = (xloc[i][j+1] + xloc[i][j-1] + xloc[i+1][j] + xloc[i-1][j]) / 4.0;
				diffnorm += (xnew[i][j] - xloc[i][j]) *  (xnew[i][j] - xloc[i][j]);
		    }
		}
		/* Only transfer the interior points */
		for (i=i_first; i<=i_last; i++) {
		    for (j=1; j<maxn-1; j++){
				xloc[i][j] = xnew[i][j];
			}
		}
		/* on each rank, diffnorm is the sum of squared differences for that rank
		but then we reduce and find the maximum of all the ranks and save that on each*/
		MPI_Allreduce( &diffnorm, &gdiffnorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
		gdiffnorm = sqrt( gdiffnorm );
		if (rank == 0) {
			printf( "At iteration %d, diff is %e\n", iter, gdiffnorm );
		}
    } while ((gdiffnorm > 1e-3) && (iter < 1000));

	double end_time = MPI_Wtime();
	if (rank==0){
		printf("\n");
		printf("Successful completion\n");
		printf("Converged after %d iterations\n", iter);
		printf("Total time: %f seconds\n", end_time - start_time);
	}
    MPI_Finalize( );

	
    return 0;
}

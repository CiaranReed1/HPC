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
    int        j_first, j_last;
    double     diffnorm, gdiffnorm;
    double     xloc[maxn][(maxn/4)+2];
    double     xnew[maxn][(maxn/4)+2];
    double     xedge[maxn+2];
    double     tol = 1e-3;
	MPI_Status status;
	
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
    || : Rank boundary
    
               ||            ||            || 					.--> second index
    {}()()(){} || {}()()(){} || {}()()(){} || {}()()(){}		|
    {}()[][]{} || {}[][][]{} || {}[][][]{} || {}[][](){}		v
    {}()[][]{} || {}[][][]{} || {}[][][]{} || {}[][](){}		first index
    {}()[][]{} || {}[][][]{} || {}[][][]{} || {}[][](){}
    {}()[][]{} || {}[][][]{} || {}[][][]{} || {}[][](){}
    {}()[][]{} || {}[][][]{} || {}[][][]{} || {}[][](){}
    {}()[][]{} || {}[][][]{} || {}[][][]{} || {}[][](){}
    {}()[][]{} || {}[][][]{} || {}[][][]{} || {}[][](){}
    {}()[][]{} || {}[][][]{} || {}[][][]{} || {}[][](){}
    {}()[][]{} || {}[][][]{} || {}[][][]{} || {}[][](){}
    {}()[][]{} || {}[][][]{} || {}[][][]{} || {}[][](){}
    {}()()(){} || {}()()(){} || {}()()(){} || {}()()(){}
               ||            ||            || 
    
    Top:				 	T=100
    Left:					T=20
    Right:					T=80;
    Bottom:					T=0
    Interior init:			T=rank
    
    */
	double start, end;
	start = MPI_Wtime();
    /* Note that top and bottom processes have one less column of interior points */
    j_first = 1;
    j_last  = maxn/size;
    if (rank == 0){
    	j_first++;
    }
    if (rank == size - 1){
    	j_last--;
    }

	/* Fill the interior with data */ 
    for (i=0; i<maxn; i++){ 
		for (j=j_first; j<=j_last; j++){
		    xloc[i][j] = (double)rank;
		}
	}
    /* Fill Boundaries */
    for (j=j_first; j<=j_last; j++){
		xloc[0][j] 				= 100.0;	// top
		xloc[maxn-1][j] 		= 0.0;	// bottom
	}
	if (rank==0){
		for (i=0; i<maxn; i++){
			xloc[i][j_first-1] 	= 20.0;	// left
		}
	}
	if (rank==3){
		for (i=0; i<maxn; i++){
			xloc[i][j_last+1] 	= 80.0;	// right
		}
	}
    
    /* All the Ghost points / places being recieved into are uninitialized */

	/* Print initial state */
	for (int k=0; k<size; k++){
		MPI_Barrier(MPI_COMM_WORLD);
		if (k==rank){
			printf("Rank: %d\n",rank);
			for (i=0; i<maxn; i++){
				for (j=j_first-1; j<=j_last+1; j++){
					printf("%2.1f ",xloc[i][j]);
				}
				printf("\n");
			}
		}
	}

    iter = 0; 
    do {
		/* Send right - if not final rank */
		if (rank < size - 1){
			for (i=0; i<maxn; i++){
				xedge[i] = xloc[i][j_last];
			}
			printf("Sending column %d from rank %d to rank %d\n", j_last, rank, rank+1);
		    MPI_Send( xedge, maxn, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD );
		}
		/* Recieve from left */
		if (rank > 0){
			printf("Recieving column %d from rank %d to rank %d\n", j_first-1, rank-1, rank);
		    MPI_Recv( xedge, maxn, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status );
		    for (i=0; i<maxn; i++){
		    	xloc[i][j_first-1] = xedge[i];
		    }
		}
		/* Send left */
		if (rank > 0){ 
			for (i=0; i<maxn; i++){
				xedge[i] = xloc[i][j_first];
			}
			printf("Sending column %d from rank %d to rank %d\n", j_first, rank, rank-1);
		    MPI_Send( xedge, maxn, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD );
		}
		/* Recieve from right */
		if (rank < size - 1){
			printf("Recieving column %d from rank %d to rank %d\n", j_last+1, rank+1, rank);
		    MPI_Recv( xedge, maxn, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status );
		    for (i=0; i<maxn; i++){
				xloc[i][j_last+1] = xedge[i];
			}
		}
		
		/* Compute new values (but not on boundary) */
		iter++;
		diffnorm = 0.0;
		for (i=1; i<maxn-1; i++) {
		    for (j=j_first; j<=j_last; j++) {
				xnew[i][j] = (xloc[i][j+1] + xloc[i][j-1] + xloc[i+1][j] + xloc[i-1][j]) / 4.0;
				diffnorm += (xnew[i][j] - xloc[i][j]) *  (xnew[i][j] - xloc[i][j]);
		    }
		}
		/* Only transfer the interior points */
		for (i=1; i<maxn-1; i++) {
		    for (j=j_first; j<=j_last; j++) {
				xloc[i][j] = xnew[i][j];
			}
		}
		MPI_Allreduce( &diffnorm, &gdiffnorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
		gdiffnorm = sqrt( gdiffnorm );
		if (rank == 0) {
			printf( "At iteration %d, diff is %e\n", iter, gdiffnorm );
		}
    } while ((gdiffnorm > 1e-3) && (iter < 1000));
	end = MPI_Wtime();
	if (rank==0){
		printf("\n");
		printf("Successful completion\n");
		printf("Converged after %d iterations\n", iter);
		printf("Total time: %f seconds\n", end - start);
	}
    MPI_Finalize( );
    return 0;
}
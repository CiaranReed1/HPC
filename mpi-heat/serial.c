/*
Adapted from Argonne National Lab MPI tutorial
https://www.mcs.anl.gov/research/projects/mpi/tutorial/mpiexmpl/contents.html
*/
#include <stdio.h>
#include <math.h>
#include <time.h>

/* This example handles a maxn x maxn mesh, on 4 processors only. */
#define maxn 24

int main( int argc, char **argv ){

	__clock_t start, end;
	start = clock();
    int        value, i, j, iter;
    int        i_first, i_last;
    double     diffnorm;
    double     xloc[maxn][maxn];	// Note: this is technically oversized for Rank 0 & 3
    double     xnew[maxn][maxn];

    /* Note that top and bottom processes have one less row of interior points */
    i_first = 1;
    i_last  = maxn-2;    // data runs up to maxn-1, interior is maxn-1-1
   
	/* Fill the interior with data */ 
    for (i=i_first; i<=i_last; i++){ 
		for (j=0; j<maxn; j++){
		    xloc[i][j]      = 0.0;
		}
	}
	/* Fill Boundaries */
	for (j=0; j<maxn; j++){
		xloc[i_first-1][j] 	= 100.0;	// top
	}
	for (j=0; j<maxn; j++){
		xloc[i_last+1][j] 	= 0.0;		// bottom
	}
	for (i=i_first; i<=i_last; i++){
		xloc[i][0] 			= 20.0;		// left
		xloc[i][maxn-1] 	= 80.0;		// right
	}
	

    iter = 0; 
    do {
        		
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
		diffnorm = sqrt( diffnorm );
		printf( "At iteration %d, diff is %e\n", iter, diffnorm );
    } while ((diffnorm > 1e-3) && (iter < 1000));
	end = clock();
	double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("\n");
	printf("Successful completion\n");
	printf("Converged after %d iterations\n", iter);
	printf("Total time: %f seconds\n", cpu_time_used);
    return 0;
}

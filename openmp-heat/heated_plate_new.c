#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

const int M = 512;
const int N = 512;
const double epsilon = 0.001;
const double hot = 100.0;
const double cold = 0.0;

void init(double unew[M][N]){
	int i,j;
	double mean;
	/* 
	Set the boundary values, which don't change. 
	*/
	for ( i = 1; i < M - 1; i++ ){
		unew[i][0] = hot;
	} //left side (not corners)
	for ( i = 1; i < M - 1; i++ ){
		unew[i][N-1] = hot;
	} // right side (not corners)
	for ( j = 0; j < N; j++ ){
		unew[M-1][j] = hot;
	} //bottom side (whole side)
	for ( j = 0; j < N; j++ ){
		unew[0][j] = cold;
	} //top side (whole side)

	/*
	Average the boundary values, to come up with a reasonable
	initial value for the interior.
	*/
	mean = 0.0;
	for ( i = 1; i < M - 1; i++ ){
		mean = mean + unew[i][0];
	}
	for ( i = 1; i < M - 1; i++ ){
		mean = mean + unew[i][N-1];
	}
	for ( j = 0; j < N; j++ ){
		mean = mean + unew[M-1][j];
	}
	for ( j = 0; j < N; j++ ){
		mean = mean + unew[0][j];
	}
	mean = mean / ( double ) ( 2 * M + 2 * N - 4 );
	
	/* 
	Initialize the interior solution to the mean value.
	*/
	for ( i = 1; i < M - 1; i++ ){
		for ( j = 1; j < N - 1; j++ ){
	  		unew[i][j] = mean;
		}
	}	
}

void copy(double u[M][N], double unew[M][N]){
	int i,j;
	/*
	Save the old solution in U.
	*/
	for ( i = 0; i < M; i++ ){
		for ( j = 0; j < N; j++ ){
			u[i][j] = unew[i][j];
		}
	}
}

void saveoutput(double unew[M][N], char output_file[80]){
	
	int i,j;
	FILE *fp;
	
	fp = fopen ( output_file, "w" );
	for ( i = 0; i < M; i++ ){
		for ( j = 0; j < N; j++){
			fprintf(fp, "%6.2f ", unew[i][j] );
		}
		fprintf(fp, "\n");
	}
	fclose( fp );
	printf( "\n" );
	printf("  Solution written to the output file '%s'\n", output_file );
}

int main ( int argc, char *argv[] ){
	/*
	Adapted from the code published by John Burkardt at FSU: https://people.sc.fsu.edu/~jburkardt/c_src/heated_plate/heated_plate.html
	*/
	double maxdiff;
	int i, j, iter, converged, printiter = 100;
	double u[M][N];
	double unew[M][N];

	printf ( "\n" );
	printf ( "HEATED_PLATE\n" );
	printf ( "  C version\n" );
	printf ( "  A program to solve for the steady state temperature distribution\n" );
	printf ( "  over a rectangular plate.\n" );
	printf ( "\n" );
	printf ( "  Spatial grid of %d by %d points.\n", M, N );
	
	printf ( "\n" );
	printf ( "  The iteration will be repeated until the change is <= %f\n", epsilon );
	maxdiff = 2.0*epsilon;
	clock_t start_t = clock();
	// initialize the unew field
	init(unew);

	// iterate until the unew solution differs from the old solution u by no more than EPSILON.
	iter = 0;
	printf ( "\n" );
	printf ( " Iteration  Change\n" );
	printf ( "\n" );
	
	while ( epsilon <= maxdiff ){
		
		// store the old solution in u
		copy(u, unew);
		
		/*
		Determine the new estimate of the solution at the interior points.
		The new solution W is the average of north, south, east and west neighbors.
		*/
		maxdiff = 0.0;
		for ( i = 1; i < M - 1; i++ ){
			for ( j = 1; j < N - 1; j++ ){
				unew[i][j] = ( u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] ) / 4.0;
				if (maxdiff < fabs ( unew[i][j] - u[i][j] )){
					maxdiff = fabs ( unew[i][j] - u[i][j] );
				}
			}
		}
		iter++;
		if (iter%printiter == 0){
			printf ( "  %8d  %f\n", iter, maxdiff );
		}
	}
	clock_t end_t = clock();
	double total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
	printf ( "\n" );
	printf ( "  %8d  %f\n", iter, maxdiff );
	printf ( "\n" );
	printf ( "  Error tolerance achieved.\n" );
	printf ( "  Total time = %f seconds.\n", total_t );
	/* 
	Terminate. 
	*/
	printf ( "\n" );
	printf ( "HEATED_PLATE:\n" );
	printf ( "  Normal end of execution.\n" );
	printf ( "\n" );
	
	return 0;
}

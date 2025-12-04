#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>


const double epsilon = 5e-6;
const double hot = 100.0;
const double cold = 0.0;

void init(double *unew,int N){
	double x;
	unew[0] = cold;
	#pragma omp parallel for default(none) private(x) shared(unew,N,hot,cold)
	for (int i = 1; i < N-1; i++){
		x = (float)(i)/(float)(N);
		unew[i] = cold + (hot-cold)*0.5*(1.0 + tanh(10.0*(x - 0.5)));
	}
	unew[N-1] = hot;
}

void step(double *unew, double *u, int N){
	unew[0] = cold;		//unew[0] 	= (u[0] + u[1])/2.0;			// alternative boundary conditions
	#pragma omp parallel for default(none) shared(N,unew,u)
	for (int i=1; i < N-1; i++){
		unew[i] = (u[i-1] + u[i+1])/2.0;
	}
	unew[N-1] = hot;	//unew[N-1]	= (u[N-1] + u[N-2])/2.0;	// alternative boundary conditions
}

void copy(double *unew, double *u, int N){
	#pragma omp parallel for default(none) shared(N,unew,u)
 	for (int i=0; i < N; i++){
		unew[i] = u[i];
	}
}

double diff(double *unew, double *u, int N){
	double maxdiff = 0.0;
	#pragma omp parallel for default(none) shared(N,unew,u) reduction(max:maxdiff)
	for (int i=0; i < N-1; i++ ){
		maxdiff = fmax(maxdiff, fabs(unew[i] - u[i]));
	}
	return maxdiff;
}

int main(int argc, char *argv[]){
	int nThreads = argc > 1 ? atoi(argv[1]) : 1;
	omp_set_num_threads(nThreads);
	int N = argc > 2 ? atoi(argv[2]) : 1<<10;
	//#pragma omp parallel
	//{
	//	printf("Thread %d of %d says hello!\n", omp_get_thread_num(), omp_get_num_threads());
	//}
	double maxdiff = 2.0*epsilon;
	double *u, *unew;
	u = (double*) malloc(N * sizeof(double));
	unew = (double*) malloc(N * sizeof(double));
	int iter = 0;

	init(unew, N);
	double start_t = omp_get_wtime();
	double copy_time = 0.0;
	while (maxdiff >= epsilon){
		double c_start = omp_get_wtime();
		copy(u, unew, N);
		double c_end = omp_get_wtime();
		copy_time += c_end - c_start;
		step(unew, u, N);
		maxdiff = diff(unew, u, N);
		iter++;
	}
	copy_time /= iter;
	double end_t = omp_get_wtime();
	double total_t = end_t - start_t;
	printf("%d,%d,%2.16f,%2.16f,%d\n",nThreads,N, total_t, copy_time, iter);
	free(u);
	free(unew);
	return 0;
}

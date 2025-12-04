#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

const int N = 1<<10;
const double epsilon = 5e-6;
const double hot = 100.0;
const double cold = 0.0;

void init(double unew[N]){
	double x;
	unew[0] = cold;
	for (int i = 1; i < N-1; i++){
		x = (float)(i)/(float)(N);
		unew[i] = cold + (hot-cold)*0.5*(1.0 + tanh(10.0*(x - 0.5)));
	}
	unew[N-1] = hot;
}

void step(double unew[N], double u[N]){
	unew[0] = cold;		//unew[0] 	= (u[0] + u[1])/2.0;			// alternative boundary conditions
	for (int i=1; i < N-1; i++){
		unew[i] = (u[i-1] + u[i+1])/2.0;
	}
	unew[N-1] = hot;	//unew[N-1]	= (u[N-1] + u[N-2])/2.0;	// alternative boundary conditions
}

void copy(double unew[N], double u[N]){
 	for (int i=0; i < N; i++){
		unew[i] = u[i];
	}
}

double diff(double unew[N], double u[N]){
	double maxdiff = 0.0;
	for (int i=0; i < N-1; i++ ){
		if (maxdiff < fabs(unew[i]-u[i])){
			maxdiff = fabs(unew[i]-u[i]);
		}
	}
	return maxdiff;
}

int main(){
	double maxdiff = 2.0*epsilon;
	double u[N], unew[N];
	int iter = 0;

	init(unew);
	clock_t start_t = clock();
	while (maxdiff >= epsilon){
		copy(u, unew);
		step(unew, u);
		maxdiff = diff(unew, u);
		iter++;
	}
	clock_t end_t = clock();
	double total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
	printf("[%d] Iteration time = %2.16fs (%d iterations).\n",N, total_t, iter);
}

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double time_ij(int N){

	int i, j;
	
	double A[N][N];
	double B[N][N];
	double C[N][N];
	
	clock_t t0 = clock();
	for(i=0; i<N; i++){
		for (j=0; j<N; j++){
			C[i][j] = A[i][j] + B[i][j];
		}
	}
	double t1 = ((double)clock() - t0) / CLOCKS_PER_SEC;
	return t1;
}

double time_ji(int N){

	int i, j;
	
	double A[N][N];
	double B[N][N];
	double C[N][N];
	
	clock_t t0 = clock();
	for(j=0; j<N; j++){
		for (i=0; i<N; i++){
			C[i][j] = A[i][j] + B[i][j];
		}
	}
	double t1 = ((double)clock() - t0) / CLOCKS_PER_SEC;
	return t1;
}

double time_ijk(int N){
	int i,j,k;
	
	double A[N][N][N];
	double B[N][N][N];
	double C[N][N][N];
	
	clock_t t0 = clock();
	#pragma omp parallel for 
	for(i=0; i<N; i++){
		for (j = 0;j<N;j++){
			for(k = 0;k<N;k++){
				C[i][j][k] = A[i][j][k] + B[i][j][k];
				}
		}
	}
	double t1 = ((double)clock() - t0)/CLOCKS_PER_SEC;
	return t1;
}

double time_ikj(int N){
	int i,j,k;
    
	double A[N][N][N];
	double B[N][N][N];
	double C[N][N][N];
    
	clock_t t0 = clock();
	#pragma omp parallel for
	for(i=0; i<N; i++){
		for(k=0; k<N; k++){
			for(j=0; j<N; j++){
				C[i][j][k] = A[i][j][k] + B[i][j][k];
			}
		}
	}
	double t1 = ((double)clock() - t0)/CLOCKS_PER_SEC;
	return t1;
}

double time_jik(int N){
	int i,j,k;
    
	double A[N][N][N];
	double B[N][N][N];
	double C[N][N][N];
    
	clock_t t0 = clock();
	#pragma omp parallel for
	for(j=0; j<N; j++){
		for(i=0; i<N; i++){
			for(k=0; k<N; k++){
				C[i][j][k] = A[i][j][k] + B[i][j][k];
			}
		}
	}
	double t1 = ((double)clock() - t0)/CLOCKS_PER_SEC;
	return t1;
}

double time_jki(int N){
	int i,j,k;
    
	double A[N][N][N];
	double B[N][N][N];
	double C[N][N][N];
    
	clock_t t0 = clock();
	#pragma omp parallel for
	for(j=0; j<N; j++){
		for(k=0; k<N; k++){
			for(i=0; i<N; i++){
				C[i][j][k] = A[i][j][k] + B[i][j][k];
			}
		}
	}
	double t1 = ((double)clock() - t0)/CLOCKS_PER_SEC;
	return t1;
}

double time_kij(int N){
	int i,j,k;
    
	double A[N][N][N];
	double B[N][N][N];
	double C[N][N][N];
    
	clock_t t0 = clock();
	#pragma omp parallel for
	for(k=0; k<N; k++){
		for(i=0; i<N; i++){
			for(j=0; j<N; j++){
				C[i][j][k] = A[i][j][k] + B[i][j][k];
			}
		}
	}
	double t1 = ((double)clock() - t0)/CLOCKS_PER_SEC;
	return t1;
}

double time_kji(int N){
	int i,j,k;
    
	double A[N][N][N];
	double B[N][N][N];
	double C[N][N][N];
    
	clock_t t0 = clock();
	#pragma omp parallel for
	for(k=0; k<N; k++){
		for(j=0; j<N; j++){
			for(i=0; i<N; i++){
				C[i][j][k] = A[i][j][k] + B[i][j][k];
			}
		}
	}
	double t1 = ((double)clock() - t0)/CLOCKS_PER_SEC;
	return t1;
}

typedef struct {
	const char *name;
	double time;
} perm_result_t;

void benchmark_3d_permutations(int N, int M){
	const char *names[6] = {"ijk", "ikj", "jik", "jki", "kij", "kji"};
	double (*funcs[6])(int) = {time_ijk, time_ikj, time_jik, time_jki, time_kij, time_kji};
	perm_result_t results[6];

	for(int p = 0; p < 6; p++){
		double t = 0.0;
		for(int m = 0; m < M; m++){
			t += funcs[p](N);
		}
		results[p].name = names[p];
		results[p].time = t / M;
	}

	/* simple ascending sort (selection) by time */
	for(int i = 0; i < 6; i++){
		int best = i;
		for(int j = i+1; j < 6; j++){
			if(results[j].time < results[best].time) best = j;
		}
		if(best != i){
			perm_result_t tmp = results[i];
			results[i] = results[best];
			results[best] = tmp;
		}
	}

	printf("3D permutations (fastest -> slowest):\n");
	for(int i = 0; i < 6; i++){
		printf("\t%s: %2.16f s\n", results[i].name, results[i].time);
	}
}


int main(int argc, char **argv){

	int N = atoi(argv[1]);
	int M = atoi(argv[2]);
	double t;
	
	printf("3D version:\n");
	
	/*
		Your 3D code here
	*/
	benchmark_3d_permutations(N, M);
	
	printf("2D version:\n");
	
	t=0.0;
	for (int m=0; m < M; m++){
		t += time_ij(N);
	}
	printf("\tij: %2.16f s\n",t/M);
	
	t=0.0;
	for (int m=0; m < M; m++){
		t += time_ij(N);
	}
	printf("\tji: %2.16f s\n",t/M);

	
	return 0;
}

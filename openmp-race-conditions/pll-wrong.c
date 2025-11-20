#include <stdio.h>
#include <stdint.h>
#include <omp.h>

int64_t sum_pll(int N, int n_threads){
	int64_t sum = 0; 
	#pragma omp parallel for num_threads(n_threads)
	for (int n=1; n<N+1; n++){
		sum+=n;
	}
	return sum;
}	

int main(){
	const int iterations = 1000;
	int threads[8] = {1, 2,3, 4,5,6,7, 8};
	int N[6] = {10,100,1000,10000,100000,1000000};
	int64_t desired[6]; 
	double scores[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	for (int i = 0;i < 8; i++){
		double start_time = omp_get_wtime();
		int n_threads = threads[i];
		double score = 0.0;
		for (int n_idx = 0; n_idx < 6; n_idx++){
			int n = N[n_idx];
			double n_score = 0.0;
			desired[n_idx] = (int64_t)n*(n+1)/2;
			for (int iter=0; iter<iterations; iter++){
				if (sum_pll(n, n_threads) == desired[n_idx]){
					n_score++;
				}
			}
			score += n_score / iterations;
			printf("Using %d threads, for n=%d, score is %f.\n", n_threads, n, n_score / iterations);
		}
		scores[i] = score / 6.0;
		double end_time = omp_get_wtime();
		printf("Using %d threads, average scores across range of n is %f.\n", n_threads, scores[i]);
		printf("Time taken with %d threads: %f seconds.\n", n_threads, end_time - start_time);
	}	
	return 0;
}

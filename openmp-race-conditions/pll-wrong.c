#include <stdio.h>


int sum_pll(int N, int n_threads){
	int sum = 0;
	#pragma omp parallel for num_threads(n_threads)
	for (int n=1; n<N+1; n++){
		sum+=n;
	}
	return sum;
}	

int main(){
	const int N = 100;
	const int desired = N*(N+1)/2;
	const int iterations = 1000;
	int threads[8] = {1, 2,3, 4,5,6,7, 8};
	double scores[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	for (int i = 0;i < 8; i++){
		int n_threads = threads[i];
		double score = 0.0;
		for (int iter=0; iter<iterations; iter++){
			if (sum_pll(N, n_threads) == desired){
				score++;
			}
		}
		scores[i] = score/iterations;
		printf("Using %d threads, score is %f.\n", n_threads, scores[i]);
	}	
	return 0;
}

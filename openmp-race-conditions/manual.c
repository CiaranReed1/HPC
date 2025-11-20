#include <stdio.h>
#include <omp.h>
int main(){
    const int N = 40000;
    int partialsums[4] = {0,0,0,0};
    #pragma omp parallel for num_threads(4) default(none) shared(partialsums,N)
    for (int n=1; n<N+1; n++){
        partialsums[omp_get_thread_num()]+=n;
    }
    int sum = 0;
    for (int n=0; n<4; n++){
        sum += partialsums[n];
    }
    printf("Result is %d. It should be %d.\n",
        sum, N*(N+1)/2);
    return 0;
}
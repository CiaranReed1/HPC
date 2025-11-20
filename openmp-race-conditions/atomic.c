#include <stdio.h>
#include <omp.h>
int main(){
    const int N = 40000;
    int sum =0;
    int _sum;
    #pragma omp parallel num_threads(8) default(none) shared(N,sum) private(_sum)
    {
        _sum = 0;
        #pragma omp for
        for (int n=1; n<N+1; n++){
            _sum+=n;
        }
        #pragma omp atomic
        sum += _sum;
    }
    printf("Result is %d. It should be %d.\n",
        sum, N*(N+1)/2);
    return 0;
}
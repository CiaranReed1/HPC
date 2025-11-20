#include <stdio.h>

int main(){
	const int N = 40000;
	int sum = 0;

	for (int n=1; n<N+1; n++){
		sum+=n;
	}

	printf("Result is %d. It should be %d.\n", sum, N*(N+1)/2);
	return 0;
}

#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <stdlib.h>


double decode_message(char *encoded_message, int nthreads,char *decoded_message) {
    double start_time = omp_get_wtime();
    omp_set_num_threads(nthreads);
    int length = strlen(encoded_message);
    #pragma omp parallel for
    for (int i = 0; i < length; i++) {
        decoded_message[i] = encoded_message[i] - 1;
    }
    decoded_message[length] = '\0'; // Null-terminate the decoded string
    double end_time = omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp single
        printf("Using %d threads , Decoded Message: %s\n", omp_get_num_threads(), decoded_message);
    }
    //printf("Decoding Time with %d threads: %f seconds\n", nthreads, end_time - start_time);
    return end_time - start_time;
}


int main(int argc, char *argv[]) {
    char encoded_message[256];
    strcpy(encoded_message, argv[1]);
    int threads[4] = {1, 2, 4,8};
    double times[4];
    for (int i = 0; i < 4; i++) {
        char *decoded_message = (char *)malloc((strlen(encoded_message) + 1) * sizeof(char));
        times[i] = decode_message(encoded_message, threads[i], decoded_message);
        free(decoded_message);
    }
    printf("Decoding times:\n");
    for (int i = 0; i < 4; i++) {
        printf("%d threads: %f seconds\n", threads[i], times[i]);
    }
    return 0;
}


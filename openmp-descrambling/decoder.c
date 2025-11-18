#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <stdlib.h>

int main(int *argc, char *argv[]) {
    char encoded_message[256];
    strcpy(encoded_message, argv[1]);
    int length = strlen(encoded_message);
    #pragma omp parallel for
    for (int i = 0; i < length; i++) {
        encoded_message[i] = encoded_message[i] - 1;
         printf("Thread Number: %d is decoding character position: %d, %c to %c", omp_get_thread_num(), i, encoded_message[i] + 1, encoded_message[i]);
    }
    printf("Decoded Message: %s\n", encoded_message);
    return 0;
}
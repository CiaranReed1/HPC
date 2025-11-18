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
        decoded_message[i] = encoded_message[i] - 1 - omp_get_thread_num();
    }
    decoded_message[length] = '\0'; // Null-terminate the decoded string
    double end_time = omp_get_wtime();
    // do not print decoded message here; return elapsed time to caller
    //printf("Decoding Time with %d threads: %f seconds\n", nthreads, end_time - start_time);
    return end_time - start_time;
}


int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <nthreads> <input_file>\n", argv[0]);
        return 1;
    }

    int nthreads = atoi(argv[1]);
    if (nthreads <= 0) {
        fprintf(stderr, "Invalid number of threads: %s\n", argv[1]);
        return 1;
    }

    const char *in_filename = argv[2];
    FILE *in = fopen(in_filename, "r");
    if (!in) {
        perror("fopen");
        fprintf(stderr, "Could not open input file '%s' for reading\n", in_filename);
        return 1;
    }

    // Read entire file into buffer
    if (fseek(in, 0, SEEK_END) != 0) {
        perror("fseek");
        fclose(in);
        return 1;
    }
    long fsize = ftell(in);
    if (fsize < 0) fsize = 0;
    rewind(in);

    char *filebuf = (char *)malloc(fsize + 1 + 1);
    if (!filebuf) {
        fprintf(stderr, "Out of memory reading input file\n");
        fclose(in);
        return 1;
    }
    size_t read = fread(filebuf, 1, fsize, in);
    filebuf[read] = '\0';
    fclose(in);

    // Treat entire file contents as the encoded message. Trim trailing newlines.
    char *message = filebuf;
    size_t msglen = strlen(message);
    while (msglen > 0 && (message[msglen-1] == '\n' || message[msglen-1] == '\r')) {
        message[msglen-1] = '\0';
        msglen--;
    }

    // Allocate decoded buffer
    size_t mlen = strlen(message);
    char *decoded_message = (char *)malloc(mlen + 1);
    if (!decoded_message) {
        fprintf(stderr, "Out of memory allocating decoded buffer\n");
        free(filebuf);
        return 1;
    }

    // Decode using requested thread count and capture elapsed time
    double elapsed = decode_message(message, nthreads, decoded_message);

    // Write decoded message to output file named decoded_<nthreads>_threads.txt
    char outname[64];
    snprintf(outname, sizeof(outname), "decoded_%d_threads.txt", nthreads);
    FILE *out = fopen(outname, "w");
    if (!out) {
        perror("fopen");
        fprintf(stderr, "Could not open output file '%s' for writing\n", outname);
        free(filebuf);
        free(decoded_message);
        return 1;
    }
    if (fprintf(out, "%s\n", decoded_message) < 0) {
        perror("fprintf");
        fclose(out);
        free(filebuf);
        free(decoded_message);
        return 1;
    }
    fclose(out);

    printf("Decoded message written to %s\n", outname);
    printf("Decoding time with %d threads: %f seconds\n", nthreads, elapsed);

    free(filebuf);
    free(decoded_message);
    return 0;
}


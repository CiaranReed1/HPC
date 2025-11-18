#include <stdio.h>  // Import the printf function
#include <omp.h>    // Import the OpenMP library functions
#include <string.h>
#include <stdlib.h>

char* encoder(char *message,int nthreads){
   int length = strlen(message);
   char *encoded_message = (char *)malloc((length + 1) * sizeof(char));
   omp_set_num_threads(nthreads);
   #pragma omp parallel for
   for (int i = 0; i < length; i++){
      encoded_message[i] = message[i] + omp_get_thread_num() + 1;
   }
   encoded_message[length] = '\0'; // Null-terminate the encoded string
   return encoded_message;
}

// Start of program
int main(int argc, char* argv[]){
   double start_time = omp_get_wtime();
   if (argc < 3) {
      fprintf(stderr, "Usage: %s <nthreads> <message>\n", argv[0]);
      return 1;
   }

   /* parse number of threads */
   char *endptr = NULL;
   long nthreads = strtol(argv[1], &endptr, 10);
   if (endptr == argv[1] || *endptr != '\0' || nthreads <= 0) {
      fprintf(stderr, "Invalid number of threads: %s\n", argv[1]);
      return 1;
   }

   /* copy message (argv[2] expected to be quoted if containing spaces) */
   char mymessage[1024];
   strncpy(mymessage, argv[2], sizeof(mymessage) - 1);
   mymessage[sizeof(mymessage) - 1] = '\0';

   /* encode with requested number of threads */
   char *encoded = encoder(mymessage, (int)nthreads);
   if (!encoded) {
      fprintf(stderr, "Encoding failed\n");
      return 1;
   }

   /* write to file named encoded_<n>_threads.txt */
   char out_filename[64];
   snprintf(out_filename, sizeof(out_filename), "encoded_%ld_threads.txt", nthreads);
   FILE *out = fopen(out_filename, "w");
   if (!out) {
      perror("fopen");
      fprintf(stderr, "Could not open output file '%s' for writing\n", out_filename);
      free(encoded);
      return 1;
   }

   if (fprintf(out, "%s\n", encoded) < 0) {
      perror("fprintf");
      fclose(out);
      free(encoded);
      return 1;
   }

   fclose(out);
   double end_time = omp_get_wtime();
   printf("Encoding Time with %ld threads: %f seconds\n", nthreads, end_time - start_time);
   printf("Encoded message written to %s\n", out_filename);
   free(encoded);
   return 0; // Exit program
}


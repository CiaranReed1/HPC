#include <stdio.h>  // Import the printf function
#include <omp.h>    // Import the OpenMP library functions
#include <string.h>
#include <stdlib.h>

// Stream-encode a buffer using OpenMP and the specified number of threads.
static void encode_buffer(unsigned char *inbuf, unsigned char *outbuf, size_t len, int nthreads) {
    omp_set_num_threads(nthreads);
    #pragma omp parallel for
    for (size_t i = 0; i < len; ++i) {
        outbuf[i] = inbuf[i] + omp_get_thread_num() + 1;
    }
}

// Start of program
int main(int argc, char* argv[]){
   double start_time = omp_get_wtime();
   if (argc < 3) {
      fprintf(stderr, "Usage: %s <nthreads> <input_file>\n", argv[0]);
      return 1;
   }

   /* parse number of threads */
   char *endptr = NULL;
   long nthreads = strtol(argv[1], &endptr, 10);
   if (endptr == argv[1] || *endptr != '\0' || nthreads <= 0) {
      fprintf(stderr, "Invalid number of threads: %s\n", argv[1]);
      return 1;
   }

   const char *in_filename = argv[2];
   FILE *in = fopen(in_filename, "rb");
   if (!in) {
      perror("fopen");
      fprintf(stderr, "Could not open input file '%s' for reading\n", in_filename);
      return 1;
   }

   char out_filename[128];
   snprintf(out_filename, sizeof(out_filename), "encoded_%ld_threads.txt", nthreads);
   FILE *out = fopen(out_filename, "wb");
   if (!out) {
      perror("fopen");
      fprintf(stderr, "Could not open output file '%s' for writing\n", out_filename);
      fclose(in);
      return 1;
   }

   // Read the entire input file into memory (user requested). This uses more RAM
   // but simplifies processing for cases where the user prefers one-shot reads.
   if (fseek(in, 0, SEEK_END) != 0) {
      perror("fseek");
      fclose(in);
      fclose(out);
      return 1;
   }
   long fsize = ftell(in);
   if (fsize < 0) fsize = 0;
   rewind(in);

   size_t usize = (size_t)fsize;
   // allow encoding of empty files
   unsigned char *inbuf = (unsigned char *)malloc(usize > 0 ? usize : 1);
   unsigned char *outbuf = (unsigned char *)malloc(usize > 0 ? usize : 1);
   if (!inbuf || !outbuf) {
      fprintf(stderr, "Out of memory allocating buffers\n");
      fclose(in);
      fclose(out);
      free(inbuf);
      free(outbuf);
      return 1;
   }

   size_t nread = fread(inbuf, 1, usize, in);
   if (nread != usize) {
      if (ferror(in)) perror("fread");
      // proceed with the number of bytes actually read
      usize = nread;
   }

   // encode entire buffer
   if (usize > 0) encode_buffer(inbuf, outbuf, usize, (int)nthreads);

   size_t nwritten = fwrite(outbuf, 1, usize, out);
   if (nwritten != usize) {
      perror("fwrite");
      fprintf(stderr, "Failed to write encoded data to '%s'\n", out_filename);
      free(inbuf);
      free(outbuf);
      fclose(in);
      fclose(out);
      return 1;
   }

   free(inbuf);
   free(outbuf);
   fclose(in);
   fclose(out);

   double end_time = omp_get_wtime();
   printf("Encoding time with %ld threads: %f seconds\n", nthreads, end_time - start_time);
   printf("Encoded file written to %s\n", out_filename);
   return 0; // Exit program
}


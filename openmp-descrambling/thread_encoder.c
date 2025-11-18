#include <stdio.h>  // Import the printf function
#include <omp.h>    // Import the OpenMP library functions
#include <string.h>
#include <stdlib.h>

char* encoder(char *message,int nthreads){
   double start_time = omp_get_wtime();
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

   // This is a "character array", consisting of 105 letters in byte form.
   char mymessage[256];
   strcpy(mymessage, argv[1]); 
   // The message is encoded by shifting all letters by one, i.e. 
   // A becomes B, c becomes d and so forth. This shift is easy to implement
   // in byte form as all you need to do is to add 1 to a character.
   // To printf a single character, use e.g. 
   //       printf("%c\n", mymessage[0]);
   // To printf the whole array do e.g.   (subtle difference between %c and %s)
   //       printf("%s\n", mymessage);

   /* MISSING CODE: for loop over the 105 characters and adding the value 1 to each entry; */
   int nthreads[4] = {1,2,4,8};
   char encoded_messages[4][256];
   for (int i = 0; i< 4; i++){
      char *result = encoder(mymessage,nthreads[i]);
      strcpy(encoded_messages[i], result);
      free(result);
      printf("Encoded message with %d threads: %s\n", nthreads[i], encoded_messages[i]);
   }
   return 0; // Exit program
}


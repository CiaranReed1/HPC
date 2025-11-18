#include <stdio.h>  // Import the printf function
#include <omp.h>    // Import the OpenMP library functions
#include <string.h>

// Start of program
int main(){

   // This is a "character array", consisting of 105 letters in byte form.
   char mymessage[] = "This is a placeholder of the original message";
   // The message is encoded by shifting all letters by one, i.e. 
   // A becomes B, c becomes d and so forth. This shift is easy to implement
   // in byte form as all you need to do is to add 1 to a character.
   // To printf a single character, use e.g. 
   //       printf("%c\n", mymessage[0]);
   // To printf the whole array do e.g.   (subtle difference between %c and %s)
   //       printf("%s\n", mymessage);

   /* MISSING CODE: for loop over the 105 characters and adding the value 1 to each entry; */
   int length = strlen(mymessage);
   for (int i = 0; i < length; i++){
      mymessage[i] = mymessage[i] + 1;
   }
   printf("%s\n", mymessage); // This prints a character array to screen

   return 0; // Exit program
}
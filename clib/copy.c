

#include <stdio.h>
#include <stdlib.h>

int copy(const char* fn_in, const char* fn_out) {
   
   FILE *fd1 = fopen(fn_in, "r");
   if(!fd1) return -1; // cannot open input
     
   FILE *fd2 = fopen(fn_out, "w");
   if(!fd2) {  // cannot open output
     fclose(fd1);
     return -2;
   }

   size_t l1;
   unsigned char buffer[8192]; 

   while((l1 = fread(buffer, 1, sizeof buffer, fd1)) > 0) {
     size_t l2 = fwrite(buffer, 1, l1, fd2);
     if(l2 == 0 || l2 < l1) {
       fclose(fd1);
       fclose(fd2);
       if(l2==0) return -3; // output error
       return -4; // disk full
     }
   }
   fclose(fd1);
   fclose(fd2);
   return 0;
}


#include "c_defs.h"
/* 
  Print the stack trace
*/
#ifdef __PTRACE
#include <execinfo.h>
#include <stdlib.h>
#endif
#include <stdio.h>

void F77_FUNC(ptrace,PTRACE)(int *kilobytes)
{
#ifdef __PTRACE
  void *array[12];
  size_t size;
  char **strings;
  size_t i;

  size = backtrace (array, 12);
  strings = backtrace_symbols (array, size);

  printf ("Obtained %zd stack frames.\n", size);
  printf ("Use 'addr2line -e /where/is/code.x 0x12345' to get the source line number\n");

  for (i = 0; i < size; i++)
     printf ("%s\n", strings[i]);

  free (strings);
#else
  printf ("No stack trace available.\n");
#endif
}


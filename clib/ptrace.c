#include "c_defs.h"
/* 
  Print the stack trace
*/
#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>

void F77_FUNC(ptrace,PTRACE)(int *kilobytes)
{
#ifdef __PTRACE
  void *array[10];
  size_t size;
  char **strings;
  size_t i;

  size = backtrace (array, 10);
  strings = backtrace_symbols (array, size);

  printf ("Obtained %zd stack frames.\n", size);

  for (i = 0; i < size; i++)
     printf ("%s\n", strings[i]);

  free (strings);
#else
  printf ("No stack trace available.\n");
#endif
}


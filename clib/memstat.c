/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include <stdlib.h>
#include <malloc.h>
#include "c_defs.h"

void MEMSTAT(int *me)
{
#if defined __ALPHA
  /* printf("\n   MEMORY STATISTICS FROM PE(%d) : \n",*me); */
#endif
#if defined __LINUX
  printf("\n   MEMORY STATISTICS FROM PE(%d) : \n",*me);
#endif
#if defined __T3E
  printf("\n   MEMORY STATISTICS FROM PE(%d) : \n",*me);
  malloc_stats (0);
#endif
#if defined __SGI || defined __ORIGIN || defined __ALTIX
  printf("\n   MEMORY STATISTICS FROM PE(%d) : \n",*me);
  /* malloc_stats (0); */
#endif
#if defined __AIX
  struct mallinfo info;  info = mallinfo();
  printf("\n   PE(%d) MEMORY STATISTICS ",*me);
  printf("\n   total space in arena :%d\n",info.arena);
#endif
}

void print_ptr( void * ptr ) {
  printf("\n   PRTVAL: %p", ptr );
}

#if defined __T3E

void MEMORY_STATISTICS(int *me)
{

    struct mallinfo Memory;
    Memory = mallinfo();

    if((*me)==1) {
      printf("\n\n ==========================================\n");
      printf(" MEMORY STATISTICS\n");
      printf(" total space in arena            %10d\n",Memory.arena);
      printf(" number of ordinary blocks       %10d\n",Memory.ordblks);
      printf(" number of small blocks          %10d\n",Memory.smblks);
      printf(" number of holding blocks        %10d\n",Memory.hblks);
      printf(" space in holding block headers  %10d\n",Memory.hblkhd);
      printf(" space in small blocks in use    %10d\n",Memory.usmblks);
      printf(" space in free small blocks      %10d\n",Memory.fsmblks);
      printf(" space in ordinary blocks in use %10d\n",Memory.uordblks);
      printf(" space in free ordinary blocks   %10d\n",Memory.fordblks);
      printf(" cost of enabling keep option    %10d\n",Memory.keepcost);
      printf(" ==========================================\n\n");
    }

}

void FREE_MEMORY(int *free)
{
    static int total_memory = 16000000;
    struct mallinfo Memory;
    Memory = mallinfo();
    *free = total_memory - Memory.arena/8;
}
 
#endif

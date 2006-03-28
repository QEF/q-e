/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include <stdlib.h>

#if !defined(__MAC)
#include <malloc.h>
#endif

#include <stdio.h>

#include "c_defs.h"

void MEMSTAT(int *me)
{
  /* T3E: malloc_stats (0); */
#if defined __AIX
  struct mallinfo info;  info = mallinfo();
  printf("\n   PE(%d) MEMORY STATISTICS ",*me);
  printf("\n   total space in arena :%ld\n",info.arena);
#endif
}


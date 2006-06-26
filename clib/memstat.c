/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#if !defined(__MAC)
#include <malloc.h>
#endif

#include "c_defs.h"

/* 
  This function return the numer of kilobytes allocated
  by the calling process. 
  Auhor: Carlo Cavazzoni.
*/

void MEMSTAT(int *kilobytes)
{
#if defined __AIX || defined __LINUX || defined __LINUX64
  struct mallinfo info;  
  info = mallinfo();
  *kilobytes = info.arena / 1024 ;
#else
  *kilobytes = -1;
#endif
}


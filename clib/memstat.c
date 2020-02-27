/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

/* 
  This function returns the number of kilobytes
  allocated by the calling process. 
  Author: Carlo Cavazzoni.
  Obsolete AIX case and F77-C binding removed by P. Giannozzi (2017)
*/

#include "qe_cdefs.h"

#if defined (__SVR4) && defined (__sun)
#define SUN_MALLINFO
#endif

#if defined(HAVE_MALLINFO) && !defined(__QK_USER__) && !defined(SUN__MALLINFO) 

#include <malloc.h>
int c_memstat( )
{
  struct mallinfo info;  
  info = mallinfo();
  return (info.arena + info.hblkhd) / 1024 ;
#else
int c_memstat( )
{
  return -1;
#endif
}

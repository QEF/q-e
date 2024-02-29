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
  Obsolete cases and F77-C binding removed by P. Giannozzi (2017)
*/

#if defined(HAVE_MALLINFO)

#include <malloc.h>
int c_memstat( )
{
#if defined(__GLIBC__) && (__GLIBC__ > 2 || (__GLIBC__ == 2  && __GLIBC_MINOR__ >= 33))
  struct mallinfo2 info;
  info = mallinfo2();
#else
  struct mallinfo info;
  info = mallinfo();
#endif
  return (info.arena + info.hblkhd) / 1024 ;
#else
int c_memstat( )
{
  return -1;
#endif
}

/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include<stdio.h>
#include<time.h>
#include<ctype.h>
#include<sys/types.h>
#include<sys/time.h>

#include"c_defs.h"


double ELAPSED_SECONDS()
{
  static time_t tstart, tend;
  static int first = 1;
  double sec;
  time(&tend);
  if( first ) {
    tstart = tend;
    first = 0;
  }
  sec = difftime( tend, tstart );
  return sec;
}


double CCLOCK()
/* Restituisce i secondi trascorsi dalla chiamata al timer rest */
{

#if defined __T3E

/*  return (double)(rtclock() * 3.333e-6 / 2.); */
    return (double)( ( _rtc() / (double)CLK_TCK ) );

#else

    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;

#endif

}

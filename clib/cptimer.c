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
#include"cp.h"


static struct tms {
        clock_t tms_utime;              /* user time */
        clock_t tms_stime;              /* system time */
        clock_t tms_cutime;             /* user time, children */
        clock_t tms_cstime;             /* system time, children */
} Start,End;


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

#if defined __CRAY

/*  return (double)(rtclock() * 3.333e-6 / 2.); */
    return (double)( ( _rtc() / (double)CLK_TCK ) );

#else

# if ! defined __USER_TIME

    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;

# else

    double sec;
    static struct tms tmp;
    times(&tmp);
    sec = (double)(tmp.tms_utime) / (double)CLK_TCK;
    return sec;

#  endif

#endif

}


int CPTIMER ( double * user, double * system)
{
     static int cptimer_swtch = 0;
     static time_t tstart, tend;
     double sec;

     time(&tend);
     times(&End);

     /* fprintf(stderr,"\n ++++ %d ++++\n",cptimer_swtch); */
     if( ! cptimer_swtch ) {
        Start = End;
        tstart = tend;
        *user = 0.0;
        *system = 0.0;
        cptimer_swtch = 1 ;
     } else {
        cptimer_swtch = 0 ;
     }
     *user   = (double)(End.tms_utime-Start.tms_utime);
     *system = (double)(End.tms_stime-Start.tms_stime);
     *user    /= (double)CLK_TCK ;
     *system  /= (double)CLK_TCK ;
     sec = difftime( tend, tstart );
     return sec;
}


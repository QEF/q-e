/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include <time.h>
#include <sys/time.h>
#include <sys/times.h>

#include "c_defs.h"


double CCLOCK()

/* Return the second elapsed since Epoch (00:00:00 UTC, January 1, 1970)
*/

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


double SCNDS ( )

/* Return the cpu time associated to the current process 
*/

{
	static struct tms T;
	static int first = 1;
	static double init_cputime = 0.0;
	double cputime;

        times(&T);

        cputime   = (double)(T.tms_utime);
        cputime  /= (double)CLK_TCK ;

	if( first ) {
		first = 0;
		init_cputime = cputime;
	}

	return cputime - init_cputime;
}


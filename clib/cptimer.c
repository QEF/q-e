/*
  Copyright (C) 2002-2006 Quantum-Espresso group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include <time.h>
#include <sys/time.h>
#include <sys/times.h>

#include "c_defs.h"

/* This should catch systems that do not define the older CLK_TCK symbol, 
   but do define the newer CLOCKS_PER_SEC constant */

#ifndef CLK_TCK
#define CLK_TCK CLOCKS_PER_SEC
#endif

double F77_FUNC(cclock,CCLOCK)()

/* Return the second elapsed since Epoch (00:00:00 UTC, January 1, 1970)
*/

{

    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;

}

#ifdef __QK_USER__
#include <catamount/dclock.h>
double F77_FUNC(scnds,SCNDS) ( )

/* cray xt3 has no times() on the compute nodes.
 * Return the wall time associated to the current process.
*/
{
	static int first = 1;
	static double init_cputime = 0.0;
	double cputime;

        cputime   = dclock();

	if( first ) {
		first = 0;
		init_cputime = cputime;
	}

	return cputime - init_cputime;
}
#else

double F77_FUNC(scnds,SCNDS) ( )

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
#endif

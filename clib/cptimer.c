/*
  Copyright (C) 2002-2006 Quantum ESPRESSO group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#if defined(_WIN32)
#include <windows.h>
#include <sys/time.h>
#include <stdint.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif
#include <unistd.h>

#include "c_defs.h"

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

double F77_FUNC(scnds,SCNDS) ( )

/* Return the cpu time associated to the current process 
*/

{
    double sec=0.0;

#if defined(_WIN32)
    // from MSDN docs.
    FILETIME ct,et,kt,ut;
    union { FILETIME ft; uint64_t ui; } cpu;
    if (GetProcessTimes(GetCurrentProcess(),&ct,&et,&kt,&ut)) {
        cpu.ft = ut;
        sec = cpu.ui * 0.0000001;
    }
#else
    static struct rusage T;

    getrusage(RUSAGE_SELF, &T);

    sec = ((double)T.ru_utime.tv_sec + ((double)T.ru_utime.tv_usec)/1000000.0);
#endif
    return sec;
}


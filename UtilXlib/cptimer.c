/*
  Copyright (C) 2002-2006 Quantum ESPRESSO group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#if defined(_WIN32)
#include <windows.h>
#include <stdint.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

#if defined(_WIN32)

int qe_gettimeofday(struct timeval * tp, struct timezone * tzp)
{
    // Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
    // This magic number is the number of 100 nanosecond intervals since January 1, 1601 (UTC)
    // until 00:00:00 January 1, 1970 
    static const uint64_t EPOCH = ((uint64_t) 116444736000000000ULL);

    SYSTEMTIME  system_time;
    FILETIME    file_time;
    uint64_t    time;

    GetSystemTime( &system_time );
    SystemTimeToFileTime( &system_time, &file_time );
    time =  ((uint64_t)file_time.dwLowDateTime )      ;
    time += ((uint64_t)file_time.dwHighDateTime) << 32;

    tp->tv_sec  = (long) ((time - EPOCH) / 10000000L);
    tp->tv_usec = (long) (system_time.wMilliseconds * 1000);
    return 0;
}

#endif


double cclock()

/* Return the second elapsed since Epoch (00:00:00 UTC, January 1, 1970)
*/

{

    struct timeval tmp;
    double sec;
#if defined(_WIN32)
	qe_gettimeofday( &tmp, (struct timezone *)0 );
#else
    gettimeofday( &tmp, (struct timezone *)0 );
#endif
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;

}

double scnds ( )

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


#ifdef __ICC
/* Workaround for icc and incompatible includes (glibc 2.3) */
#include <time.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include "cp.h"

int C_MKDIR( char * dirname )
{
   mode_t mode = 0777 ;
   return mkdir( dirname, mode );
}

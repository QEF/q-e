#include <sys/stat.h>
#include <sys/types.h>
#include "cp.h"

int C_MKDIR( char * dirname )
{
   mode_t mode = 0777 ;
   return mkdir( dirname, mode );
}

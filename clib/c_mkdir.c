#include <sys/stat.h>
#include <sys/types.h>

int c_mkdir( char * dirname )
{
   mode_t mode = 0777 ;
   return mkdir( dirname, mode );
}

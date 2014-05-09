/*
  Copyright (C) 2003-2014 Quantum ESPRESSO group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <time.h>
#include "c_defs.h"
#include <unistd.h>

#if defined(_WIN32)
#include <direct.h>
#endif

int check_writable_dir(const char *filename) {

    struct stat sb;
    if (stat(filename, &sb) == -1) {
      return -3; /* does not exist */
      /* note: this happens also if looking for "dir/" when there is a file called "dir" */
    }
    if ( (sb.st_mode & S_IFMT) != S_IFDIR) {
      fprintf( stderr , "\ncheck_writable_dir fail: file '%s' exists but is NOT a directory\n", filename ) ;
      return -2; /* not a directory */
    }
    /* if ( ! (sb.st_mode & S_IWUSR) )        return -4; */ /* not writeble by owner */
    /* return 0 if I can read, write and execute (enter) this directory, -1 otherwise
       note: we do not actually need R_OK in Quantum-ESPRESSO;
             W_OK is definitely needed, about X_OK I'm not sure */

#if !defined(_WIN32)
    if ( access(filename, W_OK|R_OK|X_OK ) ) {
      fprintf( stderr , "\ncheck_writable_dir fail: insufficient permissions to access '%s'\n", filename ) ;
      return -1; /* no permissions  */
    }
#endif

    return 0;
}  /* check_writable_dir */

int c_mkdir_safe( const char * dirname )
{
   int retval = -1 ;

   /* return directly -1 if directory exists and is writable */
   if ( check_writable_dir(dirname) == 0) return -1;

#if defined(_WIN32)
   retval = _mkdir( dirname ) ;
#else
   mode_t mode = 0777 ;
   retval = mkdir( dirname , mode ) ;
#endif
   if ( retval == -1  && errno != EEXIST ) {
     fprintf( stderr , "\nmkdir fail: [%d] %s\n" , errno , strerror( errno ) ) ;
     retval = 1 ;
     }
   /* double check that the directory is a directory and has the good permissions */
   if ( check_writable_dir(dirname) < 0) retval = 1;  
   return retval ;
}

/* EOF */

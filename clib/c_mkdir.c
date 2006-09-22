/*
  Copyright (C) 2003 PWSCF group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <time.h>
#include "c_defs.h"

static void fatal ( const char * msg )
{

   fprintf( stderr , "fatal: %s" , *msg ? msg : "Oops!" ) ;
   exit( -1 ) ;

} /* fatal */


static void * xcmalloc ( size_t size )
{

  register void * ptr = malloc( size ) ;

  if ( ptr == NULL )
    fatal( "c_mkdir: virtual memory exhausted" ) ;
  else
    memset( ptr , 0 , size ) ;

  return ptr ;

} /* xcmalloc */


int F77_FUNC_(c_mkdir,C_MKDIR)( const char * dirname , const int * length )
{

   int retval = -1 ;

   mode_t mode = 0777 ;

   char * ldir = ( char * ) xcmalloc( (*length) + 1 ) ;

   memcpy( ldir , dirname , *length ) ;

   ldir[*length] = '\0' ;	/* memset() in xcmalloc() already do this */

   retval = mkdir( ldir , mode ) ;

   if ( retval == -1  && errno != EEXIST )
     fprintf( stderr , "mkdir fail: [%d] %s\n" , errno , strerror( errno ) ) ;

   free( ldir ) ;

   return retval ;

} /* c_mkdir_ */


/* EOF */

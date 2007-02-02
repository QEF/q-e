/*
  Copyright (C) 2003-2007 Quantum-Espresso group
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

} /* c_mkdir */

/* call from fortran as
   ios = c_remame ( TRIM(old-file-name), TRIM_LEN(old-file-name), &
                    TRIM(new-file-name), TRIM_LEN(new-file-name) )
   renames file old-file-name into new-file-name (don't try this on
   open files!)
   ios should return 0 if everything is ok, -1 otherwise.
   Written by PG by imitating "c_mkdir" without really understanding it */

int F77_FUNC_(c_rename,C_RENAME)( const char * oldname, const int * oldlength ,
                                  const char * newname, const int * newlength )
{

   int retval = -1 ;

   char * oldname_ = ( char * ) xcmalloc( (*oldlength) + 1 ) ;
   char * newname_ = ( char * ) xcmalloc( (*newlength) + 1 ) ;

   memcpy( oldname_ , oldname , *oldlength ) ;
   memcpy( newname_ , newname , *newlength ) ;

   oldname_[*oldlength] = '\0' ;
   newname_[*newlength] = '\0' ;

   retval = rename( oldname_, newname_ ) ;

   if ( retval == -1 )
     fprintf( stderr , "mv fail: [%d] %s\n" , errno , strerror( errno ) ) ;

   free( oldname_ ) ;
   free( newname_ ) ;

   return retval ;

} /* c_rename */

/* EOF */

/*
 Copyright (C) 2005-2008 Quantum ESPRESSO group
 This file is distributed under the terms of the
 GNU General Public License. See the file `License'
 in the root directory of the present distribution,
 or http://www.gnu.org/copyleft/gpl.txt . 

------------------------------------------------------
*/


#include <stdio.h>
#include <stdlib.h> 
#include "c_defs.h"
#include "md5.h"

#define MAX_BUF 1024 

char *strCopy( char * str ) 
{ 
  if( !str ) 
  { 
    return NULL; 
  } 
   
  char *n = ( char * )malloc( strlen( str ) + 1 ); 
  memcpy( n, str, strlen( str ) + 1 ); 
   
  return n; 
} 

char *readFile( FILE *fp ) 
{ 
  if( !fp ) 
  { 
    return NULL; 
  } 

  char *output = NULL; 
  char buf[ MAX_BUF ]; 
  int bytes, byteCount = 0; 
   
  memset( buf, 0, sizeof( buf ) ); 
   
  while( ( bytes = fread( buf, sizeof( char ), sizeof( buf ) - 1, fp ) ) > 0 ) 
  { 
    if( !output ) 
    { 
      /* alloca */ 
      output = ( char * )malloc( bytes + 1 ); 
      memcpy( output, buf, bytes + 1 ); 
    } 
    else 
    { 
      /* rialloca se esiste */

      char *oldStr = strCopy( output ); 
      output = ( char * )realloc( output, byteCount + bytes + 1 ); 
      memcpy( output, oldStr, byteCount ); 
       
      free( oldStr ); 
      memcpy( output + byteCount, buf, bytes + 1 ); 
    } 
   
    byteCount += bytes; 
    memset( buf, 0, sizeof( buf ) ); 
  } 
   
  return output; 
} 


void get_md5(const char *file, char md5[32], int err)
{
 
     FILE *fp;
     char *data;
     //char THISMD5[31];
     md5_state_t state;
     md5_byte_t digest[16];

     if(file==NULL) { 
	err = 1;
        return;
     }

     fp=fopen(file,"rb");
     if(fp==NULL) {
	err = 2;
	return;
     }
     
     data=readFile(fp);
     if(data==NULL) {
	err = 3;
	return;
     }

     md5_init(&state);
     md5_append(&state,(const md5_byte_t *)data,strlen(data));
     md5_finish(&state,digest);

     int i=0;
     for(i;i<16;i++){
           //snprintf(THISMD5+i*2,sizeof(THISMD5),"%02x",digest[i]);
	   snprintf(md5+i*2,sizeof(md5),"%02x",digest[i]);
     }
     fclose(fp);

     //strcpy(md5, THISMD5);
     //md5 = &THISMD5;
     err = 0;
     return;
}

int F77_FUNC_(file_md5,FILE_MD5)( const int * f_name, const int * f_len, int* out )
{
     int i, err = -1 ;
     char md5[32];
     char *f = ( char * ) malloc( (*f_len) + 1 ) ;

     for( i = 0; i < * f_len; i++ ) f[ i ] = (char)f_name[ i ];

     f[*f_len] = '\0' ;
     
     get_md5( f ,  md5, err) ;
     
     for( i = 0; i < 32; i++ ) out[ i ] = md5[ i ]; 
    
     return err;
} 



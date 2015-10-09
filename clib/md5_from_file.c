/*
 Copyright (C) 2011-2013 Quantum ESPRESSO group
 This file is distributed under the terms of the
 GNU General Public License. See the file `License'
 in the root directory of the present distribution,
 or http://www.gnu.org/copyleft/gpl.txt . 

------------------------------------------------------
*/


#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include "c_defs.h"
#include "md5.h"

#define MAX_BUF 1024 


static void fatal ( const char * msg )
{
  
   fprintf( stderr , "fatal: %s" , *msg ? msg : "Oops!" ) ;
   exit( -1 ) ;

} /* fatal */

static void * xcmalloc ( size_t size )
{

  register void * ptr = malloc( size ) ;

  if ( ptr == NULL )
    fatal( "md5_from_file: virtual memory exhausted" ) ;
  else
    memset( ptr , 0 , size ) ;

  return ptr ;

} /* xcmalloc */

char *readFile( FILE *file ) 
{ 

	char *out;
	unsigned long fileLen;

	if (!file)
	{
		exit(1);
	}

	fseek(file, 0, SEEK_END);
	fileLen=ftell(file);
	fseek(file, 0, SEEK_SET);

	out=(char *)xcmalloc(fileLen+1);

	if (!out)
	{
	fprintf(stderr, "Memory error!");
	fclose(file);
	exit(1);
	}

	fread(out, fileLen, 1, file);

	return out; 

} 


void get_md5(const char *file, char *md5, int err)
{
 
     FILE *fp;
     char *data;
     int i;
     md5_state_t state;
     md5_byte_t digest[16];
     
     if(file==NULL) { 
        strcpy(md5,"Not computed, file name is empty\n");
	err = 1;
        return;
     }

     fp=fopen(file,"rb");
     if(fp==NULL) {
        strcpy(md5,"Not computed, couldn't open file\n");
	err = 2;
	return;
     }

     data=readFile(fp);
     if(data==NULL) {
        strcpy(md5,"Not computed, file is empty     \n");
	err = 3;
	return;
     }

     md5_init(&state);
     md5_append(&state,(const md5_byte_t *)data,strlen(data));
     md5_finish(&state,digest);

     for(i=0;i<16;i++){
	   snprintf(md5+i*2,sizeof(md5),"%02x",digest[i]);
     }
     fclose(fp);

     free(data);
     err = 0;
     return;
}


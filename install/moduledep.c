/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<strings.h>

#define MAXFILES 5000
#define NAMELEN 80
#define MODLEN 32
#define MAXMOD 500
#define DEBUG 0

char str[BUFSIZ];

struct module_list {
  char filename[NAMELEN];
  int nmodules;
  char modules[MAXMOD][MODLEN];
};

int used_name(char *str, char *module)
{
  int i=0,j=0,k;

  /* skip blanks */
  while( isspace( str[i] ) ) i++;                         

  /* convert to lowercase */
  for( k = i; k < i+4; k++ ) str[k] = tolower( str[k] );  

  /* check if we found a module */
  if( ! strncmp( &str[i], "use ", 4 ) ) {                 
    i += 4;                                               /* skip "use " */
    while( isspace(str[i]) ) i++;                         /* skip blanks */
    while( str[i]!=' ' && str[i]!= ',' && str[i]!='\n' && j<(MODLEN-1) ) {
       module[j++] = tolower(str[i++]);                   /* get module name */
    }
    module[j] = '\0';
    return 1;
  }
  return 0;
}
    
int module_unit(char *str, char *module)
{
  int i=0,j=0,k;
  while(isspace(str[i])) i++;
  for(k=i;k<i+7;k++) str[k] = tolower(str[k]);
  if(!strncmp(&str[i],"module ",7)) {
    i = i+7;
    while(isspace(str[i])) i++;
    while( str[i]!=' ' && str[i]!='\n' && j<(MODLEN-1) ) module[j++] = tolower(str[i++]);
    module[j] = '\0';
    return 1;
  }
  return 0;
}
    

int used_modules(char *filename, char modules[][MODLEN])
{
   FILE *fp;
   char module_name[MODLEN];
   int i=0, nused = 0, ok;
   fp=fopen(filename,"r");
   while(fgets(str,BUFSIZ-1,fp)) {
     if( used_name( str, module_name ) ) {
       ok = 1;
       for( i = 0; i < nused; i ++ ) {
         if( ! strcmp( module_name, modules[ i ] ) ) {
           ok = 0;
           break;
         }
       }
       if( ok ) {
         strcpy( modules[ nused ], module_name );
         nused ++;
       }
       if( nused == MAXMOD ) {
         fprintf(stderr,"\n maximum number of modules exceded\n");
         exit(1);
       }
     }
   }
   fclose(fp);
   return nused;
}

int defined_modules(char *filename, char modules[][MODLEN])
{
   FILE *fp;
   char module_name[MODLEN];
   int i=0, ndef = 0, ok;

   fp = fopen(filename,"r");

   while(fgets(str,BUFSIZ-1,fp)) {
     if( module_unit( str, module_name ) ) {
       ok = 1;
       for( i = 0; i < ndef; i ++ ) {
         if( ! strcmp( module_name, modules[ i ] ) ) {
           ok = 0;
           break;
         }
       }
       if( ok ) {
         strcpy( modules[ ndef ], module_name );
         ndef ++;
       }
       if( ndef == MAXMOD ) {
         fprintf(stderr,"\n maximum number of modules exceded\n");
         exit(1);
       }
     }
   }
   fclose(fp);
   return ndef;
}


int change_extension(char sour[], char dest[], char new_extension[])
{
  int i,len;
  len = (int)strlen(sour);
  for(i=0; i<len; i++ ) { dest[i] = sour[i]; } 
  for(i=len-1; i>0 && dest[i]!= '.' ; i-- ) ;
  dest[i] = '\0';
  strcat(dest,new_extension);
  return strlen(dest);
}

int main(int argc,char *argv[])
{
  int i, nm, j, k, len, lobj, ntarget, ndepend;
  char object[BUFSIZ];
  FILE * target_list;
  FILE * depend_list;
  char targets[MAXFILES][NAMELEN]; /* list of targets file Â*/
  char depends[MAXFILES][NAMELEN]; /* list of files upon which the targets may depend Â*/
  struct module_list * use;
  struct module_list * module;


  target_list = fopen( argv[1], "r" );
  depend_list = fopen( argv[2], "r" );

  i = 0;
  while( ( fscanf( target_list, "%s", targets[i] ) == 1 ) && ( i < MAXFILES ) ) i++ ;
  ntarget = i;
  use = ( struct module_list * ) malloc( sizeof( struct module_list ) * ntarget );

  i = 0;
  while( ( fscanf( depend_list, "%s", depends[i] ) == 1 ) && ( i < MAXFILES ) ) i++ ;
  ndepend = i;
  module = ( struct module_list * ) malloc( sizeof( struct module_list ) * ndepend );
 
  for( i = 0; i < ntarget; i++ ) {
     change_extension( targets[i], use[i].filename, ".o");
     use[i].nmodules =  used_modules( targets[i], use[i].modules );
     if( DEBUG ) printf("%s %d\n", targets[i], use[i].nmodules );
  }
  if( DEBUG ) printf("ntarget = %d",ntarget);

  for( i = 0; i < ndepend; i++ ) {
     change_extension( depends[i], module[i].filename, ".o");
     module[i].nmodules =  defined_modules( depends[i], module[i].modules );
     if( DEBUG ) printf("%s %d\n", depends[i], module[i].nmodules );
  }
  if( DEBUG ) printf("ndepend = %d",ndepend);

  printf("\n");

  for( i = 0; i < ntarget; i++ ) {
    printf("\n%s : ", use[i].filename );
    len = strlen( use[i].filename ) + 3;
    for( nm = 0; nm < use[i].nmodules; nm ++ ) {
      for( j = 0; j < ndepend; j++ ) {
        for( k = 0; k < module[j].nmodules; k ++ ) {
          if( DEBUG ) printf("%s %s\n",use[i].modules[nm], module[j].modules[k]);
          if( ! strcmp( use[i].modules[nm], module[j].modules[k] ) ) {
            lobj = strlen( module[j].filename );
            if( (len+lobj) > 70 ) { printf(" \\\n    "); len = 4; }
            printf("%s ", module[j].filename );
            len += lobj;
          }
        }
      }
    }
    printf("\n");
  }
 

  fclose( target_list );
  fclose( depend_list );
}

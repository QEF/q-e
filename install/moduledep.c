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

#define MMAX 10000

char str[BUFSIZ];
char module[MMAX][BUFSIZ];

int used_name(char *str, char *module)
{
  int i=0,j=0,k;
  while(isspace(str[i])) i++;
  for(k=i;k<i+4;k++) str[k] = tolower(str[k]);
  if(!strncmp(&str[i],"use ",4)) {
    i = i+4;
    while(isspace(str[i])) i++;
    while( str[i]!=' ' && str[i]!= ',' && str[i]!='\n') {
       module[j++] = tolower(str[i++]);
    }
    module[j] = '\0';
    return 1;
  }
  return 0;
}
    
int module_name(char *str, char *module)
{
  int i=0,j=0,k;
  while(isspace(str[i])) i++;
  for(k=i;k<i+7;k++) str[k] = tolower(str[k]);
  if(!strncmp(&str[i],"module ",7)) {
    i = i+7;
    while(isspace(str[i])) i++;
    while( str[i]!=' ' && str[i]!='\n') module[j++] = tolower(str[i++]);
    module[j] = '\0';
    return 1;
  }
  return 0;
}
    

int subroutine_modules(char *filename, char modules[][BUFSIZ])
{
   FILE *fp;
   int i=0;
   fp=fopen(filename,"r");
   while(fgets(str,BUFSIZ-1,fp)) {
     if( used_name(str,modules[i]) ) i++;
     if( i == MMAX ) {
       fprintf(stderr,"\n maximum number of modules exceeded\n");
       exit(1);
     }
   }
   fclose(fp);
   return i;
}

int module_file(char *filename, char modules[][BUFSIZ], int nm)
{
   FILE *fp;
   char tmp_module[BUFSIZ];
   int i=0,ok=0;
   fp=fopen(filename,"r");
   while(fgets(str,BUFSIZ-1,fp)) {
     if( module_name(str,tmp_module) ) { ok = 1; break; }
   }
   if(ok) {
     ok = 0;
     for(i=0;i<nm;i++) {
        if(!strcmp(modules[i],tmp_module) ) { ok = 1; break; }
     }
   }
   fclose(fp);
   return ok;
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
  int i,nm,len,lobj;
  char object[BUFSIZ];

  nm =  subroutine_modules(argv[1],module);
 
  change_extension(argv[1],object,".o");
  printf("\n%s : ",object);
  len = strlen(object)+3;

  for(i=2;i<argc;i++) {
    if( strcmp(argv[1],argv[i]) ) {
       if(module_file(argv[i],module,nm)) {
          lobj = change_extension(argv[i],object,".o");
          if( (len+lobj) > 70 ) {
            printf(" \\\n    ");
            len = 4;
          }
          printf("%s ",object);
          len += lobj;
       }
    }
  }
 
  printf("\n");
}

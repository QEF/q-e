/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include"cp.h"


void cc_getenv(char *name,char *value)
{
  char *ptr;

  ptr = getenv(name);
  if(ptr) {
    while ( *ptr ) *value++ = *ptr++ ;
  }
  *value = '\0';

}

int CP_GETENV(int *out_nruns,char infile[][256],
		char outfile[][256],char cdir[],char input_dir[])
{
  char runs[10][3];
  char value[256];
  char filename[256];
  char output_dir[256];
  int  i,j,k;
  int  nruns;


  cc_getenv("PWD",cdir);
  cc_getenv("CP_INPUT_DIR",input_dir);
/*  printf("%s\n",input_dir);*/ 
  cc_getenv("CP_OUTPUT_DIR",output_dir);
/*  printf("%s\n",output_dir);*/ 
  cc_getenv("CP_FILENAME",filename);
/*  printf("%s\n",filename); */  
  cc_getenv("CP_RUNS",value);
/*  printf("%s\n",value); */ 

  for(k=0,i=0; k<10 && value[i] ;k++) {
    for(j=0; j<2; j++, i++) runs[k][j] = value[i];
    runs[k][j] = '\0';
    i++;
  }
  nruns = k;

  for(k=0;k<nruns;k++) {

     printf("Run : %s\n",runs[k]);

     strcpy(value,input_dir);
     strcat(value,"/");
     strcat(value,filename);
     strcat(value,".in.");
     strcat(value,runs[k]);

     /* printf("input  file : %s\n",value); */
     for(i=0; value[i]; i++) infile[k][i] = value[i];
     infile[k][i]=' ';

     strcpy(value,output_dir);
     strcat(value,"/");
     strcat(value,filename);
     strcat(value,".out.");
     strcat(value,runs[k]);

     /*  printf("output file : %s\n",value); */
     for(i=0; value[i]; i++) outfile[k][i] = value[i];
     outfile[k][i]=' ';


  }

  for(i=0; cdir[i] && i<255 ; i++) ;
  if(i>0) cdir[i++] = '/';
  cdir[i] = ' ';
  for(i=0; input_dir[i] && i<255 ; i++) ;
  if(i>0) input_dir[i++] = '/';
  input_dir[i] = ' ';

  *out_nruns = nruns;

  return 0;

}



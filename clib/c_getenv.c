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


void use_c_getenv( char *name, char *value, int nchmax )
{
  char *ptr;
  int  ic;

  ptr = getenv( name );
  ic  = 0;
  if( ptr ) {
    while ( *ptr && ic < (nchmax-1) ) { *value++ = *ptr++ ; ic++; }
  }
  *value = '\0';

}

int CP_GETENV( /* WAITING FOR ARGUMENTS */ )
{
  char pwd[256];
  int  nchmax = 256;

  use_c_getenv( "PWD", pwd, nchmax );

  return 0;
}

/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include<stdio.h>
#include"cp.h"

static char S[BUFSIZ];

void READOCC(int *nx,double *fi)
{
 int i;
 for(i=0;i< (*nx);i++) fscanf(stdin,"%lf",fi+i);
 gets(S);
}

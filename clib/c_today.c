/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include<time.h>
#include<stdio.h>
#include<strings.h>
#include"cp.h"

int CP_DATE(char * today, int * len)
{
 time_t tempo;
 char * pt;
 int i;
 time(&tempo);
 pt = ctime(&tempo);
 for(i=0; i<(*len) && pt[i]!='\0'; i++) today[i]=pt[i];
 today[i] = '\0';
 return strlen(today);
}

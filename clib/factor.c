/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include<stdio.h>
#include"cp.h"

/* this routine is used in c3fft (t3e version).
   It determines if a given number "num" is factorizable in power of
   2, 3, 5 ( output isys = 0 ) or not ( output isys = 1 ).
*/

/* Carlo Cavazzoni , 12 feb. 97 */ 

void FACTOR235(int *num,int *isys)
{
 int i;
 i = *num;
 while((i%2)==0) i/=2; 
 while((i%3)==0) i/=3; 
 while((i%5)==0) i/=5; 
 
 if(i==1) *isys = 0;
 else *isys = 1;
}

void FACTOR2(int *num,int *isys)
{
 int i;
 i = *num;
 while((i%2)==0) i/=2; 
 if(i==1) *isys = 0;
 else *isys = 1;
}

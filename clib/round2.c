/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include"cp.h"

/* this routine given a number "n" find the smallest power of
   2 "i" greater or equal to "n"
*/

/* Carlo Cavazzoni , 2 dec. 98 */

void ROUND2(int *n,int *i)
{
   if((*n)<=1) {
      (*i) = 1;
   } else {
      register int j = (*n)-1;
      (*i) = 2;
      while(j>>=1) (*i)<<=1;
   }
}

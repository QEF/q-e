/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include<stdio.h>
#include<stdlib.h>

/* qsort - quick sort

   qsort(n,comp,swap)
   unsigned n;
   int (*comp)();
   int (*swap)();
                                  ***** see bsort for parameters

*/

static unsigned _rearr(unsigned lb,unsigned ub);
static void _quick(unsigned lb,unsigned ub);
static int (*_comp)(unsigned,unsigned), (*_swap)(unsigned,unsigned);

void Qsort(unsigned n,int (*comp)(),int (*swap)())
{
        _comp = comp;
        _swap = swap;
        _quick(0,n-1);
}


static void _quick(unsigned lb,unsigned ub)
{
unsigned j;

if(lb<ub)
        {
        if((j = _rearr(lb,ub)))
        _quick(lb,j-1);
        _quick(j+1,ub);
        }
}


static unsigned _rearr(unsigned lb,unsigned ub)
{
do
        {
        while(ub > lb && (*_comp)(ub,lb) >=0)   ub--;

        if(ub != lb)
                {
                (*_swap)(ub,lb);
                while(lb<ub && (*_comp)(lb,ub)<=0) lb++;
                if(lb != ub)    (*_swap)(lb,ub);
                }
        } while(lb != ub);

return lb;
}


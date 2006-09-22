/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include <stdio.h>
#include <stdlib.h>

#include "c_defs.h"

#define MAX_INDEX 32768

struct Index { unsigned char i[8]; } ;

static struct Index * P_Index;
static int * P_IndexIndex;

static struct Index * LN;
static int          * IG;
static int            LN_SIZE;

int IndexCmp( struct Index * A, struct Index * B)
{
    int i;

    for(i = 7; i>=0 ; i--) {
        if(A->i[i] > B->i[i] ) {
            return +1;
        }
        else if(A->i[i] < B->i[i]) { 
            return -1;
        }
    }

    return 0;
}
    

int index_comp(unsigned i,unsigned j)
{
int cmp;
cmp = IndexCmp(P_Index + i, P_Index + j);
if      ( cmp > 0  ) return 1;
else if ( cmp == 0 ) return 0;
return -1;
}

int index_swap(unsigned i,unsigned j)
{
static struct Index tmp;
static int itmp;

tmp         = P_Index[j] ;
P_Index[j]  = P_Index[i] ;
P_Index[i]  = tmp        ;

itmp             = P_IndexIndex[j] ;
P_IndexIndex[j]  = P_IndexIndex[i] ;
P_IndexIndex[i]  = itmp        ;

return 1;
}


int IndexSort(struct Index * A, int * IndexIndex, int n)
{
   void Qsort(unsigned n,int (*comp)(),int (*swap)());
   P_Index = A;
   P_IndexIndex = IndexIndex;
   Qsort((unsigned)n,index_comp,index_swap);
   return 1;
}


int IndexSet( struct Index * A, int I1, int I2, int I3 )
{
    unsigned int himask = 0xFF00;
    unsigned int lomask = 0x00FF;

    if(abs(I1)>=MAX_INDEX || abs(I2)>=MAX_INDEX || abs(I3)>=MAX_INDEX ) {
      return -1;
    }

    if(I1<0) I1 += MAX_INDEX;
    if(I2<0) I2 += MAX_INDEX;
    if(I3<0) I3 += MAX_INDEX;


    A->i[7] = (unsigned char ) 0; 
    A->i[6] = (unsigned char ) 0; 
    A->i[5] = (unsigned char ) ((himask & (unsigned int) I1)>>8); 
    A->i[4] = (unsigned char ) ( lomask & (unsigned int) I1); 
    A->i[3] = (unsigned char ) ((himask & (unsigned int) I2)>>8); 
    A->i[2] = (unsigned char ) ( lomask & (unsigned int) I2); 
    A->i[1] = (unsigned char ) ((himask & (unsigned int) I3)>>8); 
    A->i[0] = (unsigned char ) ( lomask & (unsigned int) I3); 
    return 0;
}

int IndexShow(struct Index A)
{
    int i;
    for(i=7;i>=0;i--) printf("%2x",A.i[i]);
    printf("\n");
    return 0;
}

int IndexFind(struct Index * A, int n, struct Index * B)
{
    int lb, ub, i, cmp;

    lb = 0;
    ub = n-1;
    i  = lb;

    while(lb<(ub-1)) {
       i = lb + (ub - lb)/2;
       cmp = IndexCmp(B,&A[i]);
       if(cmp>0) {
          lb = i;
       } else if(cmp<0) { 
          ub = i;
       } else {
          ub = lb = i;
       }
    } 
    if(lb<ub) {
       cmp = IndexCmp(B,&A[lb]);
       if(cmp) {
          i = ub;
       } else {
          i = lb;
       }
    }

    if ( IndexCmp(B,&A[i]) ) return -1;

    return i;
} 

void F77_FUNC_(ln_alloc,LN_ALLOC)(int * LN_DIM)
{
    LN_SIZE = * LN_DIM;
    LN = ( struct Index *) malloc ( LN_SIZE * sizeof( struct Index ));
    IG = ( int          *) malloc ( LN_SIZE * sizeof( int          ));
}

void F77_FUNC_(ln_dealloc,LN_DEALLOC)(void )
{
    free((void *)LN);
    free((void *)IG);
}

void F77_FUNC_(ln_set,LN_SET)(int * IRI1, int * IRI2, int * IRI3, int * ig)
{
    if( *ig<1 || *ig > LN_SIZE) {
       exit(*ig);
    }
    IndexSet( &LN[*ig-1], *IRI1, *IRI2, *IRI3 );
    IG[*ig-1] = *ig;

}

int F77_FUNC_(ln_activate,LN_ACTIVATE)()
{
    IndexSort(LN,IG,LN_SIZE);
    return 0;   
}

int F77_FUNC_(ln_ind,LN_IND)(int * IRI1, int * IRI2, int * IRI3)
{
    static struct Index B;
    static int          ib;

    IndexSet(&B,*IRI1,*IRI2,*IRI3);
    ib = IndexFind(LN,LN_SIZE,&B);
    if(ib>=0) return IG[ib];
    return -1;
}

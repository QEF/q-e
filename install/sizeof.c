#include<stdio.h>
#include<stdlib.h>
#if defined __FFTW 
#include"../clib/fftw.h"
#endif


int main()
{
#if defined __FFTW
  printf("SIZEOF ( fftw_plan * ). . = %d\n",(int)sizeof(fftw_plan *));
  printf("SIZEOF ( FFTW_COMPLEX * ) = %d\n",(int)sizeof(FFTW_COMPLEX *));
  printf("SIZEOF ( FFTW_COMPLEX ) . = %d\n",(int)sizeof(FFTW_COMPLEX));
#endif
  printf("SIZEOF ( void * ) . . . . = %d\n",(int)sizeof(void *));
  printf("SIZEOF ( double ) . . . . = %d\n",(int)sizeof(double));
  printf("SIZEOF ( float ). . . . . = %d\n",(int)sizeof(float));
  printf("SIZEOF ( int ). . . . . . = %d\n",(int)sizeof(int));
  printf("SIZEOF ( char ) . . . . . = %d\n",(int)sizeof(char));
}

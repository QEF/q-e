#include<stdio.h>
#include<stdlib.h>
#if defined __FFTWDRV
#include"fftw.h"
#endif


int main()
{
#if defined __FFTWDRV | defined __FFTW
  printf("SIZEOF ( fftw_plan * ). . = %d\n",sizeof(fftw_plan *));
  printf("SIZEOF ( FFTW_COMPLEX * ) = %d\n",sizeof(FFTW_COMPLEX *));
  printf("SIZEOF ( FFTW_COMPLEX ) . = %d\n",sizeof(FFTW_COMPLEX));
#endif
  printf("SIZEOF ( void * ) . . . . = %d\n",sizeof(void *));
  printf("SIZEOF ( double ) . . . . = %d\n",sizeof(double));
  printf("SIZEOF ( float ). . . . . = %d\n",sizeof(float));
  printf("SIZEOF ( int ). . . . . . = %d\n",sizeof(int));
  printf("SIZEOF ( char ) . . . . . = %d\n",sizeof(char));
}

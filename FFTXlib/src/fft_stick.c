/*
 * Copyright (C) Quantum ESPRESSO group
 *
 * This file is distributed under the terms of the
 * GNU General Public License. See the file `License'
 * in the root directory of the present distribution,
 * or http://www.gnu.org/copyleft/gpl.txt .
 */

/* Always compile FFTW beside external driver to be used as fallback
 * hook for higher level FFT drivers in the FFTX library */

#include "fftw_dp.h"
#include "fftw_sp.h"

int create_plan_1d (fftw_plan *p, int *n, int *idir)
{
   fftw_direction dir = ( (*idir < 0) ? FFTW_FORWARD : FFTW_BACKWARD ); 
   *p = qe_fftw_create_plan(*n, dir, FFTW_ESTIMATE | FFTW_IN_PLACE);
   if( *p == NULL ) fprintf(stderr," *** CREATE_PLAN: warning empty plan ***\n");
/*   printf(" pointer size = %d, value = %d\n", sizeof ( *p ), *p ); */
   return 0;
}

int float_create_plan_1d(float_fftw_plan* p, int* n, int* idir)
{
	fftw_direction dir = ((*idir < 0) ? FFTW_FORWARD : FFTW_BACKWARD);
	*p = qe_float_fftw_create_plan(*n, dir, FFTW_ESTIMATE | FFTW_IN_PLACE);
	if (*p == NULL) fprintf(stderr, " *** CREATE_PLAN: warning empty plan ***\n");
	/*   printf(" pointer size = %d, value = %d\n", sizeof ( *p ), *p ); */
	return 0;
}

int destroy_plan_1d (fftw_plan *p)
{
   if ( *p != NULL ) qe_fftw_destroy_plan(*p);
   else fprintf(stderr," *** DESTROY_PLAN: warning empty plan ***\n");
   return 0;
}

int float_destroy_plan_1d(float_fftw_plan* p)
{
	if (*p != NULL) qe_float_fftw_destroy_plan(*p);
	else fprintf(stderr, " *** DESTROY_PLAN: warning empty plan ***\n");
	return 0;
}

int create_plan_2d (fftwnd_plan *p, int *n, int *m, int *idir)
{
   fftw_direction dir = ( (*idir < 0) ? FFTW_FORWARD : FFTW_BACKWARD );
   *p = qe_fftw2d_create_plan(*m, *n, dir, FFTW_ESTIMATE | FFTW_IN_PLACE);
   if( *p == NULL ) fprintf(stderr," *** CREATE_PLAN_2D: warning empty plan ***\n");
/*   printf(" pointer size = %d, value = %d\n", sizeof ( *p ), *p ); */
   return 0;
}

int destroy_plan_2d (fftwnd_plan *p)
{
   if ( *p != NULL ) qe_fftwnd_destroy_plan(*p);
   else fprintf(stderr," *** DESTROY_PLAN_2D: warning empty plan ***\n");
   return 0;
}

int create_plan_3d (fftwnd_plan *p, int *n, int *m, int *l, int *idir)
{
   fftw_direction dir = ( (*idir < 0) ? FFTW_FORWARD : FFTW_BACKWARD );
   *p = qe_fftw3d_create_plan(*l, *m, *n, dir, FFTW_ESTIMATE | FFTW_IN_PLACE);
   if( *p == NULL ) {
	fprintf(stderr," *** CREATE_PLAN_3D: warning empty plan ***\n");
	fprintf(stderr," *** input was (n,m,l,dir): %d %d %d %d ***\n", *l, *m, *n, *idir);
   }
/*   printf(" pointer size = %d, value = %d\n", sizeof ( *p ), *p ); */
   return 0;
}

int destroy_plan_3d (fftwnd_plan *p)

{
   if ( *p != NULL ) qe_fftwnd_destroy_plan(*p);
   else fprintf(stderr," *** DESTROY_PLAN_3D: warning empty plan ***\n");
   return 0;
}


int fft_x_stick
(fftw_plan *p, FFTW_COMPLEX *a, int *nx, int *ny, int *nz, int *ldx, int *ldy )
{

   int i;
   int xstride, bigstride;
   int xhowmany, xidist;

/* trasform  along x and y */
   bigstride = (*ldx) * (*ldy);

   xhowmany = (*ny);
   xstride  = 1;
   xidist   = (*ldx);

   /* ptr = (double *)a; */

   for(i = 0; i < *nz ; i++) {
     /* trasform  along x */
     fftw(*p,xhowmany,&a[i*bigstride],xstride,xidist,0,0,0);
   }
   return 0;
}

int fft_y_stick
   (fftw_plan *p, FFTW_COMPLEX *a, int *ny, int *ldx )
{
   fftw(*p, 1, a, (*ldx), 1, 0, 0, 0);
   return 0;
}


int fft_z_stick
   (fftw_plan *p, FFTW_COMPLEX *zstick, int *ldz, int *nstick_l)
{
   int howmany, idist;
   howmany = (*nstick_l) ;
   idist   = (*ldz);
   fftw(*p, howmany, zstick, 1, idist, 0, 0, 0);
   return 0;
}

int fftw_inplace_drv_1d
   (fftw_plan *p, int *nfft, FFTW_COMPLEX *a, int *inca, int *idist)
{
   fftw(*p, (*nfft), a, (*inca), (*idist), 0, 0, 0);
   return 0;
}

int fftw_inplace_drv_2d
  ( fftwnd_plan *p, int *nfft, FFTW_COMPLEX *a, int *inca, int *idist)
{
   fftwnd( *p, (*nfft), a, (*inca), (*idist), 0, 0, 0 );
   return 0;
}

int fftw_inplace_drv_3d
   ( fftwnd_plan *p, int *nfft, FFTW_COMPLEX *a, int *inca, int *idist)
{
   fftwnd( *p, (*nfft), a, (*inca), (*idist), 0, 0, 0 );
   return 0;
}

int fft_x_stick_single
(fftw_plan *p, FFTW_COMPLEX *a, int *nx, int *ny, int *nz, int *ldx, int *ldy )
{

   int xstride;
   int xhowmany, xidist;

/* trasform  along x and y */
   xhowmany = (*ny);
   xstride  = 1;
   xidist   = (*ldx);

   fftw(*p,xhowmany,a,xstride,xidist,0,0,0);

   return 0;
}


int fft_z_stick_single (fftw_plan *p, FFTW_COMPLEX *a, int *ldz)
{
  fftw(*p, 1,a, 1, 0, 0, 0, 0);

  return 0;
}

int float_fft_z_stick_single(float_fftw_plan* p, FFTW_FLOAT_COMPLEX* a, int* ldz)
{
	float_fftw(*p, 1, a, 1, 0, 0, 0, 0);

	return 0;
}

/* Computing the N-Dimensional FFT 
void fftwnd(fftwnd_plan plan, int howmany,
            FFTW_COMPLEX *in, int istride, int idist,
            FFTW_COMPLEX *out, int ostride, int odist);
*/


/*
void fftw(fftw_plan plan, int howmany,
          fftw_complex *in, int istride, int idist,
          fftw_complex *out, int ostride, int odist);

The function fftw computes the one-dimensional Fourier transform, 
using a plan created by fftw_create_plan (See Section Plan Creation for 
One-dimensional Transforms.) The function fftw_one provides a simplified 
interface for the common case of single input array of stride 1. 

Arguments
plan    is the plan created by fftw_create_plan 
howmany is the number of transforms fftw will compute. It is faster to 
   tell FFTW to compute many transforms, instead of simply calling fftw 
   many times. 
in, istride and idist describe the input array(s). 
  There are howmany input arrays; the first one is pointed to by in, 
  the second one is pointed to by in + idist, and so on, up to 
  in + (howmany - 1) * idist. Each input array consists of complex numbers, 
  which are not necessarily contiguous in memory. Specifically, in[0] is 
  the first element of the first array, in[istride] is the second element 
  of the first array, and so on. In general, the i-th element of the j-th 
  input array will be in position in[i * istride + j * idist]. 
out, ostride and odist describe the output array(s). 
  The format is the same as for the input array. 

In-place transforms: If the plan specifies an in-place transform, ostride 
and odist are always ignored. If out is NULL, out is ignored, too. Otherwise, 
out is interpreted as a pointer to an array of n complex numbers, that FFTW 
will use as temporary space to perform the in-place computation. out is used 
as scratch space and its contents destroyed. In this case, out must be an 
ordinary array whose elements are contiguous in memory (no striding). 

*/

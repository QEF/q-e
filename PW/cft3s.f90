!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

# if defined __AIX || __FFTW || __SGI
#  define __FFT_MODULE_DRV
# endif


#ifdef __PARA

!
!----------------------------------------------------------------------
subroutine cft3s (f, n1, n2, n3, nx1, nx2, nx3, sign)
  !----------------------------------------------------------------------
  !
  !   sign = +-1 : parallel 3d fft for rho and for the potential
  !   sign = +-2 : parallel 3d fft for wavefunctions
  !
  !   sign = + : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !              fft along z using pencils        (cft_1)
  !              transpose across nodes           (fft_scatter)
  !                 and reorder
  !              fft along y (using planes) and x (cft_2)
  !   sign = - : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  !              fft along x and y(using planes)  (cft_2)
  !              transpose across nodes           (fft_scatter)
  !                 and reorder
  !              fft along z using pencils        (cft_1)
  !
  !   The array "planes" signals whether a fft is needed along y :
  !     planes(i)=0 : column f(i,*,*) empty , don't do fft along y
  !     planes(i)=1 : column f(i,*,*) filled, fft along y needed
  !   "empty" = no active components are present in f(i,*,*)
  !             after (sign>0) or before (sign<0) the fft on z direction
  !
  !   Note that if sign=+/-1 (fft on rho and pot.) all fft's are needed
  !   and all planes(i) are set to 1
  !
#include "machine.h"

#if defined __FFT_MODULE_DRV
  use fft_scalar, only: cft_1z, cft_2xy
#endif

  use fft_base, only: fft_scatter

  use parameters, only: DP
  use para, only: ncts, ncplanes, ncp0s, nkcp, nprocp, nxxs, npps, ncps, me
  use sticks, only: dffts

  implicit none
  integer :: n1, n2, n3, nx1, nx2, nx3, sign
  complex (kind=DP) :: f ( nxxs )
  !
  integer :: mc, i, j, ii, iproc, k, nppx
  complex (kind=DP), allocatable :: aux (:)
  integer :: planes ( nx1 )

#if defined(__FFTW)
# define CFT_1S cft_1s
#else
# define CFT_1S cft_1
#endif

#if defined(__FFTW) || defined(__AIX)
# define CFT_2S cft_2s
#else
# define CFT_2S cft_2
#endif
  !

  call start_clock ('cft3s')

  allocate( aux( nxxs ) )

  !
  ! see comments in cft3.F for the logic (or lack of it) of the following
  !

  if (nprocp.eq.1) then
     nppx = nx3
  else
     nppx = npps (me)
  endif

  if (sign.gt.0) then
     if (sign.ne.2) then
#if defined __FFT_MODULE_DRV
        call cft_1z (f, ncps (me), n3, nx3, sign, aux)
#else
        call CFT_1S (f, ncps (me), n3, nx3, sign, aux)
#endif
        call fft_scatter (aux, nx3, nxxs, f, ncps, npps, sign)
        f(:) = (0.d0,0.d0)
        do i = 1, ncts
           mc = dffts%ismap (i)
           do j = 1, npps (me)
              f (mc + (j - 1) * ncplanes) = aux (j + (i - 1) * nppx)
           enddo
        enddo
        do i = 1, nx1
           planes (i) = 1
        enddo
     else
#if defined __FFT_MODULE_DRV
        call cft_1z (f, nkcp (me), n3, nx3, sign, aux)
#else
        call CFT_1S (f, nkcp (me), n3, nx3, sign, aux)
#endif
        call fft_scatter (aux, nx3, nxxs, f, nkcp, npps, sign)
        f(:) = (0.d0,0.d0)
        ii = 0
        do i = 1, nx1
           planes (i) = 0
        enddo
        do iproc = 1, nprocp
           do i = 1, nkcp (iproc)
              mc = dffts%ismap (i + ncp0s (iproc) )
              ii = ii + 1
              k = mod (mc - 1, nx1) + 1
              planes (k) = 1
              do j = 1, npps (me)
                 f (mc + (j - 1) * ncplanes) = aux (j + (ii - 1) * nppx)
              enddo
           enddo
        enddo
     endif
#if defined __FFT_MODULE_DRV
     call cft_2xy (f, npps (me), n1, n2, nx1, nx2, sign, planes)
#else
     call CFT_2S (f, npps (me), n1, n2, nx1, nx2, sign, planes)
#endif
  else
     if (sign.ne. - 2) then
        do i = 1, nx1
           planes (i) = 1
        enddo
     else
        do i = 1, nx1
           planes (i) = 0
        enddo
        do iproc = 1, nprocp
           do i = 1, nkcp (iproc)
              mc = dffts%ismap (i + ncp0s (iproc) )
              k = mod (mc - 1, nx1) + 1
              planes (k) = 1
           enddo
        enddo
     endif
#if defined __FFT_MODULE_DRV
     call cft_2xy (f, npps (me), n1, n2, nx1, nx2, sign, planes)
#else
     call CFT_2S (f, npps (me), n1, n2, nx1, nx2, sign, planes)
#endif
     if (sign.ne. - 2) then
        do i = 1, ncts
           mc = dffts%ismap (i)
           do j = 1, npps (me)
              aux (j + (i - 1) * nppx) = f (mc + (j - 1) * ncplanes)
           enddo
        enddo
        call fft_scatter (aux, nx3, nxxs, f, ncps, npps, sign)
#if defined __FFT_MODULE_DRV
        call cft_1z (aux, ncps (me), n3, nx3, sign, f)
#else
        call CFT_1S (aux, ncps (me), n3, nx3, sign, f)
#endif
     else
        ii = 0
        do iproc = 1, nprocp
           do i = 1, nkcp (iproc)
              mc = dffts%ismap (i + ncp0s (iproc) )
              ii = ii + 1
              do j = 1, npps (me)
                 aux (j + (ii - 1) * nppx) = f (mc + (j - 1) * ncplanes)
              enddo
           enddo
        enddo
        call fft_scatter (aux, nx3, nxxs, f, nkcp, npps, sign)
#if defined __FFT_MODULE_DRV
        call cft_1z (aux, nkcp (me), n3, nx3, sign, f)
#else
        call CFT_1S (aux, nkcp (me), n3, nx3, sign, f)
#endif
     endif
  endif

  deallocate( aux )

  call stop_clock ('cft3s')

  return
end subroutine cft3s

#else

# define NOPENCILS

#if defined __HPM
#  include "/cineca/prod/hpm/include/f_hpm.h"
#endif

!
!----------------------------------------------------------------------
subroutine cft3s (f, n1, n2, n3, nx1, nx2, nx3, sign)
  !----------------------------------------------------------------------
  !
  use parameters

  use fft_scalar, only: cfft3ds, cfft3d    !  common scalar fft driver
  use sticks, only: dffts          !  data structure for fft data layout

  implicit none

  integer :: n1, n2, n3, nx1, nx2, nx3, sign

  complex(kind=DP) :: f (nx1 * nx2 * nx3)

  call start_clock ('cft3s')

#if defined __HPM
            CALL f_hpmstart( 20, 'cft3s' )
#endif

  !
  !   sign = +-1 : complete 3d fft (for rho and for the potential)
  !

  if (sign.eq.1) then

#if defined __FFT_MODULE_DRV
     call cfft3d (f, n1, n2, n3, nx1, nx2, nx3, 1)
#else
     call cft_3 (f, n1, n2, n3, nx1, nx2, nx3, 2, 1)
#endif

  elseif (sign.eq. - 1) then

#if defined __FFT_MODULE_DRV
     call cfft3d (f, n1, n2, n3, nx1, nx2, nx3, - 1)
#else
     call cft_3 (f, n1, n2, n3, nx1, nx2, nx3, 2, - 1)
#endif

     !
     !   sign = +-2 : if available, call the "short" fft (for psi's)
     !

  elseif (sign.eq.2) then

#if defined __FFT_MODULE_DRV && ( defined __AIX || defined __FFTW )
     call cfft3ds (f, n1, n2, n3, nx1, nx2, nx3,  1, dffts%isind, dffts%iplw)
#elif defined __FFT_MODULE_DRV
     call cfft3d (f, n1, n2, n3, nx1, nx2, nx3, 1)
#elif defined NOPENCILS
     call cft_3 (f, n1, n2, n3, nx1, nx2, nx3, 2, 1)
#else
     call cfts_3 (f, n1, n2, n3, nx1, nx2, nx3, 2, 1, dffts%isind, dffts%iplw)
#endif

  elseif (sign.eq. - 2) then

#if defined __FFT_MODULE_DRV && ( defined __AIX || defined __FFTW )
     call cfft3ds (f, n1, n2, n3, nx1, nx2, nx3, -1, dffts%isind, dffts%iplw)
#elif defined __FFT_MODULE_DRV
     call cfft3d (f, n1, n2, n3, nx1, nx2, nx3, -1)
#elif defined NOPENCILS
     call cft_3 (f, n1, n2, n3, nx1, nx2, nx3, 2, - 1)
#else
     call cfts_3 (f, n1, n2, n3, nx1, nx2, nx3, 2, - 1, dffts%isind, dffts%iplw)
#endif

  else

     call errore ('cft3', 'what should i do?', 1)

  endif

#if defined __HPM
            CALL f_hpmstop( 20 )
#endif

  call stop_clock ('cft3s')

  return
end subroutine cft3s


#endif

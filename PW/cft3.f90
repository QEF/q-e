!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
# if defined __AIX || defined __FFTW || defined __SGI
#  define __FFT_MODULE_DRV
# endif
!
#ifdef __PARA
!
!----------------------------------------------------------------------
subroutine cft3 (f, n1, n2, n3, nx1, nx2, nx3, sign)
  !----------------------------------------------------------------------
  !
  !   sign = +-1 : parallel 3d fft for rho and for the potential
  !
  !   sign = +1 : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !               fft along z using pencils (cft_1)
  !               transpose across nodes    (fft_scatter)
  !                  and reorder
  !               fft along y and x         (cft_2)
  !   sign = -1 : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  !               fft along x and y         (cft_2)
  !               transpose across nodes    (fft_scatter)
  !                  and reorder
  !               fft along z using pencils (cft_1)
  !
#include "f_defs.h"

#if defined __FFT_MODULE_DRV
  use fft_scalar, only : cft_1z, cft_2xy
#endif
  use sticks, only: dfftp
  use fft_base, only: fft_scatter
  USE kinds, only : DP
  USE mp_global, ONLY : nproc_pool, me_pool
  use pfft, only: nct, ncp, ncplane, nxx, npp

  implicit none

  integer :: n1, n2, n3, nx1, nx2, nx3, sign

  complex (kind=DP) :: f (nxx)
  integer :: nxx_save, mc, i, j, ii, iproc, nppx
  complex (kind=DP), allocatable  :: aux (:)
  !
  call start_clock ('cft3')

  allocate( aux( nxx ) )

  !
  ! the following is needed if the fft is distributed over only one proces
  ! for the special case nx3.ne.n3. Not an elegant solution, but simple, f
  ! and better than the preceding one that did not work in some cases. Not
  ! that fft_scatter does nothing if nproc_pool=1. PG
  !
  if (nproc_pool.eq.1) then
     nppx = nx3
  else
     nppx = npp (me_pool+1)
  endif
  !
  if (sign.eq.1) then
#if defined __FFT_MODULE_DRV
     call cft_1z (f, ncp (me_pool+1), n3, nx3, sign, aux)
     ! call cft_1z (f, dfftp%nsp(me_pool+1), n3, nx3, sign, aux)
#else
     call cft_1 (f, ncp (me_pool+1), n3, nx3, sign, aux)
#endif
     call fft_scatter (aux, nx3, nxx, f, ncp, npp, sign)
     f(:) = (0.d0,0.d0)
     do i = 1, nct
        mc = dfftp%ismap (i)
        do j = 1, npp (me_pool+1)
           f (mc + (j - 1) * ncplane) = aux (j + (i - 1) * nppx)
        enddo
     enddo
#if defined __FFT_MODULE_DRV
     call cft_2xy (f, npp (me_pool+1), n1, n2, nx1, nx2, sign)
     ! call cft_2xy (f, dfftp%npp (me_pool+1), n1, n2, nx1, nx2, sign)
#else
     call cft_2 (f, npp (me_pool+1), n1, n2, nx1, nx2, sign)
#endif
  elseif (sign.eq. - 1) then
#if defined __FFT_MODULE_DRV
     call cft_2xy (f, npp (me_pool+1), n1, n2, nx1, nx2, sign)
     ! call cft_2xy (f, dfftp%npp (me_pool+1), n1, n2, nx1, nx2, sign)
#else
     call cft_2 (f, npp (me_pool+1), n1, n2, nx1, nx2, sign)
#endif
     do i = 1, nct
        mc = dfftp%ismap (i)
        do j = 1, npp (me_pool+1)
           aux (j + (i - 1) * nppx) = f (mc + (j - 1) * ncplane)
        enddo
     enddo
     call fft_scatter (aux, nx3, nxx, f, ncp, npp, sign)
#if defined __FFT_MODULE_DRV
     call cft_1z (aux, ncp (me_pool+1), n3, nx3, sign, f)
     ! call cft_1z (aux, dfftp%nsp (me_pool+1), n3, nx3, sign, f)
#else
     call cft_1 (aux, ncp (me_pool+1), n3, nx3, sign, f)
#endif
  else
     call errore ('cft3', 'not allowed', abs (sign) )

  endif
  deallocate( aux )
  call stop_clock ('cft3')
  return
end subroutine cft3
#else
!
!----------------------------------------------------------------------
subroutine cft3 (f, n1, n2, n3, nx1, nx2, nx3, sign)
  !----------------------------------------------------------------------
  !
#if defined __FFT_MODULE_DRV
  use fft_scalar, only : cfft3d
#endif
  USE kinds
  implicit none
  integer :: n1, n2, n3, nx1, nx2, nx3, sign

  complex(kind=DP) :: f (nx1 * nx2 * nx3)
  call start_clock ('cft3')
  !
  !   sign = +-1 : complete 3d fft (for rho and for the potential)
  !
  if (sign.eq.1) then
#if defined __FFT_MODULE_DRV
     call cfft3d (f, n1, n2, n3, nx1, nx2, nx3, 1)
#else
     call cft_3 (f, n1, n2, n3, nx1, nx2, nx3, 1, 1)
#endif
  elseif (sign.eq. - 1) then
#if defined __FFT_MODULE_DRV
     call cfft3d (f, n1, n2, n3, nx1, nx2, nx3, - 1)
#else
     call cft_3 (f, n1, n2, n3, nx1, nx2, nx3, 1, - 1)
#endif
  else
     call errore ('cft3', 'what should i do?', 1)
  endif

  call stop_clock ('cft3')
  return
end subroutine cft3
#endif


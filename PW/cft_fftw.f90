!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#ifdef __FFTW
!----------------------------------------------------------------------
subroutine cft_1 (f, m, n, nx, isign, fout)
  !     ===============
  !     driver routine for m 1d complex fft's (dense grid) - fftw
  !     NOTA BENE: not in-place! output in fout
  !----------------------------------------------------------------------
#include "f_defs.h"
  use fftw_mod
  USE kinds, only : DP
  implicit none
  integer :: m, n, nx, isign

  complex (kind=DP) :: f (nx * m), fout (nx * m)
  real (kind=DP) :: fac
  integer :: ibid
  !
  ! initialization variables
  !
  C_POINTER :: plan (2)
  save plan
  data plan / 0, 0 /
  !
  !
  if (isign.eq.1) then
     ibid = 1
  elseif (isign.eq. - 1) then
     ibid = 2
  else
     call errore ('cft_1', 'wrong call', isign)
  endif
  !
  if (plan (ibid) .eq.0) call FFTW_F77_CREATE_PLAN (plan (ibid), &
       n, isign, FFTW_ESTIMATE)
  !
  call FFTW_F77 (plan (ibid), m, f, 1, nx, fout, 1, nx)
  !
  if (isign.eq. - 1) then
     fac = 1.0 / dble (n)
     call DSCAL (2 * nx * m, fac, fout, 1)

  endif
  return
end subroutine cft_1
!
!----------------------------------------------------------------------
subroutine cft_1s (f, m, n, nx, isign, fout)
  !     ===============
  !     driver routine for m 1d complex fft's (dense grid) - fftw
  !     NOTA BENE: not in-place! output in fout
  !----------------------------------------------------------------------
#include "f_defs.h"
  use fftw_mod
  USE kinds, only : DP
  implicit none
  integer :: m, n, nx, isign

  complex (kind=DP) :: f (nx * m), fout (nx * m)
  real (kind=DP) :: fac
  integer :: isign1, ibid
  !
  ! initialization variables
  !
  C_POINTER :: plan (2)
  save plan
  data plan / 0, 0 /
  !
  !
  if (isign.eq.1.or.isign.eq.2) then
     isign1 = 1
     ibid = 1
  elseif (isign.eq. - 1.or.isign.eq. - 2) then
     isign1 = - 1
     ibid = 2
  else
     call errore ('cft_1s', 'wrong call', isign)
  endif
  !
  if (plan (ibid) .eq.0) call FFTW_F77_CREATE_PLAN (plan (ibid), &
       n, isign1, FFTW_ESTIMATE)
  !
  call FFTW_F77 (plan (ibid), m, f, 1, nx, fout, 1, nx)
  !
  if (isign1.eq. - 1) then
     fac = 1.0 / dble (n)
     call DSCAL (2 * nx * m, fac, fout, 1)

  endif
  return

end subroutine cft_1s
!----------------------------------------------------------------------
subroutine cft_2 (f, mplane, n1, n2, nx1, nx2, isign)
  !     ===============
  !     driver routine for mplane 2d complex fft's of lengths n1 and n2
  !     for charge density and potentials -  uses FFTW
  !     nx1=n1+1 is allowed (in order to avoid memory conflicts)
  !     for compatibility: nx2=n2, nx2 is not used
  !
  !----------------------------------------------------------------------
#include "f_defs.h"
  use fftw_mod
  USE kinds, only : DP
  implicit none
  integer :: n1, n2, mplane, nx1, nx2, isign
  complex (kind=DP) :: f (nx1 * nx2 * mplane)
  !
  integer, parameter :: nmax = 256
  complex (kind=DP) :: fout (nmax)
  integer :: ibid, isign1, i, k
  real (kind=DP) :: fac
  external DSCAL
  !
  ! initialization variables
  !
  C_POINTER :: plan1 (2), plan2 (2)
  save plan1, plan2
  data plan1 / 0, 0 /, plan2 / 0, 0 /
  !
  !
  if (n1 > nmax .or. n2 > nmax) &
       call errore ('cft_fftw.f90:cft_2', 'increase nmax', max (n1, n2) )
  if (n2.ne.nx2) call errore ('cft_2', 'no longer implemented', 1)
  if (isign.eq.1.or.isign.eq.2) then
     isign1 = 1
     ibid = 1
  elseif (isign.eq. - 1.or.isign.eq. - 2) then
     isign1 = - 1
     ibid = 2
  else
     call errore ('cft_2', 'wrong call', isign)
  endif
  !
  if (isign1.eq.1) then
     ! j-direction

     if (plan2 (ibid) .eq.0) call FFTW_F77_CREATE_PLAN (plan2 (ibid) &
          , n2, isign1, FFTW_ESTIMATE+FFTW_IN_PLACE)
     do i = 1, n1
        call FFTW_F77 (plan2 (ibid), mplane, f (i), nx1, nx1 * nx2, &
             fout, 0, 0)
     enddo
     ! i-direction

     if (plan1 (ibid) .eq.0) call FFTW_F77_CREATE_PLAN (plan1 (ibid) &
          , n1, isign1, FFTW_ESTIMATE+FFTW_IN_PLACE)
     call FFTW_F77 (plan1 (ibid), n2 * mplane, f, 1, nx1, fout, 1, &
          nx1)
  else
     ! i-direction

     if (plan1 (ibid) .eq.0) call FFTW_F77_CREATE_PLAN (plan1 (ibid) &
          , n1, isign1, FFTW_ESTIMATE+FFTW_IN_PLACE)
     call FFTW_F77 (plan1 (ibid), n2 * mplane, f, 1, nx1, fout, 1, &
          nx1)
     ! j-direction

     if (plan2 (ibid) .eq.0) call FFTW_F77_CREATE_PLAN (plan2 (ibid) &
          , n2, isign1, FFTW_ESTIMATE+FFTW_IN_PLACE)
     do i = 1, n1
        call FFTW_F77 (plan2 (ibid), mplane, f (i), nx1, nx1 * nx2, &
             fout, 0, 0)
     enddo
     !
     fac = 1.0 / dble (n1 * n2)
     call DSCAL (2 * nx1 * nx2 * mplane, fac, f, 1)
     !

  endif
  return

end subroutine cft_2
!----------------------------------------------------------------------

subroutine cft_2s (f, mplane, n1, n2, nx1, nx2, isign, planes)
  !     ===============
  !     driver routine for mplane 2d complex fft's of lenghts n1 and n2
  !     (sparse and wavefunction grid) - fftw
  !----------------------------------------------------------------------
  use fftw_mod
  USE kinds, only : DP
  implicit none
  integer :: n1, n2, mplane, nx1, nx2, isign, planes (nx1)
  complex (kind=DP) :: f (nx1 * nx2 * mplane)
#include "f_defs.h"
  integer, parameter :: nmax = 256
  complex (kind=DP) :: fout (nmax)
  real (kind=DP) :: fac
  integer :: ibid, isign1, i, k, m, istrt
  !
  ! initialization variables
  !
  C_POINTER :: plan1 (2), plan2 (2)
  save plan1, plan2
  data plan1 / 0, 0 /, plan2 / 0, 0 /
  !
  !
  if (n1 > nmax .or. n2 > nmax) &
       call errore ('cft_fftw.f90:cft_2s', 'increase nmax', max (n1, n2) )
  if (n2.ne.nx2) call errore ('cft_2s', 'not implemented', 1)
  if (isign.eq.1.or.isign.eq.2) then
     isign1 = 1
     ibid = 1
  elseif (isign.eq. - 1.or.isign.eq. - 2) then
     isign1 = - 1
     ibid = 2
  else
     call errore ('cft_2s', 'wrong call', isign)
  endif
  ! check how many columns along x are nonzero
  m = 0
  do i = 1, n1
     m = m + planes (i)
  enddo
  if (m.gt.n1.or.m.le.0) call errore ('cft_2s', 'something wrong with planes', 1)
  !
  if (isign1.eq.1) then
     ! j-direction

     if (plan2 (ibid) .eq.0) call FFTW_F77_CREATE_PLAN (plan2 (ibid) &
          , n2, isign1, FFTW_ESTIMATE+FFTW_IN_PLACE)
     do i = 1, n1
        !
        ! do only ffts on columns (i,*,k) resulting in nonzero components
        !
        if (planes (i) .eq.1) then
           call FFTW_F77 (plan2 (ibid), mplane, f (i), nx1, nx1 * nx2, &
                fout, 0, 0)
        endif
     enddo
     ! i-direction

     if (plan1 (ibid) .eq.0) call FFTW_F77_CREATE_PLAN (plan1 (ibid) &
          , n1, isign1, FFTW_ESTIMATE+FFTW_IN_PLACE)
     call FFTW_F77 (plan1 (ibid), n2 * mplane, f, 1, nx1, fout, 1, &
          nx1)
  else
     ! i-direction

     if (plan1 (ibid) .eq.0) call FFTW_F77_CREATE_PLAN (plan1 (ibid) &
          , n1, isign1, FFTW_ESTIMATE+FFTW_IN_PLACE)
     call FFTW_F77 (plan1 (ibid), n2 * mplane, f, 1, nx1, fout, 1, &
          nx1)
     ! j-direction

     if (plan2 (ibid) .eq.0) call FFTW_F77_CREATE_PLAN (plan2 (ibid) &
          , n2, isign1, FFTW_ESTIMATE+FFTW_IN_PLACE)
     do i = 1, n1
        !
        ! do only ffts on columns (i,*,k) resulting in nonzero components
        !
        if (planes (i) .eq.1) then
           call FFTW_F77 (plan2 (ibid), mplane, f (i), nx1, nx1 * nx2, &
                fout, 0, 0)
        endif
     enddo
     !
     fac = 1.0 / dble (n1 * n2)
     call DSCAL (2 * nx1 * nx2 * mplane, fac, f, 1)
     !

  endif
  return
end subroutine cft_2s
#else
subroutine bidon_fftw
  stop 'cft_fftw'
end subroutine bidon_fftw
#endif

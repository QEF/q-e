!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#ifdef __FFTW
subroutine bidon_sp
  stop 'cft_sp'
end subroutine bidon_sp
#else
#ifdef AIX
!----------------------------------------------------------------------

subroutine cft_1 (f, m, n, nx, sgn, fout)
  !     ===============
  !     driver routine for m 1d complex fft's (sp2/t3d only)
  !     nx=n+1 is allowed (in order to avoid memory conflicts)
  !     NOTA BENE: not in-place! output in fout
  !----------------------------------------------------------------------
#include "machine.h"
  implicit none
  integer :: m, n, nx, sgn

  complex (kind=8) :: f (nx * m), fout (nx * m)
  integer :: sign, isign, op (2), om (2), naux1, naux2, itype
  parameter (naux1 = 20000, naux2 = 20000)
  real (kind=8) :: aux1p (naux1, 2), aux1m (naux1, 2), aux2 (naux2), &
       scale
  external DSCAL
  data op, om / 0, 0, 0, 0 /

  save op, om, aux1p, aux1m
  isign = - sign (1, sgn)
  itype = abs (sgn)

  if (itype.le.0.or.itype.gt.2) call error ('cft_1', 'wrong call', &
       1)

  scale = 1.d0

  if (isign.eq.1) then
     if (n.ne.op (itype) ) then
        call dcft (1, f, 1, nx, fout, 1, nx, n, m, isign, scale, &
             aux1p (1, itype), naux1, aux2, naux2)
        op (itype) = n
     endif
     call dcft (0, f, 1, nx, fout, 1, nx, n, m, isign, scale, aux1p &
          (1, itype), naux1, aux2, naux2)
     call DSCAL (2 * nx * m, 1d0 / n, fout, 1)

  elseif (isign.eq. - 1) then
     if (n.ne.om (itype) ) then
        call dcft (1, f, 1, nx, fout, 1, nx, n, m, isign, scale, &
             aux1m (1, itype), naux1, aux2, naux2)
        om (itype) = n
     endif

     call dcft (0, f, 1, nx, fout, 1, nx, n, m, isign, scale, aux1m &
          (1, itype), naux1, aux2, naux2)
  else
     call error ('cft_1', 'wrong sign', 1)

  endif
  return

end subroutine cft_1
!----------------------------------------------------------------------

subroutine cft_2 (f, mplane, n1, n2, nx1, nx2, sgn)
  !     ===============
  !     driver routine for mplane 2d complex fft's of lenghts n1 and n2
  !     nx1 is the actual dimension of f (may differ from n)
  !     for compatibility: nx2=n2, nx2 is not used - sp2 version
  !
  !----------------------------------------------------------------------
#include "machine.h"
  implicit none
  integer :: n1, n2, mplane, nx1, nx2, sgn

  complex (kind=8) :: f (nx1 * nx2 * mplane)
  integer :: isign, itype, o1p, o2p, o1m, o2m, m, incx2, incx1, k, &
       istrt, naux1, naux2
  parameter (naux1 = 20000, naux2 = 20000)
  real (kind=8) :: aux1p (naux1, 2), aux1m (naux1, 2), aux2 (naux2), &
       scale
  external DSCAL
  data o1p, o2p, o1m, o2m / 0, 0, 0, 0 /


  save o1p, o2p, o1m, o2m, aux1p, aux1m
  isign = - sign (1, sgn)
  itype = abs (sgn)

  if (itype.le.0.or.itype.gt.2) call error ('cft_1', 'wrong call', &
       1)
  if (n2.ne.nx2) call error ('cft_2', 'no longer implemented', 1)

  scale = 1.d0


  if (isign.eq.1) then
     !  i - direction ...
     incx1 = 1
     incx2 = nx1
     m = n2 * mplane
     if (n1.ne.o1p) then
        call dcft (1, f, incx1, incx2, f, incx1, incx2, n1, m, &
             isign, scale, aux1p (1, 1), naux1, aux2, naux2)
        o1p = n1

     endif


     call dcft (0, f, incx1, incx2, f, incx1, incx2, n1, m, isign, &
          scale, aux1p (1, 1), naux1, aux2, naux2)
     ! ... j-direction ...
     incx1 = nx1
     incx2 = 1
     m = n1
     if (n2.ne.o2p) then
        call dcft (1, f, incx1, incx2, f, incx1, incx2, n2, m, &
             isign, scale, aux1p (1, 2), naux1, aux2, naux2)
        o2p = n2
     endif
     do k = 1, mplane
        istrt = 1 + (k - 1) * nx1 * n2
        call dcft (0, f (istrt), incx1, incx2, f (istrt), incx1, incx2, &
             n2, m, isign, scale, aux1p (1, 2), naux1, aux2, naux2)

     enddo

     call DSCAL (2 * nx1 * n2 * mplane, 1d0 / (n1 * n2), f, 1)


  elseif (isign.eq. - 1) then
     !     i - direction ...
     incx1 = 1
     incx2 = nx1
     m = n2 * mplane
     if (n1.ne.o1m) then
        call dcft (1, f, incx1, incx2, f, incx1, incx2, n1, m, &
             isign, scale, aux1m (1, 1), naux1, aux2, naux2)
        o1m = n1
     endif


     call dcft (0, f, incx1, incx2, f, incx1, incx2, n1, m, isign, &
          scale, aux1m (1, 1), naux1, aux2, naux2)
     !     ... j-direction ...
     incx1 = nx1
     incx2 = 1
     m = n1
     if (n2.ne.o2m) then
        call dcft (1, f, incx1, incx2, f, incx1, incx2, n2, m, &
             isign, scale, aux1m (1, 2), naux1, aux2, naux2)
        o2m = n2
     endif
     do k = 1, mplane
        istrt = 1 + (k - 1) * nx1 * n2
        call dcft (0, f (istrt), incx1, incx2, f (istrt), incx1, incx2, &
             n2, m, isign, scale, aux1m (1, 2), naux1, aux2, naux2)

     enddo
  else
     call error ('cft_2', 'wrong sign', 1)

  endif
  return

end subroutine cft_2
!----------------------------------------------------------------------

subroutine cft_2s (f, mplane, n1, n2, nx1, nx2, sgn, planes)
  !     ===============
  !     driver routine for mplane 2d complex fft's of lengths n1 and n2
  !     for  wavefunctions (planes used) - uses ESSL
  !     nx1 is the actual dimension of f (may differ from n)
  !     for compatibility: nx2=n2, nx2 is not used
  !
  !----------------------------------------------------------------------
#include "machine.h"
  implicit none
  integer :: n1, n2, mplane, nx1, nx2, sgn, planes (nx1)

  complex (kind=8) :: f (nx1 * nx2 * mplane)
  integer :: isign, itype, o1p, o2p, o1m, o2m, m, incx2, incx1, k, &
       i, istrt, naux1, naux2
  parameter (naux1 = 20000, naux2 = 20000)
  real (kind=8) :: aux1p (naux1, 2), aux1m (naux1, 2), aux2 (naux2), &
       scale
  external DSCAL
  data o1p, o2p, o1m, o2m / 0, 0, 0, 0 /


  save o1p, o2p, o1m, o2m, aux1p, aux1m
  isign = - sign (1, sgn)
  itype = abs (sgn)

  if (itype.le.0.or.itype.gt.2) call error ('cft_2', 'wrong call', &
       1)
  if (n2.ne.nx2) call error ('cft_2', 'no longer implemented', 1)


  scale = 1.d0
  ! check how many columns along x are nonzero
  m = 0
  do i = 1, n1
     m = m + planes (i)
  enddo
  if (m.gt.n1.or.m.le.0) call error ('cft_2', 'something wrong with planes', 1)
  !


  if (isign.eq.1) then
     ! ... i - direction
     incx1 = 1

     incx2 = nx1
     m = n2 * mplane
     if (n1.ne.o1p) then
        call dcft (1, f, incx1, incx2, f, incx1, incx2, n1, m, &
             isign, scale, aux1p (1, 1), naux1, aux2, naux2)
        o1p = n1

     endif


     call dcft (0, f, incx1, incx2, f, incx1, incx2, n1, m, isign, &
          scale, aux1p (1, 1), naux1, aux2, naux2)
     ! j-direction ...
     incx1 = nx1

     incx2 = 1
     m = 1
     if (n2.ne.o2p) then
        call dcft (1, f, incx1, incx2, f, incx1, incx2, n2, m, &
             isign, scale, aux1p (1, 2), naux1, aux2, naux2)
        o2p = n2
     endif
     do i = 1, n1
        !
        ! do only ffts on columns (i,*,k) resulting in nonzero components
        !
        if (planes (i) .eq.1) then
           do k = 1, mplane
              istrt = i + (k - 1) * nx1 * n2
              call dcft (0, f (istrt), incx1, incx2, f (istrt), incx1, &
                   incx2, n2, m, isign, scale, aux1p (1, 2), naux1, aux2, &
                   naux2)
           enddo
        endif

     enddo

     call DSCAL (2 * nx1 * n2 * mplane, 1d0 / (n1 * n2), f, 1)


  elseif (isign.eq. - 1) then
     !     ... j-direction
     incx1 = nx1

     incx2 = 1
     m = 1
     if (n2.ne.o2m) then
        call dcft (1, f, incx1, incx2, f, incx1, incx2, n2, m, &
             isign, scale, aux1m (1, 2), naux1, aux2, naux2)
        o2m = n2
     endif
     do i = 1, n1
        !
        ! do only ffts for columns (i,*,k) having nonzero components
        !
        if (planes (i) .eq.1.or.itype.eq.1) then
           do k = 1, mplane
              istrt = i + (k - 1) * nx1 * n2
              call dcft (0, f (istrt), incx1, incx2, f (istrt), incx1, &
                   incx2, n2, m, isign, scale, aux1m (1, 2), naux1, aux2, &
                   naux2)
           enddo
        endif


     enddo
     !     i - direction ...
     incx1 = 1
     incx2 = nx1
     m = n2 * mplane
     if (n1.ne.o1m) then
        call dcft (1, f, incx1, incx2, f, incx1, incx2, n1, m, &
             isign, scale, aux1m (1, 1), naux1, aux2, naux2)
        o1m = n1
     endif

     call dcft (0, f, incx1, incx2, f, incx1, incx2, n1, m, isign, &
          scale, aux1m (1, 1), naux1, aux2, naux2)
  else
     call error ('cft_2', 'wrong sign', 1)

  endif
  return
end subroutine cft_2s
#else
subroutine bidon_sp
  stop 'cft_sp'
end subroutine bidon_sp
#endif
#endif

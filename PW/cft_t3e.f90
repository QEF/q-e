!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#ifdef __FFTW
subroutine bidon_t3e
  stop 'cft_t3e'
end subroutine bidon_t3e
#else
#ifdef __T3E
!----------------------------------------------------------------------

subroutine cft_1 (f, m, n, nx, sgn, fout)
  !     ===============
  !     driver routine for m 1d complex fft's of lenght n
  !     nx is the actual dimension of f (may differ from n)
  !     t3d/t3e version
  !     NOTA BENE: not in-place! output in fout
  !----------------------------------------------------------------------
#include "machine.h"
  implicit none
  integer :: m, n, nx, sgn

  complex (kind=8) :: f (nx * m), fout (nx * m)
  integer :: sign, isign, on (2), naux1, naux2, itype, i
  parameter (naux1 = 20000, naux2 = 10000)
  !
  ! NOTA BENE: naux1 and naux2 should be checked
  !
  real (kind=8) :: aux1 (naux1, 2), aux2 (naux2), scale
  external DSCAL
  data on / 0, 0 /


  save on, aux1
  isign = sign (1, sgn)
  itype = abs (sgn)

  if (itype.le.0.or.itype.gt.2) call errore ('cft_1', 'wrong call', &
       1)

  scale = 1.d0
  if (n.ne.on (itype) ) then
     call ccfft (0, n, scale, f, fout, aux1 (1, itype), aux2, 0)
     on (itype) = n

  endif
  do i = 1, m
     call ccfft (isign, n, scale, f (1 + (i - 1) * nx), fout (1 + &
          (i - 1) * nx), aux1 (1, itype), aux2, 0)

  enddo

  if (sgn.lt.0) call DSCAL (2 * nx * m, 1d0 / n, fout, 1)
  return

end subroutine cft_1
!----------------------------------------------------------------------

subroutine cft_2 (f, mplane, n1, n2, nx1, nx2, sgn)
  !     ===============
  !     driver routine for mplane 2d complex fft's of lenghts n1 and n2
  !     nx1=n1+1 is allowed (in order to avoid memory conflicts)
  !     for compatibility: nx2=n2, nx2 is not used - t3d/t3e version
  !
  !----------------------------------------------------------------------
#include "machine.h"
  implicit none
  complex (kind=8) :: f (nx1 * nx2 * mplane)
  integer :: n1, n2, mplane, nx1, nx2, sgn
  integer :: sign, isign, o1p, o2p, o1m, o2m, m, k, istrt, naux1, &
       naux2, i
  parameter (naux1 = 20000, naux2 = 10000)
  !
  ! NOTA BENE: naux1 and naux2 should be checked
  !
  real (kind=8) :: aux1p (naux1, 2), aux1m (naux1, 2), aux2 (naux2), &
       scale
  complex (kind=8) :: fj (naux2)
  ! NOTA BENE: fj should be dimensioned fj(2*n2)
  external DSCAL
  data o1p, o2p, o1m, o2m / 0, 0, 0, 0 /

  save o1p, o2p, o1m, o2m, aux1p, aux1m
  isign = - sign (1, sgn)
  if (n2.ne.nx2) call errore ('cft_2', 'no longer implemented', 1)


  scale = 1.d0
  !      WRITE( stdout,*)'in cft_2  ',n1,n2,nx1,nx2


  if (isign.gt.0) then
     !  i - direction ...
     m = n2 * mplane
     if (n1.ne.o1p) then
        call ccfft (0, n1, scale, f, f, aux1p (1, 1), aux2, 0)
        o1p = n1
     endif
     do i = 1, m
        call ccfft ( - isign, n1, scale, f (1 + (i - 1) * nx1), &
             f (1 + (i - 1) * nx1), aux1p (1, 1), aux2, 0)


     enddo
     ! ... j-direction ...
     m = n1
     if (n2.ne.o2p) then
        call ccfft (0, n2, scale, f, f, aux1p (1, 2), aux2, 0)
        o2p = n2
     endif
     do k = 1, mplane
        istrt = 1 + (k - 1) * nx1 * n2
        do i = 1, m
           call ZCOPY (n2, f (istrt + i - 1), nx1, fj, 1)
           call ccfft ( - isign, n2, scale, fj, fj, aux1p (1, 2), aux2, 0)
           call ZCOPY (n2, fj, 1, f (istrt + i - 1), nx1)
        enddo

     enddo

     call DSCAL (2 * nx1 * n2 * mplane, 1d0 / (n1 * n2), f, 1)


  else
     !     i - direction ...
     m = n2 * mplane
     if (n1.ne.o1m) then
        call ccfft (0, n1, scale, f, f, aux1m (1, 1), aux2, 0)
        o1m = n1
     endif
     do i = 1, m
        call ccfft ( - isign, n1, scale, f (1 + (i - 1) * nx1), &
             f (1 + (i - 1) * nx1), aux1m (1, 1), aux2, 0)


     enddo
     !     ... j-direction ...
     m = n1
     if (n2.ne.o2m) then
        call ccfft (0, n2, scale, f, f, aux1m (1, 2), aux2, 0)
        o2m = n2
     endif
     do k = 1, mplane
        istrt = 1 + (k - 1) * nx1 * n2
        do i = 1, m
           call ZCOPY (n2, f (istrt + i - 1), nx1, fj, 1)
           call ccfft ( - isign, n2, scale, fj, fj, aux1m (1, 2), aux2, 0)
           call ZCOPY (n2, fj, 1, f (istrt + i - 1), nx1)
        enddo

     enddo

  endif
  return
end subroutine cft_2
#else
subroutine bidon_t3e
  stop 'cft_t3e'
end subroutine bidon_t3e
#endif
#endif

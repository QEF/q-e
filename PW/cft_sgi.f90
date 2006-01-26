!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

#if defined (__FFTW)
subroutine bidon_sgi
  stop 'cft_sgi'
end subroutine bidon_sgi
#elif defined (__COMPLIB)
!----------------------------------------------------------------------

subroutine cft_1 (f, m, n, nx, sgn, fout)
  !     ===============
  !     driver routine for m 1d complex fft's of lenght n
  !     nx is the actual dimension of f (may differ from n)
  !     sgi origin-2000 version
  !     NOTA BENE: not in-place! output in fout
  !----------------------------------------------------------------------
  USE kinds, only : DP
  implicit none

  integer :: m, n, nx, sgn
  integer :: on (2), naux1, isign, itype, i
  parameter (naux1 = 20000)
  real (DP) :: aux1 (naux1, 2)
  complex (DP) :: f (nx * m), fout (nx * m)
  external DSCAL
  data on / 0, 0 /


  save on, aux1
  isign = sign (1, sgn)
  itype = abs (sgn)

  if (itype.le.0.or.itype.gt.2) call errore ('cft_1', 'wrong call', 1)
  if (n.ne.on (itype) ) then
     call zfft1di (n, aux1 (1, itype) )
     on (itype) = n

  endif
  do i = 1, m
     call zfft1d (isign, n, f (1 + (i - 1) * nx), 1, aux1 (1, itype) )
  enddo
  !
  ! Argh, these sgi fft are only in-place...
  !
  call DCOPY (2 * nx * m, f, 1, fout, 1)

  if (sgn.lt.0) call DSCAL (2 * nx * m, 1d0 / n, fout, 1)
  return
end subroutine cft_1
!
!----------------------------------------------------------------------

subroutine cft_2 (f, mplane, n1, n2, nx1, nx2, sgn)
  !     ===============
  !     driver routine for mplane 2d complex fft's of lenghts n1 and n2
  !     nx1=n1+1 is allowed (in order to avoid memory conflicts)
  !     for compatibility: nx2=n2, nx2 is not used - sgi origin version
  !
  !----------------------------------------------------------------------
  !
  USE kinds, only : DP
  implicit none
  integer :: n1, n2, mplane, nx1, nx2, sgn
  complex (DP) :: f (nx1 * nx2 * mplane)
  !
  integer :: isign, itype, on1(2), on2(2), m, i, k, istrt, naux1, naux2
  parameter (naux1 = 20000, naux2 = 10000)
  real (DP) :: aux1 (naux1, 2, 2), fj (naux2)
  !
  ! NOTA BENE: aux1 should be dimensioned aux1(??)
  !              fj should be dimensioned fj(2*n2)
  !
  external DSCAL
  data on1 / 0, 0 /, on2 / 0, 0 /
  save on1, on2, aux1
  !
  !
  isign = sign (1, sgn)
  if (isign.ne. - 1.and.isign.ne.1) call errore ('cft_2', 'wrong call', 1)
  itype = abs (sgn)

  if (itype.le.0.or.itype.gt.2) call errore ('cft_2', 'wrong call', 2)

  if (n2.ne.nx2) call errore ('cft_2', 'no longer implemented', 1)
  if (n1.ne.on1 (itype) ) then
     call zfft1di (n1, aux1 (1, 1, itype) )
     on1 (itype) = n1
  endiF
  if (n2.ne.on2 (itype) ) then
     call zfft1di (n2, aux1 (1, 2, itype) )
     on2 (itype) = n2


  endif
  !  i - direction ...
  m = n2 * mplane
  do i = 1, m
     call zfft1d (isign, n1, f (1 + (i - 1) * nx1), 1, aux1 (1, 1, &
          itype) )
  enddo
  ! ... j-direction ...
  m = n1
  do k = 1, mplane
     istrt = 1 + (k - 1) * nx1 * n2
     do i = 1, m
        call ZCOPY (n2, f (istrt + i - 1), nx1, fj, 1)
        call zfft1d (isign, n2, fj, 1, aux1 (1, 2, itype) )
        call ZCOPY (n2, fj, 1, f (istrt + i - 1), nx1)
     enddo

  enddo

  if (isign.eq. - 1) call DSCAL (2 * nx1 * n2 * mplane, 1d0 / &
       (n1 * n2), f, 1)
  return
end subroutine cft_2
#elif defined (__SCSL)
! ----------------------------------------------------------------------
     SUBROUTINE cft_1 (f, m, n, nx, sgn, fout)
! ----------------------------------------------------------------------
! Silicon Graphics driver routine for m 1d complex ffts of length n
! nx is the actual dimension of f (may differ from n)
! NOTA BENE: not in-place! output in fout
!
! This version uses the SCSL scienfic library which is both supported
! on Origin/IRIX and Altix/Linux machines. 
!
! Contributed by Martin Hilgeman <hilgeman@sgi.com>, October 2004
! ----------------------------------------------------------------------
!
#include "f_defs.h"
      USE kinds, ONLY : DP
      IMPLICIT NONE
!
      INTEGER                :: m, n, nx, sgn
      COMPLEX(DP)       :: f(nx * m), fout(nx * m)
!
! Local parameters
!
      INTEGER                :: isign, itype, isys(0:1), i
      INTEGER, PARAMETER     :: naux1=20000, nwork=20000
      INTEGER, SAVE          :: on(2)
      REAL(DP), SAVE    :: aux1(naux1, 2)
      REAL(DP)          :: dummy, scale, work(nwork)
      DATA                      on / 0, 0 / 
!
      isys(0) = 1
      isign = sign(1, sgn)
      itype = ABS(sgn)
      IF (itype .LE. 0 .OR. itype .GT. 2) THEN 
         CALL errore('cft_1', 'wrong call', 1)
      END IF
!
      IF ( n .NE. on(itype) ) THEN
         CALL ZZFFTM(0, n, 0, 0.0D0, dummy, 1, dummy, 1,                &
    &               aux1(1, itype), dummy, isys)
         on(itype) = n
      END IF

      IF (isign .GT. 0) THEN
         scale = 1.0D0
      ELSE
         scale = 1.0d0 / n
      END IF
!
      CALL ZZFFTM(isign, n, m, scale, f, nx, fout, nx, aux1(1, itype),  &
     &            work, isys)
!
      RETURN
      END SUBROUTINE cft_1
!
! ----------------------------------------------------------------------
     SUBROUTINE cft_2 (f, mplane, n1, n2, nx1, nx2, sgn)
! ----------------------------------------------------------------------
! Silicon Graphics driver routine for mplane 2d complex ffts of 
! lengths n1 and n2
! nx1=n1+1 is allowed (in order to avoid memory conflicts)
! for compatibility: nx2=n2, nx2 is not used
! NOTA BENE: not in-place! output in fout
!
! This version uses the SCSL scienfic library which is both supported
! on Origin/IRIX and Altix/Linux machines. 
!
! Contributed by Martin Hilgeman <hilgeman@sgi.com>, October 2004
! ----------------------------------------------------------------------
!
#include "f_defs.h"
      USE kinds, ONLY : DP
      IMPLICIT NONE
! 
      INTEGER                :: n1, n2, mplane, nx1, nx2, sgn
      COMPLEX(DP)       :: f(nx1 * nx2 * mplane),                  &
     &                          fout(nx1 * nx2 * mplane)
! 
! Local parameters
!
      INTEGER                :: isign, itype, i, k, istrt, isys(0:1)
      INTEGER, PARAMETER     :: naux1 = 20000, nwork=20000
      INTEGER, SAVE          :: on1(2), on2(2)
      REAL(DP), SAVE    :: aux1(naux1, 2, 2)
      REAL(DP)          :: dummy, scale, work(nwork)
      DATA                      on1 / 0, 0 /, on2 / 0, 0 /
!
! Statements
!
      isys(0) = 1
      isign = sign(1, sgn)
      IF (isign .NE. -1 .AND. isign .NE. 1)                             &
     &    CALL errore('cft_2', 'wrong call', 1)
      itype = ABS(sgn)
      IF (itype .LE. 0 .OR. itype .GT. 2)                               &
     &    CALL errore('cft_2', 'wrong call', 2)
      IF (n2 .NE. nx2) CALL errore('cft_2', 'no longer implemented', 1)
!
      IF ( n1 .NE. on1(itype) ) THEN
         CALL ZZFFTM(0, n1, 0, 0.0D0, dummy, 1, dummy, 1,               &
     &               aux1(1, 1, itype), dummy, isys)
         on1(itype) = n1
      END IF
      IF (n2 .NE. on2(itype)) THEN
         CALL ZZFFTMR(0, n2, 0, 0.0D0, dummy, 1, dummy, 1,              &
     &               aux1(1, 2, itype), dummy, isys)
         on2(itype) = n2
      END IF
!
      IF (isign .GT. 0) THEN
         scale = 1.0D0
      ELSE 
         scale = 1.0D0 / n1
      END IF
!
! i - direction
!
      CALL ZZFFTM(isign, n1, n2 * mplane, scale, f, nx1, f, nx1,        &
     &            aux1(1, 1, itype), work, isys)
!
! j - direction
!

      IF (isign .GT. 0) THEN
         scale = 1.0D0
      ELSE 
         scale = 1.0D0 / n2
      END IF
!
      DO k = 1, mplane
         istrt = 1 + (k - 1) * nx1 * n2
         CALL ZZFFTMR(isign, n2, nx1, scale, f(istrt), nx1, f(istrt),   &
     &                nx1, aux1(1, 2, itype), work, isys)
      END DO
!
      RETURN
      END SUBROUTINE cft_2
#else
subroutine bidon_sgi
  stop 'cft_sgi'
end subroutine bidon_sgi
#endif

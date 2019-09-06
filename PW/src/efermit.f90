!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
FUNCTION efermit( et, nbnd, nks, nelec, nspin, ntetra, tetra, is, isk )
  !--------------------------------------------------------------------
  !! Finds the Fermi energy - tetrahedron method
  !! (see P. E. Bloechl et al, PRB49, 16223 (1994))
  !
  USE io_global,  ONLY: stdout
  USE kinds,      ONLY: DP
  USE constants,  ONLY: rytoev
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: et(nbnd, nks)
  !! the eigenvalues
  INTEGER, INTENT(IN)  :: nbnd
  !! the number of bands
  INTEGER, INTENT(IN)  :: nks
  !! the number of k points
  REAL(DP), INTENT(IN) :: nelec
  !! the number of electrons
  INTEGER, INTENT(IN)  :: nspin
  !! the number of spin components
  INTEGER, INTENT(IN)  :: ntetra
  !! the number of tetrahedra
  INTEGER, INTENT(IN)  :: tetra(4, ntetra)
  !! the vertices of a tetrahedron
  INTEGER, INTENT(IN)  :: is
  !! spin component
  INTEGER, INTENT(IN)  :: isk(nks)
  !! spin components for each k-point
  REAL(DP):: efermit
  !! output: the Fermi energy
  !
  !  ... two parameters
  !
  INTEGER,  PARAMETER :: maxiter = 300
  ! the maximum number of iterations in bisection
  REAL(DP), PARAMETER :: eps= 1.0d-10
  ! a small quantity
  !
  !  ... local variables
  !
  INTEGER :: nlw, ik, iter
  ! the minimum energy band
  ! counter on k points
  ! counter on iterations
  !
  REAL(DP) :: ef, elw, eup, sumkup, sumklw, sumkmid
  ! elw, eup: lower and upper bounds for fermi energy (ef)
  ! sumklw, sumkup: number of states for ef=elw, ef=eup resp.
  ! sumkmid:        number of states for ef=(elw+eup)/2
  REAL(DP), EXTERNAL :: sumkt
  !
  REAL(DP) :: efbetter, better
  !
  !      find bounds for the Fermi energy.
  !
  nlw = MAX( 1, NINT(nelec / 2.0d0 - 5.0d0) )
  elw = et(nlw, 1)
  eup = et(nbnd, 1)
  DO ik = 2, nks
     elw = MIN( elw, et(nlw, ik) )
     eup = MAX( eup, et(nbnd, ik) )
  ENDDO
  !
  !      Bisection method
  !
  sumkup = sumkt( et, nbnd, nks, nspin, ntetra, tetra, eup, is, isk )
  sumklw = sumkt( et, nbnd, nks, nspin, ntetra, tetra, elw, is, isk )
  better = 1.0d+10
  IF ( (sumkup-nelec)<-eps .OR. (sumklw - nelec)>eps )  THEN
     !
     ! this is a serious error and the code should stop here
     ! we don't stop because this may occasionally happen in nonscf
     ! calculations where it may be completely irrelevant
     !
     CALL infomsg( 'efermit', 'internal error, cannot bracket Ef' )
     efermit = better
     RETURN
  END IF
  DO iter = 1, maxiter
     ef = (eup + elw) / 2.d0
     sumkmid = sumkt (et, nbnd, nks, nspin, ntetra, tetra, ef, is, isk)
     IF (ABS(sumkmid-nelec) < better) THEN
        better = ABS(sumkmid-nelec)
        efbetter = ef
     ENDIF
     ! converged
     IF (ABS(sumkmid-nelec) < eps) THEN
        GOTO 100
     ELSEIF ( (sumkmid-nelec) < -eps) THEN
        elw = ef
     ELSE
        eup = ef
     ENDIF
     !
  ENDDO
  !     unconverged exit:
  !     the best available ef is used . Needed in some difficult cases
  ef = efbetter
  sumkmid = sumkt (et, nbnd, nks, nspin, ntetra, tetra, ef, is, isk )
  !
  IF (is /= 0) WRITE(stdout, '(5x,"Spin Component #",i3)') is
  WRITE( stdout, 9010) ef * rytoev, sumkmid
  !     converged exit:
100 CONTINUE
  !     Check if Fermi level is above any of the highest eigenvalues
  DO ik = 1, nks
     IF (is /= 0) THEN
        IF (isk(ik) /= is ) CYCLE
     END IF
     IF (ef > et (nbnd, ik) + 1.d-4) &
          WRITE( stdout, 9020) ef * rytoev, ik, et (nbnd, ik) * rytoev
  ENDDO
  !
  efermit = ef
  !
  RETURN
  !
9010 format (/5x,'Warning: too many iterations in bisection'/ &
       &          5x,'ef = ',f10.6,' sumk = ',f10.6,' electrons')

9020 format (/5x,'Warning: ef =',f10.6, &
       &     ' is above the highest band at k-point',i4,/5x,9x, &
       &     'e  = ',f10.6)
END FUNCTION efermit


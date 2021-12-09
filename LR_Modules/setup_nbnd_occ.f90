!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE setup_nbnd_occ
  !-----------------------------------------------------------------------
  !
  ! This subroutine computes the number of occupied bands for each k point
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : degspin, pi
  USE klist,            ONLY : xk, ltetra, lgauss, degauss, ngauss, nks, &
                               nelec, nelup, neldw, two_fermi_energies, wk
  USE ener,             ONLY : ef
  USE wvfct,            ONLY : nbnd, et
  USE control_lr,       ONLY : nbnd_occ, nbnd_occx
  USE io_global,        ONLY : stdout
  USE noncollin_module, ONLY : noncolin
  USE lsda_mod,         ONLY : lsda, isk
  USE ktetra,           ONLY : tetra_type
  !
  IMPLICIT NONE
  !
  REAL(DP) :: target, small, fac, xmax
  ! auxiliary variables used
  ! to set nbnd_occ in the metallic case
  INTEGER :: ik, ibnd, ipol
  !
  CALL start_clock ('setup_nbnd_occ')
  !
  ALLOCATE ( nbnd_occ(nks) )
  nbnd_occ(:) = 0
  IF (lgauss) THEN
     !
     ! Discard conduction bands such that w0gauss(x,n) < small
     !
     ! hint:
     !   small = 1.0333492677046d-2  ! corresponds to 2 gaussian sigma
     !   small = 6.9626525973374d-5  ! corresponds to 3 gaussian sigma
     !   small = 6.3491173359333d-8  ! corresponds to 4 gaussian sigma
     !
     small = 6.9626525973374d-5
     !
     ! - appropriate limit for gaussian broadening (used for all ngauss)
     !
     xmax = sqrt ( - log (sqrt (pi) * small) )
     !
     ! - appropriate limit for Fermi-Dirac
     !
     IF (ngauss.eq. - 99) THEN
        fac = 1.d0 / sqrt (small)
        xmax = 2.d0 * log (0.5d0 * (fac + sqrt (fac * fac - 4.d0) ) )
     ENDIF
     !
     target = ef + xmax * degauss
     !
     DO ik = 1, nks
        DO ibnd = 1, nbnd
           IF (et(ibnd, ik) .lt.target) nbnd_occ(ik) = ibnd
        ENDDO
        IF (nbnd_occ(ik) .eq. nbnd) WRITE( stdout, '(5x,/,&
             &"Possibly too few bands at point ", i4,3f10.5)') &
             ik,  (xk (ipol, ik) , ipol = 1, 3)
     ENDDO
     !
  ELSE IF (ltetra) THEN
     IF (tetra_type /= 1 .and. tetra_type /= 2) CALL errore &
          ('setup_nbnd_occ','Optimized or linear tetrahedra only', 1)
  ELSE
     !
     IF (noncolin) THEN
        nbnd_occ = nint (nelec)
     ELSE
        IF ( two_fermi_energies ) THEN
           DO ik = 1, nks
              IF (isk(ik)==1) THEN
                 nbnd_occ (ik) = nint (nelup)
              ELSE
                 nbnd_occ (ik) = nint (neldw)
              ENDIF
           ENDDO
        ELSE
           IF (lsda) CALL infomsg('setup_nbnd_occ', &
                                 'Occupation numbers probably wrong')
           DO ik = 1, nks
              nbnd_occ (ik) = nint (nelec) / degspin
           ENDDO
        ENDIF
     ENDIF
     !
  ENDIF
  !
  nbnd_occx = MAXVAL(nbnd_occ)
  !
  CALL stop_clock ('setup_nbnd_occ')
  !
  RETURN
  !
END SUBROUTINE setup_nbnd_occ

!
! Copyright (C) 2016 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine h_prec_freq(ik, omega, h_diag)
   !-----------------------------------------------------------------------
   !
   ! ... Compute the precondition vector h_diag used in the solution of the
   ! ... linear system - On input:
   ! ... ik     index of k-point
   !...  h_diag must be allocated
   ! ... g2kin  contains kinetic energy (k+q+G)^2 for the current k+q point
   !
   ! Similar to h_prec, but works for the finite-frequency Sternheimer equation.
   ! Also, use a different form of preconditioning. It is taken from the
   ! TDDFPT code (lr_sternheimer.f90).
   !
   USE kinds,            ONLY : DP
   USE klist,            ONLY : ngk
   USE wvfct,            ONLY : npwx, nbnd, et
   USE noncollin_module, ONLY : noncolin, npol
   USE qpoint,           ONLY : ikqs, ikks
   USE control_lr,       ONLY : nbnd_occ
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: ik
   !! k point index
   COMPLEX(DP), INTENT(IN) :: omega
   !! Frequency for the Sternheimer equation
   COMPLEX(DP), INTENT(OUT) :: h_diag(npwx*npol, nbnd)
   !
   INTEGER :: ibnd, ig, ikk, ikq, npwq
   !! Counters
   REAL(DP) :: denom
   !! Temporary variable for the denominator of the preconditioner
   REAL(DP), ALLOCATABLE :: h_dia(:,:)
   !! Diagonal part of the Hamiltonian
   REAL(DP), ALLOCATABLE :: s_dia(:,:)
   !! Diagonal part of the overlap matrix (for USPP and PAW)
   !
   ALLOCATE(h_dia(npwx, npol))
   ALLOCATE(s_dia(npwx, npol))
   !
   ikk = ikks(ik)
   ikq = ikqs(ik)
   npwq = ngk(ikq)
   !
   CALL start_clock('h_prec_freq')
   !
   ! Compute the diagonal part of the Hamiltonian and overlap matrix
   !
   CALL usnldiag(npwq, npol, h_dia, s_dia)
   !
   h_diag  = (0.0_DP, 0.0_DP)
   !
   DO ibnd = 1, nbnd_occ(ikk)
      !
      DO ig = 1, npwq
         denom = h_dia(ig, 1) - (et(ibnd, ikk) + omega) * s_dia(ig, 1)
         IF (ABS(denom) < 1.0_DP) denom = 1.0_DP
         h_diag(ig, ibnd) = CMPLX(1.0d0, 0.d0, KIND=DP) / denom
      ENDDO
      !
      IF (noncolin) THEN
         DO ig = 1, npwq
            denom = h_dia(ig, 2) - (et(ibnd, ikk) + omega) * s_dia(ig, 2)
            IF (ABS(denom) < 1.0_DP) denom = 1.0_DP
            h_diag(ig+npwx, ibnd) = CMPLX(1.d0, 0.d0, KIND=DP) / denom
         ENDDO
      ENDIF
   ENDDO
   !
   DEALLOCATE(h_dia)
   DEALLOCATE(s_dia)
   !
   CALL stop_clock( 'h_prec_freq' )
   !
END SUBROUTINE h_prec_freq

!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define TEST_NEW_PRECONDITIONING
!
!----------------------------------------------------------------------------
SUBROUTINE g_psi( lda, n, m, psi, e )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes an estimate of the inverse Hamiltonian
  ! ... and applies it to m wavefunctions
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps4
  USE g_psi_mod, ONLY : h_diag, s_diag
  !
  IMPLICIT NONE
  !
  INTEGER :: lda, n, m
    ! input: the leading dimension of psi
    ! input: the real dimension of psi
    ! input: the number of bands
  REAL(KIND=DP) :: e(m)
    ! input: the eigenvalues
  COMPLEX(KIND=DP) :: psi(lda,m)
    ! input/output: the psi vector
  !
  ! ... Local variables
  !
  REAL(KIND=DP) :: x, scala, denm
  INTEGER :: k, i
    ! counter on psi functions
    ! counter on G vectors
  !
  !
  CALL start_clock( 'g_psi' )
  !
#if defined (TEST_NEW_PRECONDITIONING)
  !
  scala = 1.D0
  !
  DO k = 1, m
     !
     DO i = 1, n
        !
        x = ( h_diag(i) - e(k) * s_diag(i) ) * scala
        !
        denm = ( 1.D0 + x + &
                 SQRT( 1.D0 + ( x - 1.D0 )*( x - 1.D0 ) ) ) / scala
        !
      ! denm = 1.D0 + 16.D0 * x*x*x*x / &
      !        ( 27.D0 + 18.D0 * x + 12.D0 * x*x + 8.D0 * x*x*x )
        !
        psi(i,k) = psi(i,k) / denm
        !
     END DO
     !
  END DO
  !
#else
  !
  DO k = 1, m
     !
     DO i = 1, n
        !
        ! ... denm = g2+v(g=0) - e(k)
        !        
        denm = h_diag(i) - e(k) * s_diag(i)
        !
        denm = SIGN( MAX( ABS( denm ), eps4 ), denm )
        !
        psi(i,k) = psi(i,k) / denm
        !
     END DO
     !
  END DO
  !
#endif
  !
  CALL stop_clock( 'g_psi' )
  !
  RETURN
  !
END SUBROUTINE g_psi

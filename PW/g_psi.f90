!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE g_psi( lda, n, m, psi, e )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes an estimate of the inverse Hamiltonian
  ! ... and applies it to m wavefunctions
  !
  USE kinds, ONLY : DP
  USE g_psi_mod,  ONLY : h_diag, s_diag, test_new_preconditioning
  !
  IMPLICIT NONE
  !
  INTEGER :: lda, n, m
    ! input: the leading dimension of psi
    ! input: the real dimension of psi
    ! input: the number of bands
  REAL(KIND=DP) :: e(m)
    ! input: the eigenvectors
  COMPLEX(KIND=DP) :: psi(lda,m)
    ! inp/out: the psi vector
  !
  ! ... Local variables
  !
  REAL(KIND=DP), PARAMETER :: eps4 = 1.D-4
    ! a small number
  REAL(KIND=DP) :: x, scala, denm
  INTEGER :: k, i
    ! counter on psi functions
    ! counter on G vectors
  !
  !
  CALL start_clock( 'g_psi' )
  !
  IF ( test_new_preconditioning ) THEN
     !
     scala = 1.D0
     !
     DO k = 1, m
        DO i = 1, n
           !
           x = ( h_diag(i) - e(k) * s_diag(i) ) * scala
           !
           denm = ( 1.D0 + x + &
                    SQRT( 1.D0 + ( x - 1.D0 )*( x - 1.D0 ) ) ) / scala
           !
           !         denm = 1.d0 + 16*x*x*x*x/(27.d0+18*x+12*x*x+8*x*x*x)
           !
           psi(i,k) = psi(i,k) / denm
           !
        END DO
     END DO
     !
  ELSE
     !
     DO k = 1, m
        DO i = 1, n
           !
           denm = h_diag(i) - e(k) * s_diag(i)
           !
           ! ... denm = g2+v(g=0) - e(k)
           !
           IF ( ABS( denm ) < eps4 ) denm = SIGN( eps4, denm )
           !
           ! ... denm = sign( max( abs(denm),eps ), denm )
           !
           psi(i,k) = psi(i,k) / denm
           !
        END DO
     END DO
     !
  END IF
  !
  CALL stop_clock( 'g_psi' )
  !
  RETURN
  !
END SUBROUTINE g_psi

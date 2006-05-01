!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE cgramg1( lda, nvecx, n, start, finish, psi, spsi, hpsi )
  !----------------------------------------------------------------------------
  !
  ! ... This routine orthogonalizes several vectors with the method of
  ! ... Gram-Schmidt and imposing that  <psi_i|S|psi_j> = delta_ij.
  ! ... It receives on input the psi and the spsi.
  ! ... It updates also the Hamiltonian so that it contains the new hpsi.
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps8
  USE io_global, ONLY : stdout
  USE wvfct,     ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: &
        lda, n, nvecx, start, finish
    ! input: leading dimension of the vectors
    ! input: physical dimension
    ! input: dimension of psi
    ! input: first vector to orthogonalize
    ! input: last vector to orthogonalize
  COMPLEX(DP), INTENT(INOUT) :: &
        psi(lda,nvecx), spsi(lda,nvecx), hpsi(lda,nvecx)
    ! input/output: the vectors to be orthogonalized
  !
  INTEGER :: vec, vecp
    ! counter on vectors
    ! counter on vectors
  REAL(DP) :: psi_norm
    ! the norm of a vector
  REAL(DP), EXTERNAL :: DDOT
    ! function computing the dot product of two vectros
  !
  !
  CALL start_clock( 'cgramg1' )
  !
  IF ( gamma_only ) THEN
     !
     CALL cgramg1_gamma()
     !
  ELSE
     !
     CALL cgramg1_k()
     !
  END IF
  !
  CALL stop_clock( 'cgramg1' )
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE cgramg1_gamma()
       !-----------------------------------------------------------------------
       !
       USE gvect, ONLY : gstart
       !
       IMPLICIT NONE
       !
       REAL(DP), ALLOCATABLE :: ps(:)
       !
       !
       ALLOCATE( ps( finish ) )
       !
       DO vec = start, finish
          !
          IF ( vec > 1 ) THEN
             !
             CALL DGEMV( 'T', 2*n, vec-1, 2.D0, &
                         psi, 2*lda, spsi(1,vec), 1, 0.D0, ps, 1 )
             !
             IF ( gstart == 2 ) &
                ps(1:vec-1) = ps(1:vec-1) - psi(1,1:vec-1)*spsi(1,vec)
             !
             CALL reduce( ( vec - 1 ), ps )
             !
             DO vecp = 1, ( vec - 1 )
                !
                psi(:,vec)  = psi(:,vec)  - ps(vecp) * psi(:,vecp)
                hpsi(:,vec) = hpsi(:,vec) - ps(vecp) * hpsi(:,vecp)
                spsi(:,vec) = spsi(:,vec) - ps(vecp) * spsi(:,vecp)
                !
             END DO
             !
          END IF
          !
          psi_norm = 2.D0*DDOT( 2*n, psi(1,vec), 1, spsi(1,vec), 1 )
          !
          IF ( gstart == 2 ) psi_norm = psi_norm - psi(1,vec) * spsi(1,vec)
          !
          CALL reduce( 1, psi_norm )
          !
          IF ( psi_norm < eps8 ) THEN
             !
             WRITE( stdout, '(/,5X,"norm = ",F16.10,I4,/)' ) psi_norm, vec
             !
             CALL errore( 'cgramg1_gamma', 'negative or zero norm in S ', 1 )
             !
          END IF
          !
          psi_norm = 1.D0 / SQRT( psi_norm )
          !
          psi(:,vec)  = psi_norm * psi(:,vec)
          hpsi(:,vec) = psi_norm * hpsi(:,vec)
          spsi(:,vec) = psi_norm * spsi(:,vec)
          !
       END DO
       !
       DEALLOCATE( ps )
       !
       RETURN
       !
     END SUBROUTINE cgramg1_gamma
     !
     !-----------------------------------------------------------------------
     SUBROUTINE cgramg1_k()
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       COMPLEX(DP), ALLOCATABLE :: ps(:)
       !
       !
       ALLOCATE( ps( finish ) )
       !
       DO vec = start, finish
          !
          IF ( vec > 1 ) THEN
             !
             CALL ZGEMV( 'C', n, vec-1, ONE, &
                         psi(1,1), lda, spsi(1,vec), 1, ZERO, ps, 1 )
             !
             CALL reduce( 2*( vec - 1 ), ps )
             !
             DO vecp = 1, vec - 1
                !
                psi(:,vec)  = psi(:,vec)  - ps(vecp)*psi(:,vecp)
                hpsi(:,vec) = hpsi(:,vec) - ps(vecp)*hpsi(:,vecp)
                spsi(:,vec) = spsi(:,vec) - ps(vecp)*spsi(:,vecp)
                !
             END DO
             !
          END IF
          !
          psi_norm = DDOT( 2*n, psi(1,vec), 1, spsi(1,vec), 1 )
          !
          CALL reduce( 1, psi_norm )
          !
          IF ( psi_norm < eps8 ) THEN
             !
             WRITE( stdout, '(/,5X,"norm = ",F16.10,I4,/)' ) psi_norm, vec
             !
             CALL errore( 'cgramg1_k', 'negative or zero norm in S', 1 )
             !
          END IF
          !
          psi_norm = 1.D0 / SQRT( psi_norm )
          !
          psi(:,vec)  = psi_norm * psi(:,vec)
          hpsi(:,vec) = psi_norm * hpsi(:,vec)
          spsi(:,vec) = psi_norm * spsi(:,vec)
          !
       END DO
       !
       DEALLOCATE( ps )
       !
       RETURN
       !
     END SUBROUTINE cgramg1_k
     !
END SUBROUTINE cgramg1

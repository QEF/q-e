!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
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
  !
  IMPLICIT NONE
  !
  ! ... first the dummy variables
  !
  INTEGER :: lda, n, nvecx, start, finish
    ! input: leading dimension of the vectors
    ! input: physical dimension
    ! input: dimension of psi
    ! input: first vector to orthogonalize
    ! input: last vector to orthogonalize
  COMPLEX(KIND=DP) :: psi(lda,nvecx), spsi(lda,nvecx), hpsi(lda,nvecx)
    ! input/output: the vectors to be orthogonalized
  !
  ! ... parameters
  !
  INTEGER, PARAMETER :: ierrx = 3
    ! maximum number of errors
  !
  ! ... here the local variables
  !
  INTEGER :: vec, vecp, ierr
    ! counter on vectors
    ! counter on vectors
    ! counter on errors
  COMPLEX(KIND=DP), ALLOCATABLE :: ps(:)
    ! the scalar products
  COMPLEX(KIND=DP) ::  ZDOTC
    ! function which computes scalar products
  REAL(KIND=DP) :: norm, DDOT
    ! the norm of a vector
    ! function computing the dot product of two v
  !
  !
  CALL start_clock( 'cgramg1' )
  !
  ALLOCATE( ps( finish ) )    
  !
  ierr = 0
  !
  DO vec = start, finish
     !
     DO vecp = 1, vec
        !
        ps(vecp) = ZDOTC( n, psi(1,vecp), 1, spsi(1,vec), 1 )
        !
     END DO
     !
     CALL reduce( 2 * vec, ps )
     !
     DO vecp = 1, ( vec - 1 )
        !
        psi(:,vec)  = psi(:,vec)  - ps(vecp) * psi(:,vecp)
        hpsi(:,vec) = hpsi(:,vec) - ps(vecp) * hpsi(:,vecp)
        spsi(:,vec) = spsi(:,vec) - ps(vecp) * spsi(:,vecp)
        !
     END DO
     !
     norm = ps(vec) - DDOT( 2 * ( vec - 1 ), ps, 1, ps, 1 )
     !
     IF ( norm < 0.D0 ) THEN
        !
        WRITE( stdout, '(/,5X,"norm = ",F16.10,I4,/)' ) norm, vec
        !
        CALL errore( 'cgramg1', ' negative norm in S ', 1 )
        !
     END IF
     !
     IF ( norm < eps8 ) THEN
        !
        norm = 1.D0 / SQRT( norm )
        !
        psi(:,vec)  = norm * psi(:,vec)
        hpsi(:,vec) = norm * hpsi(:,vec)
        spsi(:,vec) = norm * spsi(:,vec)
        !
        ierr = ierr + 1
        !
        IF ( ierr <= ierrx ) CYCLE
        !
        CALL errore( 'cgramg1', ' absurd correction vector', vec )
        !
     ELSE
        !
        norm = 1.D0 / SQRT( norm )
        !
        psi(:,vec)  = norm * psi(:,vec)
        hpsi(:,vec) = norm * hpsi(:,vec)
        spsi(:,vec) = norm * spsi(:,vec)
        !   
     END IF
     !
  END DO
  !
  DEALLOCATE( ps )
  !  
  CALL stop_clock( 'cgramg1' )
  !
  RETURN
  !
END SUBROUTINE cgramg1

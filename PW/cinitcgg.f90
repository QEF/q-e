!
! Copyright (C) 2001-2005 PWSCF group
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
SUBROUTINE cinitcgg( ndmx, ndim, nstart, nbnd, psi, evc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Hamiltonian diagonalization in the subspace spanned
  ! ... by nstart states psi (atomic or random wavefunctions).
  ! ... Produces on output nbnd eigenvectors (nbnd <= nstart) in evc.
  ! ... Minimal memory use - evc and psi may overlap
  ! ... Calls h_1psi to calculate H|psi>, S|psi
  !
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : noncolin, npol
  USE bp,               ONLY : lelfield
  !
  IMPLICIT NONE
  !
  INTEGER :: ndim, ndmx, nstart, nbnd
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! k-point considered
  COMPLEX(DP) :: psi(ndmx,npol,nstart), evc(ndmx,npol,nbnd)
    ! input and output eigenvectors (may overlap) 
  REAL(DP) :: e(nbnd)
    ! eigenvalues
  !
  ! ... local variables
  !
  INTEGER                        :: m, i, j, ipol
  INTEGER                        :: kdim, kdmx, kdim2, kdmx2
  REAL (DP)                 :: rtmp(2)
  COMPLEX (DP)              :: ctmp(2)
  COMPLEX (DP), ALLOCATABLE :: aux(:,:,:), hc(:,:,:), sc(:,:)
  REAL (DP),    ALLOCATABLE :: en(:)
  !
  !
  CALL start_clock( 'wfcrot1' )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = ndim
     kdmx = ndmx
     !
  ELSE
     !
     kdim = ndmx * npol
     kdmx = ndmx * npol
     !
  END IF
  !
  kdim2 = 2 * kdim
  kdmx2 = 2 * kdmx
  !
  ALLOCATE( aux( ndmx, npol, 2 ) )
  ALLOCATE( hc( nstart, nstart, 2 ) )
  ALLOCATE( sc( nstart, nstart ) )
  ALLOCATE( en( nstart ) )
  !
  aux = ZERO
  hc  = ZERO
  sc  = ZERO
  en  = 0.D0
  !
  ! ... Set up the Hamiltonian and Overlap matrix
  !
  DO m = 1, nstart
     !
     CALL h_1psi( ndmx, ndim, psi(1,1,m), aux(1,1,1), aux(1,1,2) )
     !
     IF(lelfield) call h_epsi_her(ndmx,ndim,1, psi(1,1,m), aux(1,1,1))
     !
     CALL DGEMV( 'T', kdim2, 2, 1.D0, aux(1,1,1), &
                 kdmx2, psi(1,1,m), 1, 0.D0, rtmp, 1 )
     !
     hc(m,m,1) = CMPLX( rtmp(1), 0.D0 )
     sc(m,m)   = CMPLX( rtmp(2), 0.D0 )
     !
     DO j = m + 1, nstart
        !
        CALL ZGEMV( 'C', kdim, 2, ONE, aux(1,1,1), &
                    kdmx, psi(1,1,j), 1, ZERO, ctmp, 1 )
        !
        hc(m,j,1) = ctmp(1)
        sc(m,j)   = ctmp(2)
        !
        hc(j,m,1) = CONJG( ctmp(1) )
        sc(j,m)   = CONJG( ctmp(2) )
        !
     END DO
     !
  END DO
  !
  CALL reduce( 2 * nstart * nstart, hc(1,1,1) )
  CALL reduce( 2 * nstart * nstart, sc(1,1) )
  !
  ! ... diagonalize
  !
  CALL cdiaghg( nstart, nbnd, hc, sc, nstart, en, hc(1,1,2) )
  !
  e(1:nbnd) = en(1:nbnd)
  !
  ! ... update the basis set
  !
  DO ipol = 1, npol
     !
     DO i = 1, ndim
        !
        DO m = 1, nbnd
           !
           hc(m,1,1) = SUM( hc(:,m,2) * psi(i,ipol,:) )
           !
        END DO
        !
        evc(i,ipol,1:nbnd) = hc(1:nbnd,1,1)
        !
     END DO
     !
  END DO
  !
  DEALLOCATE( en )
  DEALLOCATE( sc )
  DEALLOCATE( hc )
  DEALLOCATE( aux )
  !
  CALL stop_clock( 'wfcrot1' )
  !
  RETURN
  !
END SUBROUTINE cinitcgg

!
! Copyright (C) 2002-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE rinitcgg( npwx, npw, nstart, nbnd, psi, evc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Hamiltonian diagonalization in the subspace spanned
  ! ... by nstart states psi (atomic or random wavefunctions).
  ! ... Produces on output nbnd eigenvectors (nbnd <= nstart) in evc.
  ! ... Minimal memory use - evc and psi may overlap
  ! ... Calls h_1psi to calculate H|psi>, S|psi
  !
  USE kinds, ONLY : DP
  USE gvect, ONLY : gstart
  !
  IMPLICIT NONE
  !
  INTEGER :: npw, npwx, nstart, nbnd
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
  COMPLEX (DP) :: psi(npwx,nstart), evc(npwx,nbnd)
    ! input and output eigenvectors (may overlap) 
  REAL(DP) :: e(nbnd)
    ! eigenvalues
  !
  INTEGER                   :: m, i, j, npw2, npwx2
  REAL (DP)                 :: rtmp(2)
  COMPLEX (DP), ALLOCATABLE :: aux(:,:)
  COMPLEX (DP), ALLOCATABLE :: ctmp(:)
  REAL (DP),    ALLOCATABLE :: hr(:,:,:), sr(:,:)
  REAL (DP),    ALLOCATABLE :: en(:)
  !
  !
  CALL start_clock( 'wfcrot1' )
  !
  npw2  = 2*npw
  npwx2 = 2*npwx
  !
  ALLOCATE( aux( npwx, 2 ) )
  ALLOCATE( ctmp( nbnd ) )
  ALLOCATE( hr( nstart, nstart, 2 ) )
  ALLOCATE( sr( nstart, nstart ) )
  ALLOCATE( en( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix
  !
  DO m = 1, nstart
     !
     CALL h_1psi( npwx, npw, psi(1,m), aux(1,1), aux(1,2) )
     !
     CALL DGEMV( 'T', npw2, 2, 2.D0, aux, npwx2, psi(1,m), 1, 0.D0, rtmp, 1 )
     !
     IF ( gstart == 2 ) rtmp(:) = rtmp(:) - psi(1,m)*aux(1,:)
     !
     hr(m,m,1) = rtmp(1)
     sr(m,m)   = rtmp(2)
     !
     DO j = m + 1, nstart
        !
        CALL DGEMV( 'T', npw2, 2, 2.D0, aux, npwx2, psi(1,j), 1, 0.D0, rtmp, 1 )
        !
        IF ( gstart == 2 ) rtmp(:) = rtmp(:) - psi(1,j)*aux(1,:)
        !
        hr(j,m,1) = rtmp(1)
        sr(j,m)   = rtmp(2)
        !
        hr(m,j,1) = rtmp(1)
        sr(m,j)   = rtmp(2)
        !
     END DO
     !
  END DO
  !
  CALL reduce( nstart*nstart, hr(1,1,1) )
  CALL reduce( nstart*nstart, sr(1,1) )
  !
  ! ... diagonalize
  !
  CALL rdiaghg( nstart, nbnd, hr, sr, nstart, en, hr(1,1,2) )
  !
  e(1:nbnd) = en(1:nbnd)
  !
  ! ... update the basis set
  !
  DO i = 1, npw
     !
     DO m = 1, nbnd
        !
        ctmp(m) = SUM( hr(:,m,2)*psi(i,:) )
        !
     END DO
     !
     evc(i,1:nbnd) = ctmp(1:nbnd)
     !
  END DO
  !
  DEALLOCATE( en )
  DEALLOCATE( sr )
  DEALLOCATE( hr )
  DEALLOCATE( ctmp )
  DEALLOCATE( aux )
  !
  CALL stop_clock( 'wfcrot1' )
  !
  RETURN
  !
END SUBROUTINE rinitcgg

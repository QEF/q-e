!
! Copyright (C) 2001-2003 PWSCF group
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
  COMPLEX (KIND=DP) :: psi(npwx,nstart), evc(npwx,nbnd)
    ! input and output eigenvectors (may overlap) 
  REAL(KIND=DP) :: e(nbnd)
    ! eigenvalues
  !
  !... local variables
  !
  INTEGER                        :: m, ibnd, i, j, npw2
  COMPLEX (KIND=DP), ALLOCATABLE :: hpsi(:), spsi(:)
  REAL (KIND=DP),    ALLOCATABLE :: hr(:,:,:), sr(:,:)
  REAL (KIND=DP),    ALLOCATABLE :: en(:)
  !
  REAL (KIND=DP), EXTERNAL :: DDOT
  !
  !
  CALL start_clock( 'wfcrot1' )
  !
  npw2 = 2 * npw
  !
  ALLOCATE( spsi( npwx ) )    
  ALLOCATE( hpsi( npwx ) )    
  ALLOCATE( hr( nstart, nstart, 2 ) )   
  ALLOCATE( sr( nstart, nstart ) )
  ALLOCATE( en( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix
  !
  DO m = 1, nstart
     !
     CALL h_1psi( npwx, npw, psi(1,m), hpsi, spsi )
     !
     hr(m,m,1) = 2.D0 * DDOT( npw2, psi(1,m), 1, hpsi, 1 )
     sr(m,m)   = 2.D0 * DDOT( npw2, psi(1,m), 1, spsi, 1 )
     !
     IF ( gstart == 2 ) THEN
        !
        hr(m,m,1) = hr(m,m,1) - psi(1,m) * hpsi(1)
        sr(m,m)   = sr(m,m)   - psi(1,m) * spsi(1)
        !
     END IF
     !
     DO j = m + 1, nstart
        !
        hr(j,m,1) = 2.D0 * DDOT( npw2, psi(1,j), 1, hpsi, 1 )
        sr(j,m)   = 2.D0 * DDOT( npw2, psi(1,j), 1, spsi, 1 )
        !
        IF ( gstart == 2 ) THEN
           !
           hr(j,m,1) = hr(j,m,1) - psi(1,j) * hpsi(1)
           sr(j,m)   = sr(j,m)   - psi(1,j) * spsi(1)
           !
        END IF
        !
        hr(m,j,1) = hr(j,m,1) 
        sr(m,j)   = sr(j,m)
        !
     END DO
     !
  END DO
  !
  CALL reduce( 2 * nstart * nstart, hr(1,1,1) )
  CALL reduce( 2 * nstart * nstart, sr(1,1) )
  !
  ! ... diagonalize
  !
  CALL rdiaghg( nstart, nbnd, hr, sr, nstart, en, hr(1,1,2) )
  !
  e(1:nbnd) = en(1:nbnd)
  !
  ! ... update the basis set (hpsi is used as workspace)
  !
  DO m = 1, nbnd
     !
     hpsi = ( 0.D0, 0.D0 )
     !
     DO ibnd = 1, nbnd
        !
        hpsi = hpsi + hr(ibnd,m,2) * psi(:,ibnd)
        !
     END DO
     !
     evc(:,m) = hpsi
     !
  END DO
  !
  DEALLOCATE( en )
  DEALLOCATE( sr )
  DEALLOCATE( hr )
  DEALLOCATE( hpsi )
  DEALLOCATE( spsi )
  !
  CALL stop_clock( 'wfcrot1' )
  !
  RETURN
  !
END SUBROUTINE rinitcgg

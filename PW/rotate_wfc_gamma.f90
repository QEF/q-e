!
! Copyright (C) 2003-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#include "f_defs.h"
!----------------------------------------------------------------------------
SUBROUTINE rotate_wfc_gamma( npwx, npw, nstart, gstart, nbnd, &
                             psi, overlap, evc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Hamiltonian diagonalization in the subspace spanned
  ! ... by nstart states psi (atomic or random wavefunctions).
  ! ... Produces on output nbnd eigenvectors (nbnd <= nstart) in evc.
  ! ... (nbnd <= nvec). Calls h_psi, s_psi to calculate H|psi> and S|psi>.
  ! ... This version assumes real wavefunctions (k=0) with only
  ! ... half plane waves stored: psi(-G)=psi*(G), except G=0
  ! ... It only uses an auxiliary array of the same size as psi.
  !
  USE kinds, ONLY : DP
  USE wvfct, ONLY : gamma_only 
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER :: npw, npwx, nstart, nbnd, gstart
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
  LOGICAL :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP) :: psi(npwx,nstart), evc(npwx,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP) :: e(nbnd)
    ! eigenvalues
  !
  ! ... auxiliary variables:
  !
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  REAL(DP),    ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:), en(:)
  !
  !
  ALLOCATE( aux(  npwx, nstart ) )    
  ALLOCATE( hr( nstart, nstart ) )    
  ALLOCATE( sr( nstart, nstart ) )    
  ALLOCATE( vr( nstart, nstart ) )    
  ALLOCATE( en( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  CALL h_psi( npwx, npw, nstart, psi, aux )
  !
  CALL DGEMM( 'T', 'N', nstart, nstart, 2 * npw, 2.D0 , psi, &
              2 * npwx, aux, 2 * npwx, 0.D0, hr, nstart )
  !  
  IF ( gstart == 2 ) &
     call DGER( nstart, nstart, -1.D0, psi, 2 * npwx, aux, &
                2 * npwx, hr, nstart )
  !     
#if defined (__PARA)
  CALL reduce( nstart * nstart, hr )
#endif
  !
  IF ( overlap ) THEN
     !
     CALL s_psi( npwx, npw, nstart, psi, aux )
     !
     CALL DGEMM( 'T', 'N', nstart, nstart, 2 * npw, 2.D0 , psi, &
                 2 * npwx, aux, 2 * npwx, 0.D0, sr, nstart )
     !            
     IF ( gstart == 2 ) &
        CALL DGER( nstart, nstart, -1.D0, psi, 2 * npwx, &
                   aux, 2 * npwx, sr, nstart )
     !              
  ELSE
     !
     CALL DGEMM( 'T', 'N', nstart, nstart, 2 * npw, 2.D0, psi, &
                 2 * npwx, psi, 2 * npwx, 0.D0, sr, nstart )
     !
     IF ( gstart == 2 ) &
        CALL DGER( nstart, nstart, -1.D0, psi, 2 * npwx, &
                   psi, 2 * npwx, sr, nstart )
     !
  END IF
  !
#if defined (__PARA)
  CALL reduce( nstart * nstart, sr )
#endif
  !
  ! ... Diagonalize
  !
  CALL rdiaghg( nstart, nbnd, hr, sr, nstart, en, vr )
  !
  e(:) = en(1:nbnd)
  !
  ! ... update the basis set
  !
  CALL DGEMM( 'N', 'N', 2 * npw, nbnd, nstart, 1.D0, psi, 2 * npwx, &
              vr, nstart, 0.D0, aux, 2 * npwx )
  !   
  evc(:,:) = aux(:,1:nbnd)
  !
  DEALLOCATE( en )
  DEALLOCATE( vr )
  DEALLOCATE( sr )
  DEALLOCATE( hr )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE rotate_wfc_gamma

!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
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
  !
  USE parameters, ONLY : DP
  USE wvfct,      ONLY : gamma_only 
  !
  IMPLICIT NONE
  !
  INTEGER          :: npw, npwx, nstart, nbnd, gstart
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix psi, as declared in the calling pgm unit
  ! input number of states
  ! output number of states
  LOGICAL          :: overlap
  ! if .FALSE. : S|psi> not needed
  COMPLEX(KIND=DP) :: psi(npwx,nstart), evc(npwx,nbnd)
  ! input and output eigenvectors (may overlap)
  REAL(KIND=DP)    :: e(nbnd)
  ! eigenvalues
  !
  ! ... auxiliary variables:
  !
  COMPLEX(KIND=DP), ALLOCATABLE :: hpsi(:,:), spsi(:,:)
  REAL(KIND=8), ALLOCATABLE     :: hr(:,:), sr(:,:), vr(:,:), en(:)
  !
  !
  ALLOCATE( hpsi(  npwx, nstart ) )    
  IF ( overlap ) &
       ALLOCATE( spsi(  npwx, nstart ) )    
  ALLOCATE( hr( nstart, nstart ) )    
  ALLOCATE( sr( nstart, nstart ) )    
  ALLOCATE( vr( nstart, nstart ) )    
  ALLOCATE( en( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix
  !
  CALL h_psi( npwx, npw, nstart, psi, hpsi )
  !
  IF ( overlap ) CALL s_psi( npwx, npw, nstart, psi, spsi )
  !
  CALL DGEMM( 'T', 'N', nstart, nstart, 2 * npw, 2.D0 , psi, &
              2 * npwx, hpsi, 2 * npwx, 0.D0, hr, nstart )
  !  
  IF ( gstart == 2 ) &
     call DGER( nstart, nstart, -1.D0, psi, 2 * npwx, hpsi, &
                2 * npwx, hr, nstart )
  !     
#ifdef __PARA
  CALL reduce( nstart * nstart, hr )
#endif
  !
  IF ( overlap ) THEN
     CALL DGEMM( 'T', 'N', nstart, nstart, 2 * npw, 2.D0 , psi, &
                 2 * npwx, spsi, 2 * npwx, 0.D0, sr, nstart )
     IF ( gstart == 2 ) &
        CALL DGER( nstart, nstart, -1.D0, psi, 2 * npwx, &
                   spsi, 2 * npwx, sr, nstart )
  ELSE
     CALL DGEMM( 'T', 'N', nstart, nstart, 2 * npw, 2.D0, psi, &
                 2 * npwx, psi, 2 * npwx, 0.D0, sr, nstart )
     IF ( gstart == 2 ) &
        CALL DGER( nstart, nstart, -1.D0, psi, 2 * npwx, &
                   psi, 2 * npwx, sr, nstart )
  END IF
  !
#ifdef __PARA
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
              vr, nstart, 0.D0, hpsi, 2 * npwx )
  !   
  evc(:,:) = hpsi(:,1:nbnd)
  !
  DEALLOCATE( en )
  DEALLOCATE( vr )
  DEALLOCATE( sr )
  DEALLOCATE( hr )
  IF ( overlap ) &
     DEALLOCATE( spsi )
  DEALLOCATE( hpsi )
  !
  RETURN
  !
END SUBROUTINE rotate_wfc_gamma  


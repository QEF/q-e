!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_wfc( npwx, npw, nstart, nbnd, psi, overlap, evc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Hamiltonian diagonalization in the subspace spanned
  ! ... by nstart states psi ( atomic or random wavefunctions ).
  ! ... produces on output nbnd eigenvectors ( nbnd <= nstart ) in evc.
  ! ... Calls h_psi, s_psi to calculate H|psi> ans S|psi>
  ! ... It only uses an auxiliary array of the same size as psi.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER :: npw, npwx, nstart, nbnd
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
  LOGICAL :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(KIND=DP) :: psi(npwx,nstart), evc(npwx,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(KIND=DP) :: e(nbnd)
    ! eigenvalues
  !
  ! ... local variables
  !
  COMPLEX(KIND=DP), ALLOCATABLE :: aux(:,:), hc(:,:), sc(:,:), vc(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: en(:)
  !
  !
  CALL start_clock( 'wfcrot' )
  !
  ALLOCATE( aux(  npwx, nstart ) )    
  ALLOCATE( hc( nstart, nstart) )    
  ALLOCATE( sc( nstart, nstart) )    
  ALLOCATE( vc( nstart, nstart) )    
  ALLOCATE( en( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  CALL h_psi( npwx, npw, nstart, psi, aux )
  !
  call ZGEMM( 'C', 'N', nstart, nstart, npw, ( 1.D0, 0.D0 ), psi, npwx, &
              aux, npwx, ( 0.D0, 0.D0 ), hc, nstart )
  !            
#if defined (__PARA)
  CALL reduce( 2 * nstart * nstart, hc )
#endif
  !
  IF ( overlap ) THEN
     !
     CALL s_psi( npwx, npw, nstart, psi, aux )
     !
     CALL ZGEMM( 'C', 'N', nstart, nstart, npw, ( 1.D0, 0.D0 ), psi, npwx, &
                 aux, npwx, ( 0.D0, 0.D0 ), sc, nstart )
     !
  ELSE
     !
     CALL ZGEMM( 'C', 'N', nstart, nstart, npw, ( 1.D0, 0.D0 ), psi, npwx, &
                 psi, npwx, ( 0.D0, 0.D0 ), sc, nstart )
     !  
  END IF
  !
#if defined (__PARA)
  CALL reduce( 2 * nstart * nstart, sc )
#endif
  !
  ! ... Diagonalize
  !
  CALL cdiaghg( nstart, nbnd, hc, sc, nstart, en, vc )
  !
  e(:) = en(1:nbnd)
  !
  ! ...  update the basis set
  !  
  CALL ZGEMM( 'N', 'N', npw, nbnd, nstart, ( 1.D0, 0.D0 ), psi, npwx, &
              vc, nstart, ( 0.D0, 0.D0 ), aux, npwx ) 
  !     
  evc(:,:) = aux(:,1:nbnd)
  !
  DEALLOCATE( en )
  DEALLOCATE( vc )
  DEALLOCATE( sc )
  DEALLOCATE( hc )
  DEALLOCATE( aux )
  !
  CALL stop_clock( 'wfcrot' )
  !
  RETURN
  !
END SUBROUTINE rotate_wfc

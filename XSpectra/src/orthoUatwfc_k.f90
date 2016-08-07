!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  this routine initializes some variables for DFT+U
!  calculation, in particular atomic wfc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



!-----------------------------------------------------------------------
SUBROUTINE orthoUwfc_k(ik)
  !-----------------------------------------------------------------------
  !
  ! This routine is meant to orthogonalize all the atomic wfcs. This is
  ! useful when we want to compute the occupation of the atomic orbitals
  ! in order to make lda+U calculations
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunhub, nwordwfcU, diropn
  USE ions_base,  ONLY : nat
  USE basis,      ONLY : natomwfc, swfcatom
  USE klist,      ONLY : nks, xk, ngk, igk_k
  USE ldaU,       ONLY : U_projection, wfcU, nwfcU, copy_U_wfc
  USE wvfct,      ONLY : npwx
  USE uspp,       ONLY : nkb, vkb
  USE becmod,     ONLY : allocate_bec_type, deallocate_bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol

  ! 
  IMPLICIT NONE
  !
  INTEGER :: ik, ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
       l, lm, ltot, ntot, ipol, npw
  ! the k point under consideration
  ! counter on bands
  LOGICAL :: orthogonalize_wfc, normalize_only
  COMPLEX(DP) , ALLOCATABLE :: wfcatom (:,:)
  LOGICAL :: exst

  ALLOCATE (wfcatom( npwx*npol, natomwfc),swfcatom(npwx*npol, natomwfc))    

  IF ( U_projection == "pseudo" ) THEN
     WRITE(stdout,'(8x,a)') &
       '! Beta functions used for LDA+U Projector'
     RETURN
  ELSE IF (U_projection=="file") THEN
     WRITE(stdout,'(8x,a)') &
       '! LDA+U Projector read from file '
     RETURN
  ELSE IF (U_projection=="atomic") THEN
     orthogonalize_wfc = .FALSE.
     normalize_only = .FALSE.
     WRITE(stdout,'(8x,a)') &
       '! Atomic wfc used for LDA+U Projector are NOT orthogonalized'
  ELSE IF (U_projection=="ortho-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .FALSE.    
     WRITE(stdout,'(8x,a)') &
        '! Atomic wfc used for LDA+U Projector are orthogonalized'
     IF (gamma_only) CALL errore('orthoatwfc', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE IF (U_projection=="norm-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .TRUE.
     WRITE(stdout,'(8x,a)') &
       '! Atomic wfc used for LDA+U Projector are normalized but NOT orthogonalized'
     IF (gamma_only) CALL errore('orthoatwfc', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE
     WRITE(stdout,'(8x,a)') & 
        "! U_projection_type =", U_projection
     CALL errore ("orthoatwfc"," this U_projection_type is not valid",1)
  END IF

  CALL allocate_bec_type (nkb,natomwfc, becp) 
    ! Allocate the array becp = <beta|wfcatom>
  CALL atomic_wfc (ik, wfcatom)
  npw = ngk(ik)
  CALL init_us_2(npw, igk_k(1,ik), xk (1, ik), vkb)
  CALL calbec (npw, vkb, wfcatom, becp)
  CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)

  IF (orthogonalize_wfc) &
     CALL ortho_swfc ( normalize_only, natomwfc, wfcatom, swfcatom )
  !
  CALL copy_U_wfc (swfcatom)
  !
  DEALLOCATE (wfcatom)
  DEALLOCATE (swfcatom)
  CALL deallocate_bec_type (becp )
  !
  RETURN
     
END SUBROUTINE orthoUwfc_k



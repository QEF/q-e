!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE lr_orthoUwfc (lflag)
  !-----------------------------------------------------------------------
  !
  ! This routine computes and writes to file atomic wavefunctions phi 
  ! at k and k+q. phi(k) is defined in Eq. (A6) in Ref. [1]
  ! [1] Phys. Rev. B 98, 085127 (2018)
  !
  ! lflag=.TRUE.         write phi(k) and phi(k+q) to file (unit iuatwfc)
  ! lflag=.FALSE. DO NOT write phi(k) and phi(k+q) to file (unit iuatwfc)
  !
  ! Independently from the value of the lflag, this routine will always
  ! write S(k)*phi(k) and S(k+q)*phi(k+q) to file with unit iuatswfc
  ! (note that this is not the same unit as iuatwfc).
  !
  ! If lgamma = .TRUE., write phi(k) and S(k)*phi(k) also to iunhub and
  ! iunhub_noS. These are needed in commutator_Vhubx_psi.
  !
  ! In the norm-conserving case, S(k)=1 and S(k+q)=1.
  ! Note: here the array wfcU is used as a workspace.
  ! Note2: wfcU is not deallocated because it is used also in ch_psi_all 
  ! Inspired by PW/src/orthoatwfc.f90
  !
  ! Written by I. Timrov (01.10.2018)
  !
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : iunhub, iunhub_noS, nwordwfcU
  USE ions_base,        ONLY : nat, ityp
  USE basis,            ONLY : natomwfc
  USE upf_ions,         ONLY : n_atom_wfc
  USE klist,            ONLY : xk, ngk, igk_k
  USE wvfct,            ONLY : npwx
  USE control_flags,    ONLY : gamma_only
  USE uspp,             ONLY : vkb, nkb, okvan
  USE becmod,           ONLY : allocate_bec_type_acc, deallocate_bec_type_acc,&
                               becp
  USE buffers,          ONLY : save_buffer
  USE io_global,        ONLY : stdout
  USE mp,               ONLY : mp_sum
  USE mp_pools,         ONLY : intra_pool_comm
  USE mp_bands,         ONLY : use_bgrp_in_hpsi
  USE noncollin_module, ONLY : noncolin, npol
  USE qpoint,           ONLY : nksq, ikks, ikqs
  USE control_lr,       ONLY : lgamma
  USE units_lr,         ONLY : iuatwfc, iuatswfc
  USE ldaU,             ONLY : Hubbard_projectors, wfcU, nwfcU, copy_U_wfc
  ! 
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: lflag
  !
  INTEGER :: ik,  & ! dummy k point index
             ikk, & ! k point
             ikq, & ! k+q point
             npw, & ! number of plane waves at k
             npwq   ! number of plane waves at k+q
  LOGICAL :: orthogonalize_wfc, normalize_only, save_flag
  COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:), &  ! atomic wfc
                              swfcatom(:,:)    ! S * atomic wfc
  !
  CALL start_clock ('lr_orthoUwfc')
  !
  IF (Hubbard_projectors=="atomic") THEN
     orthogonalize_wfc = .FALSE.
     normalize_only = .FALSE.
     WRITE( stdout, '(/5x,"Atomic wfc used for the projector on the Hubbard manifold are NOT orthogonalized")')
  ELSEIF (Hubbard_projectors=="ortho-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .FALSE.
     WRITE( stdout, '(/5x,"Atomic wfc used for the projector on the Hubbard manifold are orthogonalized")')
     IF (gamma_only) CALL errore('lr_orthoUwfc', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSEIF (Hubbard_projectors=="norm-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .TRUE.
     WRITE( stdout, '(/5x,"Atomic wfc used for the projector on the Hubbard manifold are normalized but NOT orthogonalized")')
     IF (gamma_only) CALL errore('lr_orthoUwfc', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE
     WRITE(stdout,*) "Hubbard_projectors =", Hubbard_projectors
     CALL errore ("hp_sphi"," This Hubbard projectors type is not valid",1)
  ENDIF
  !
  natomwfc = n_atom_wfc( nat, ityp, noncolin )
  ALLOCATE (wfcatom(npwx*npol,natomwfc))
  ALLOCATE (swfcatom(npwx*npol,natomwfc))
  !$acc enter data create(wfcatom, swfcatom)
  IF ( .NOT. ALLOCATED(wfcU) ) THEN
     ALLOCATE (wfcU(npwx*npol,nwfcU))
     !$acc enter data create(wfcU)
  ENDIF
  !
  save_flag = use_bgrp_in_hpsi ; use_bgrp_in_hpsi=.false.
  !
  ! Allocate the array becp = <beta|wfcatom>
  IF (okvan) CALL allocate_bec_type_acc (nkb,natomwfc,becp)
  !
  DO ik = 1, nksq
     !
     ikk  = ikks(ik)
     ikq  = ikqs(ik)
     npw  = ngk(ikk)
     npwq = ngk(ikq)
     !
     ! Determine the atomic orbital at k : phi(k)
     !
     IF (noncolin) THEN
       CALL atomic_wfc_nc_updown (ikk, wfcatom)
     ELSE
       CALL atomic_wfc (ikk, wfcatom)
     ENDIF
     !
     ! Compute S(k)*phi(k) (phi means the atomic orbital)
     !
     CALL s_phi (npw, ikk, wfcatom, swfcatom)     
     !
     ! Orthonormalize or normalize the atomic orbitals (if needed) 
     !
     IF (orthogonalize_wfc) THEN
        CALL ortho_swfc (npw, normalize_only, natomwfc, wfcatom, swfcatom, lflag)
     END IF
     !
     ! If lflag=.TRUE. copy the result from (orthonormalized) wfcatom 
     ! (which uses the offset oatwfc) to wfcU (which uses the offset offsetU),
     ! and then write wfcU = phi(k) to file with unit iuatwfc.
     !
     IF (lflag) THEN
        CALL copy_U_wfc (wfcatom, noncolin)
        !$acc update host(wfcU)
        CALL save_buffer (wfcU, nwordwfcU, iuatwfc, ikk)
        IF (lgamma) CALL save_buffer (wfcU, nwordwfcU, iunhub_noS, ik)
     ENDIF
     !
     ! Copy the result from (orthonormalized) swfcatom 
     ! (which uses the offset oatwfc) to wfcU (which uses the offset offsetU),
     ! and then write wfcU = S(k)*phi(k) to file with unit iuatswfc.
     !
     CALL copy_U_wfc (swfcatom, noncolin)
     !$acc update host(wfcU)
     CALL save_buffer (wfcU, nwordwfcU, iuatswfc, ikk)
     IF (lgamma) CALL save_buffer (wfcU, nwordwfcU, iunhub, ik)
     !
     IF (.NOT.lgamma) THEN
        !
        ! Determine the atomic orbital at k+q : phi(k+q)
        !
        IF (noncolin) THEN
           CALL atomic_wfc_nc_updown (ikq, wfcatom)
        ELSE 
           CALL atomic_wfc (ikq, wfcatom)
        ENDIF
        !
        ! Compute S(k+q)*phi(k+q) 
        !
        CALL s_phi (npwq, ikq, wfcatom, swfcatom)
        !
        ! Orthonormalize or normalize the atomic orbitals (if needed)
        !
        IF (orthogonalize_wfc) THEN 
           CALL ortho_swfc (npwq, normalize_only, natomwfc, wfcatom, swfcatom, lflag)
        END IF
        !
        ! If lflag=.TRUE. copy the result from (orthonormalized) wfcatom 
        ! (which uses the offset oatwfc) to wfcU (which uses the offset offsetU),
        ! and then write wfcU = phi(k) to file with unit iuatwfc.
        !
        IF (lflag) THEN
           CALL copy_U_wfc (wfcatom, noncolin)
           !$acc update host(wfcU)
           CALL save_buffer (wfcU, nwordwfcU, iuatwfc, ikq)
        ENDIF
        !
        ! Copy the result from (orthonormalized) swfcatom 
        ! (which uses the offset oatwfc) to wfcU (which uses the offset offsetU),
        ! and then write wfcU = S(k+q)*phi(k+q) to file with unit iuatswfc.
        !
        CALL copy_U_wfc (swfcatom, noncolin)
        !$acc update host(wfcU)
        CALL save_buffer (wfcU, nwordwfcU, iuatswfc, ikq)
        !
     ENDIF
     !
  ENDDO  
  !$acc exit data delete(wfcatom, swfcatom)
  DEALLOCATE (wfcatom)
  DEALLOCATE (swfcatom)
  IF (okvan) CALL deallocate_bec_type_acc (becp)
  !
  use_bgrp_in_hpsi = save_flag
  !
  CALL stop_clock ('lr_orthoUwfc')
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE s_phi (npw_, ik_, wfc, swfc)
  !-----------------------------------------------------------------------
  !
  ! NCPP: swfc = wfc
  ! USPP: swfc = S * wfc
  !
  USE kinds,          ONLY : DP
  USE becmod,         ONLY : calbec, becp
  USE uspp_init,      ONLY : init_us_2
  USE control_flags,  ONLY : offload_type
  !
  IMPLICIT NONE
  INTEGER,     INTENT(IN)  :: npw_
  INTEGER,     INTENT(IN)  :: ik_
  COMPLEX(DP), INTENT(IN)  :: wfc  (npwx*npol, natomwfc)
  COMPLEX(DP), INTENT(OUT) :: swfc (npwx*npol, natomwfc)
  !
  ! NCPP case
  !
  IF ( nkb == 0 .OR. .NOT. okvan ) THEN
     !$acc kernels
     swfc = wfc
     !$acc end kernels
     RETURN
  ENDIF
  !
  ! USPP case
  !
  ! Compute beta functions vkb at ik_
  !
  CALL init_us_2 (npw_, igk_k(1,ik_), xk(1,ik_), vkb, .true.)
  !
  ! Compute the product of beta functions vkb
  ! with the functions wfc : becp = <vkb|wfc>
  !
  CALL calbec (offload_type, npw_, vkb, wfc, becp)
  !
  ! Calculate S*|wfc> = |wfc> + \sum qq * |vkb> * becp 
  !  
  CALL s_psi_acc (npwx, npw_, natomwfc, wfc, swfc)
  !
  RETURN
  !
END SUBROUTINE s_phi
  
END SUBROUTINE lr_orthoUwfc

!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
!-----------------------------------------------------------------------
subroutine bcast_wfc (igk_k_all, ngk_all) 
  !-----------------------------------------------------------------------
  !
  !! This routine gathers all the wavefunctions in different pools 
  !! and save these on a buffer containing all the k-points (not just 
  !! the ones belonging to the particular pool) 
  !   
  !
  USE kinds,             ONLY : DP
  USE lsda_mod,          ONLY : nspin, lsda, isk
  USE klist,             ONLY : nkstot, nks, igk_k, ngk
  USE buffers,           ONLY : open_buffer, save_buffer,  get_buffer
  USE control_kcw,       ONLY : evc0, iuwfc_wann, iuwfc_wann_allk, num_wann, spin_component
  USE wvfct,             ONLY : npwx
  USE mp,                ONLY : mp_sum
  USE mp_pools,          ONLY : inter_pool_comm
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), ALLOCATABLE :: evc0_aux(:,:,:) 
  ! ... a vector to store all the wfc at all the k-points 
  !
  INTEGER :: igk_k_all(npwx,nkstot)
  ! ... index of G corresponding to a given index of k+G
  !
  INTEGER :: ngk_all(nkstot)
  ! ... number of plane waves for each k point
  !
  INTEGER :: global_ik, ik
  ! ... global and local index of the k point 
  !
  INTEGER   :: lrwfc, ik_eff
  ! ... record length , counter 
  !
  INTEGER, EXTERNAL :: global_kpoint_index
  ! ... find the global kpoint index given the local one. 
  !
  ALLOCATE ( evc0_aux(npwx, num_wann, nkstot/nspin) )
  lrwfc = num_wann * npwx
  evc0_aux = CMPLX(0.0, 0.0, kind=DP)
  igk_k_all = 0
  ngk_all = 0
  !
  DO ik = 1, nks 
     IF (lsda .AND. isk(ik) /= spin_component) CYCLE
     global_ik = global_kpoint_index (nkstot, ik)
     global_ik = global_ik - (spin_component-1)*nkstot/nspin
     ik_eff = ik - (spin_component-1)*nkstot/nspin
     igk_k_all(:,global_ik) = igk_k(:,ik)
     ngk_all(global_ik) = ngk(ik)
     CALL get_buffer ( evc0, lrwfc, iuwfc_wann, ik_eff )
     evc0_aux(:,:,global_ik) = evc0(:,:)
  ENDDO
  !
  ! ... mp_sum used as a trick to gather all the k points
  CALL mp_sum (evc0_aux, inter_pool_comm)
  CALL mp_sum (ngk_all, inter_pool_comm)
  CALL mp_sum (igk_k_all, inter_pool_comm)
  !
  DO ik = 1, nkstot/nspin
    !ik_eff = ik + (spin_component-1)*nkstot/nspin
    !WRITE(*,'("NICOLA DOPO", 3i5, 6f20.15)'), my_pool_id, ik, ik, igk_k_all (1:3,ik)
    lrwfc = num_wann * npwx
    evc0(:,:) = evc0_aux (:,:,ik)
    CALL save_buffer (evc0, lrwfc, iuwfc_wann_allk, ik)
  ENDDO
  DEALLOCATE (evc0_aux) 
  !
  RETURN
  !
END SUBROUTINE bcast_wfc

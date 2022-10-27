!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#include "f_defs.h"
#define DEBUG
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!-----------------------------------------------------------------------
SUBROUTINE rotate_ks ()
  !-----------------------------------------------------------------------
  !
  !! This routine rotate the KS orbital to the localized representation
  !
  USE wavefunctions,         ONLY : evc
  USE kinds,                 ONLY : DP
  USE io_global,             ONLY : stdout, ionode
  USE io_files,              ONLY : nwordwfc
  USE units_lr,              ONLY : iuwfc
  USE wvfct,                 ONLY : npw, nbnd, npwx, current_k
  USE lsda_mod,              ONLY : lsda, isk, nspin, current_spin
  USE klist,                 ONLY : ngk, xk, igk_k, nkstot, nks
  USE uspp,                  ONLY : nkb, vkb
  USE uspp_init,             ONLY : init_us_2
  USE buffers,               ONLY : get_buffer, save_buffer
  USE control_kcw,           ONLY : kcw_at_ks, evc0, read_unitary_matrix, iuwfc_wann, &
                                    check_ks, spin_component, num_wann, occ_mat
  USE control_lr,            ONLY : nbnd_occ
  !
  USE wvfct,                 ONLY : wg
  USE klist,                 ONLY : wk
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, n_orb, ik_eff
  INTEGER :: lrwfc 
  !
  COMPLEX(DP),EXTERNAL :: zdotc
  INTEGER, EXTERNAL :: global_kpoint_index
  INTEGER :: global_ik
  INTEGER :: i
  REAL(DP) :: occ_mat_aux(num_wann,num_wann)
  !
  IF ( ionode )  THEN 
    IF (read_unitary_matrix) WRITE( stdout, '(/,5x,A)') &
               'INFO: Minimizing orbitals from Unitary Matrix Rotation'
    IF (kcw_at_ks) WRITE( stdout, '(/,5x,A,/,5x,A)') &
               'INFO: NO ROTATION, using the canonical KS manifold'
    IF (.not. read_unitary_matrix .and. .not. kcw_at_ks) & 
    CALL errore('rotate_ks', 'Specify how to get the minimizing orbital. From file or from unitary matrix?',1)
  ENDIF
  !
  !
  IF (check_KS .AND. .NOT. kcw_at_ks ) & 
      WRITE( stdout, '(/,8x,A)') &
               'INFO: Performing a check on the eigenvalues of the rotated KS Hamilotnian ... '
  !
  ! ... Loop over k_point
  !
  k_loop: DO ik = 1, nks
     !
     current_k = ik
     IF ( lsda ) current_spin = isk(ik)
     !
     ! ... the global kpoint index 
     global_ik = global_kpoint_index(nkstot, ik)
     IF ( lsda .AND. isk(ik) /= spin_component) CYCLE
     !
     npw = ngk(ik)
     !
     CALL get_buffer ( evc, nwordwfc, iuwfc, ik ) 
     !
     IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
     !
     evc0(:,:) = ZERO
     ik_eff = global_ik - (isk(ik)-1)*(nkstot/nspin)
     !
     IF (kcw_at_ks) THEN 
        !
        !  ... If ki@ks simply store evc in evc0 ...
        evc0 = evc 
        n_orb = nbnd
        DO i = 1, num_wann; occ_mat_aux(i,i) = wg(i,ik_eff)/wk(ik_eff); ENDDO
        !
     ELSE
        ! 
        !  ... If KI@WANN rotate the KS orbital
        n_orb = nbnd_occ(ik)
        !
        !  ... If unitary matrix provided, rotate the KS==>MLWF and store in evc0 ...
        ik_eff = global_ik - (isk(ik)-1)*(nkstot/nspin)
        !! ... In PW the spin down are treatead as k points. In W90 I have only up or dw
        !! ... and ik run always from 1 to nkstot/nspin (see read_wannier.f90). 
        !
        occ_mat_aux = 0.D0
        CALL apply_u_matrix(evc, evc0, occ_mat_aux, ik_eff,n_orb)
        !
     ENDIF
     ! 
     ! ... Store the (Rotated) occupation matrix
     occ_mat(:,:,ik_eff) = occ_mat_aux(:,:)
     !
     ! ... Store the KS states in the Wannier gauge 
     !
     lrwfc = num_wann * npwx 
     ik_eff = ik-(spin_component-1)*nkstot/nspin
     CALL save_buffer ( evc0, lrwfc, iuwfc_wann, ik_eff )
     !
     !
     ! ... Check that the rotation did not spoil the KS eigenvalues
     ! ... and store it for later Hamiltonian diagonalization
     !
     CALL ks_hamiltonian(evc0, ik, n_orb) 
     !
  ENDDO k_loop
  !
  IF (check_ks .AND. .NOT. kcw_at_ks )  WRITE( stdout, '(/,8x,A)') &
               'INFO: Performing a check on the eigenvalues of the rotated KS Hamilotnian ... DONE'
  WRITE(stdout, '(/,5X,"INFO: Minimizing orbitals DEFINED")')
  !
END subroutine rotate_ks

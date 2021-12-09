!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE compute_becsum( iflag )
  !----------------------------------------------------------------------------
  !! Compute "becsum" = \sum_i w_i <psi_i|beta_l><beta_m|\psi_i> term.
  !! Output in module uspp and (PAW only) in rho%bec (symmetrized)
  !! if iflag = 1, weights w_k are re-computed.
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE scf,                  ONLY : rho
  USE uspp,                 ONLY : nkb, vkb, becsum, becsum_d, okvan
  USE uspp_param,           ONLY : nhm
  USE wavefunctions,        ONLY : evc
  USE wvfct,                ONLY : nbnd, npwx, wg
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm 
  USE mp,                   ONLY : mp_sum, mp_get_comm_null
  USE paw_symmetry,         ONLY : PAW_symmetrize
  USE paw_variables,        ONLY : okpaw
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   bec_type, becp
  USE wavefunctions_gpum,   ONLY : using_evc
  USE becmod_subs_gpum,     ONLY : using_becp_auto
  USE uspp_init,            ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iflag
  !
  INTEGER :: ik, & ! counter on k points
             ibnd_start, ibnd_end, this_bgrp_nbnd ! first, last and number of band in this bgrp
  !
  !
  IF ( .NOT. okvan ) RETURN
  CALL using_evc(0)
  !
  CALL start_clock( 'compute_becsum' )
  !
  ! ... calculates weights of Kohn-Sham orbitals
  !
  IF ( iflag == 1) CALL weights( )
  !
  becsum(:,:,:) = 0.D0
  CALL allocate_bec_type( nkb,nbnd, becp,intra_bgrp_comm )
  CALL using_becp_auto(2)
  CALL divide( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
  this_bgrp_nbnd = ibnd_end - ibnd_start + 1
  !
  k_loop: DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     IF ( nks > 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, ik )
     IF ( nks > 1 ) CALL using_evc(2)
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
     !
     ! ... actual calculation is performed inside routine "sum_bec"
     !
     CALL sum_bec( ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd )
     !
  ENDDO k_loop
  !
  ! ... with distributed <beta|psi>, sum over bands
  !
  IF( becp%comm /= mp_get_comm_null() ) &
       CALL mp_sum( becsum, becp%comm )
  !
#if defined __CUDA
  ! GPU-sync (protecting zero-size allocations)
  if (nhm>0) becsum_d=becsum
#endif

  !
  ! ... Needed for PAW: becsum has to be symmetrized so that they reflect 
  ! ... a real integral in k-space, not only on the irreducible zone. 
  ! ... For USPP there is no need to do this as becsums are only used
  ! ... to compute the density, which is symmetrized later.
  !
  IF ( okpaw )  THEN
     rho%bec(:,:,:) = becsum(:,:,:) ! becsum is filled in sum_band_{k|gamma}
     CALL mp_sum( rho%bec, inter_pool_comm )
     call mp_sum( rho%bec, inter_bgrp_comm )
     CALL PAW_symmetrize( rho%bec )
  ENDIF
  !
  CALL deallocate_bec_type( becp )
  CALL using_becp_auto(2)
  !
  CALL stop_clock( 'compute_becsum' )
  !
END SUBROUTINE compute_becsum

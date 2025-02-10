!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine bcast_kcw_input ( )
  !-----------------------------------------------------------------------
  !
  !! This routine broadcast all the KC input variables 
  !
#ifdef __MPI
  !
  USE run_info,         ONLY : title
  USE io_global,        ONLY : ionode_id
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : intra_image_comm
  USE io_files,         ONLY : tmp_dir, prefix
  USE control_kcw
  USE control_lr,       ONLY : lrpa
  USE input_parameters, ONLY : assume_isolated
                            
  !
  IMPLICIT NONE
  !
  call mp_bcast ( title,               ionode_id, intra_image_comm ) 
  call mp_bcast ( tmp_dir,             ionode_id, intra_image_comm )
  call mp_bcast ( prefix,              ionode_id, intra_image_comm )
  call mp_bcast ( kcw_at_ks,           ionode_id, intra_image_comm )
  call mp_bcast ( fix_orb,             ionode_id, intra_image_comm )
  call mp_bcast ( homo_only,           ionode_id, intra_image_comm )
  call mp_bcast ( spread_thr,          ionode_id, intra_image_comm )
  call mp_bcast ( read_unitary_matrix, ionode_id, intra_image_comm )
  call mp_bcast ( qp_symm,             ionode_id, intra_image_comm )
  call mp_bcast ( kipz_corr,           ionode_id, intra_image_comm )
  call mp_bcast ( has_disentangle,     ionode_id, intra_image_comm )
  call mp_bcast ( have_empty,          ionode_id, intra_image_comm )
  call mp_bcast ( seedname,            ionode_id, intra_image_comm )
  call mp_bcast ( num_wann_occ,        ionode_id, intra_image_comm )
  call mp_bcast ( num_wann_emp,        ionode_id, intra_image_comm )
  call mp_bcast ( check_ks,            ionode_id, intra_image_comm )
  call mp_bcast ( kcw_iverbosity,      ionode_id, intra_image_comm )
  call mp_bcast ( spin_component,      ionode_id, intra_image_comm )
  call mp_bcast ( niter,               ionode_id, intra_image_comm )
  call mp_bcast ( alpha_mix,           ionode_id, intra_image_comm )
  call mp_bcast ( nmix,                ionode_id, intra_image_comm )
  call mp_bcast ( tr2,                 ionode_id, intra_image_comm )
  call mp_bcast ( lrpa,                ionode_id, intra_image_comm )
  call mp_bcast ( i_orb,               ionode_id, intra_image_comm )
  call mp_bcast ( mp1,                 ionode_id, intra_image_comm )
  call mp_bcast ( mp2,                 ionode_id, intra_image_comm )
  call mp_bcast ( mp3,                 ionode_id, intra_image_comm )
  call mp_bcast ( do_bands,            ionode_id, intra_image_comm )
  call mp_bcast ( use_ws_distance,     ionode_id, intra_image_comm )
  call mp_bcast ( write_hr,            ionode_id, intra_image_comm )
  call mp_bcast ( l_vcut,              ionode_id, intra_image_comm )
  call mp_bcast ( assume_isolated,     ionode_id, intra_image_comm )
  call mp_bcast ( eps_inf,             ionode_id, intra_image_comm )
  call mp_bcast ( l_alpha_corr,        ionode_id, intra_image_comm )
  call mp_bcast ( l_unique_manifold,   ionode_id, intra_image_comm )
  call mp_bcast ( check_spread,        ionode_id, intra_image_comm )
  call mp_bcast ( on_site_only,        ionode_id, intra_image_comm )
  call mp_bcast ( calculation,         ionode_id, intra_image_comm )
  call mp_bcast ( io_sp,               ionode_id, intra_image_comm )
  call mp_bcast ( io_real_space,       ionode_id, intra_image_comm )
  call mp_bcast ( irr_bz,              ionode_id, intra_image_comm )
  call mp_bcast ( use_wct,             ionode_id, intra_image_comm )
   !
#endif
  !
  return
  !
end subroutine bcast_kcw_input

!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE bcast_lr_input
  !-----------------------------------------------------------------------
  !
  !  ...The first processor sends the input to all the other processors.
  !
#if defined(__MPI)

  USE lr_variables
  USE lr_dav_variables
  USE realus,              ONLY: real_space, real_space_debug
  USE mp,                  ONLY: mp_bcast, mp_barrier
  USE io_files,            ONLY: tmp_dir, prefix, wfc_dir
  USE control_flags,       ONLY: tqr, tddfpt
  USE charg_resp,          ONLY: omeg, w_T_prefix, w_T_npol,epsil
  USE io_global,           ONLY: ionode, ionode_id
  USE mp_global,           ONLY: intra_image_comm
  USE mp_world,            ONLY: world_comm
  USE exx,                 ONLY: ecutfock
  USE qpoint,              ONLY: xq
  USE control_lr,          ONLY: lrpa

  IMPLICIT NONE
  !
  CALL mp_barrier(world_comm)
  CALL mp_bcast (lr_io_level, ionode_id, world_comm )
  CALL mp_bcast (itermax, ionode_id, world_comm )
  CALL mp_bcast (itermax_int, ionode_id, world_comm )
  CALL mp_bcast (charge_response, ionode_id, world_comm )
  CALL mp_bcast (sum_rule, ionode_id, world_comm )
  CALL mp_bcast (project, ionode_id, world_comm )
  CALL mp_bcast (restart, ionode_id, world_comm )
  CALL mp_bcast (restart_step, ionode_id, world_comm )
  CALL mp_bcast (lr_verbosity, ionode_id, world_comm )
  CALL mp_bcast (prefix, ionode_id, world_comm )
  CALL mp_bcast (tmp_dir, ionode_id, world_comm )
  CALL mp_bcast (wfc_dir, ionode_id, world_comm )
  CALL mp_bcast (LR_polarization, ionode_id, world_comm )
  CALL mp_bcast (ltammd, ionode_id, world_comm )
  CALL mp_bcast (pseudo_hermitian, ionode_id, world_comm )
  CALL mp_bcast (real_space, ionode_id, world_comm )
  CALL mp_bcast (real_space_debug, ionode_id, world_comm )
  CALL mp_bcast (tqr, ionode_id, world_comm )
  CALL mp_bcast (test_case_no, ionode_id, world_comm )
  CALL mp_bcast (omeg, ionode_id, world_comm )
  CALL mp_bcast (epsil, ionode_id, world_comm )
  CALL mp_bcast (w_T_prefix, ionode_id, world_comm )
  CALL mp_bcast (w_T_npol, ionode_id, world_comm )
  CALL mp_bcast (n_ipol, ionode_id, world_comm )
  CALL mp_bcast (plot_type, ionode_id, world_comm )
  CALL mp_bcast (no_hxc, ionode_id, world_comm )
  CALL mp_bcast (bgz_suffix, ionode_id, world_comm )
  call mp_bcast (scissor, ionode_id, world_comm)
  CALL mp_bcast (ecutfock, ionode_id, world_comm )
  CALL mp_bcast (d0psi_rs, ionode_id,world_comm )
  CALL mp_bcast (lshift_d0psi, ionode_id,world_comm )
  CALL mp_bcast (tddfpt, ionode_id, world_comm )
  CALL plugin_arguments_bcast(ionode_id, world_comm)

  ! for EELS
  call mp_bcast (eels, ionode_id, world_comm )
  call mp_bcast (q1, ionode_id, world_comm )
  call mp_bcast (q2, ionode_id, world_comm )
  call mp_bcast (q3, ionode_id, world_comm )
  call mp_bcast (xq, ionode_id, world_comm )
  call mp_bcast (approximation, ionode_id, world_comm ) 
  call mp_bcast (lrpa, ionode_id, world_comm ) 

  ! for lr_dav
  CALL mp_bcast (davidson, ionode_id, world_comm )
  CALL mp_bcast (num_eign, ionode_id, world_comm )
  CALL mp_bcast (num_init, ionode_id, world_comm )
  CALL mp_bcast (num_basis_max, ionode_id, world_comm )
  CALL mp_bcast (residue_conv_thr, ionode_id, world_comm )
  CALL mp_bcast (precondition, ionode_id, world_comm )
  CALL mp_bcast (dav_debug, ionode_id, world_comm )
  CALL mp_bcast (reference, ionode_id, world_comm )
  CALL mp_bcast (vccouple_shift, ionode_id, world_comm )
  CALL mp_bcast (single_pole, ionode_id, world_comm )
  CALL mp_bcast (sort_contr, ionode_id, world_comm )
  CALL mp_bcast (diag_of_h, ionode_id, world_comm )
  CALL mp_bcast (close_pre, ionode_id, world_comm )
  CALL mp_bcast (broadening, ionode_id, world_comm )
  CALL mp_bcast (print_spectrum, ionode_id, world_comm )
  CALL mp_bcast (start, ionode_id, world_comm )
  CALL mp_bcast (finish, ionode_id, world_comm )
  CALL mp_bcast (step, ionode_id, world_comm )
  CALL mp_bcast (if_check_orth, ionode_id, world_comm )
  CALL mp_bcast (if_random_init, ionode_id, world_comm )
  CALL mp_bcast (if_check_her, ionode_id, world_comm )
  CALL mp_bcast (p_nbnd_occ, ionode_id, world_comm )
  CALL mp_bcast (p_nbnd_virt, ionode_id, world_comm )
  CALL mp_bcast (poor_of_ram, ionode_id, world_comm )
  CALL mp_bcast (poor_of_ram2, ionode_id, world_comm )
  CALL mp_bcast (max_iter, ionode_id, world_comm )
  CALL mp_bcast (conv_assistant, ionode_id, world_comm )
  CALL mp_bcast (if_dft_spectrum, ionode_id, world_comm )
  CALL mp_bcast (lplot_drho, ionode_id, world_comm )
  CALL mp_barrier(world_comm)
  
#endif
  RETURN
END SUBROUTINE bcast_lr_input

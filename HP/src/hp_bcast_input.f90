!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE hp_bcast_input ( )
  !-----------------------------------------------------------------------
  !
  ! In this routine the first processor sends the input parameters to all
  ! the other processors
  !
#if defined (__MPI)

  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : world_comm
  USE io_files,         ONLY : tmp_dir, prefix
  USE control_flags,    ONLY : iverbosity
  USE input_parameters, ONLY : max_seconds
  USE io_global,        ONLY : meta_ionode_id
  USE control_lr,       ONLY : lrpa, ethr_nscf
  USE ldaU_hp,          ONLY : conv_thr_chi, thresh_init, find_atpert, skip_atom,      &
                               skip_type, equiv_type, alpha_mix, start_q, last_q,      &
                               background, compute_hp, sum_pertq, perturb_only_atom,   &
                               determine_num_pert_only, skip_equivalence_q, niter_max, &
                               disable_type_analysis, docc_thr, num_neigh, lmin, rmax, &
                               nmix, nq1, nq2, nq3
  !
  IMPLICIT NONE
  !
  ! Logicals
  !
  CALL mp_bcast (skip_atom, meta_ionode_id, world_comm)
  CALL mp_bcast (skip_type, meta_ionode_id, world_comm)
  CALL mp_bcast (perturb_only_atom, meta_ionode_id, world_comm)
  CALL mp_bcast (skip_equivalence_q, meta_ionode_id, world_comm)
  CALL mp_bcast (equiv_type, meta_ionode_id, world_comm)
  CALL mp_bcast (background, meta_ionode_id, world_comm)
  CALL mp_bcast (compute_hp, meta_ionode_id, world_comm)
  CALL mp_bcast (sum_pertq, meta_ionode_id, world_comm)
  CALL mp_bcast (lrpa, meta_ionode_id, world_comm)
  CALL mp_bcast (determine_num_pert_only, meta_ionode_id, world_comm)
  CALL mp_bcast (disable_type_analysis, meta_ionode_id, world_comm)
  !
  ! Integers
  !
  CALL mp_bcast (nq1, meta_ionode_id, world_comm)
  CALL mp_bcast (nq2, meta_ionode_id, world_comm)
  CALL mp_bcast (nq3, meta_ionode_id, world_comm)
  CALL mp_bcast (start_q, meta_ionode_id, world_comm)
  CALL mp_bcast (last_q, meta_ionode_id, world_comm)
  CALL mp_bcast (find_atpert, meta_ionode_id, world_comm)
  CALL mp_bcast (iverbosity, meta_ionode_id, world_comm)
  CALL mp_bcast (niter_max, meta_ionode_id, world_comm)
  CALL mp_bcast (nmix, meta_ionode_id, world_comm)
  CALL mp_bcast (num_neigh, meta_ionode_id, world_comm)
  CALL mp_bcast (lmin, meta_ionode_id, world_comm)
  !
  ! Real*8
  !
  CALL mp_bcast (conv_thr_chi, meta_ionode_id, world_comm)
  CALL mp_bcast (thresh_init, meta_ionode_id, world_comm)
  CALL mp_bcast (ethr_nscf, meta_ionode_id, world_comm)
  CALL mp_bcast (docc_thr, meta_ionode_id, world_comm)
  CALL mp_bcast (alpha_mix, meta_ionode_id, world_comm)
  CALL mp_bcast (max_seconds, meta_ionode_id, world_comm)
  CALL mp_bcast (rmax, meta_ionode_id, world_comm)
  !
  ! Characters
  !
  CALL mp_bcast (prefix, meta_ionode_id, world_comm)
  CALL mp_bcast (tmp_dir, meta_ionode_id, world_comm)
  !
#endif
  !
  RETURN
  !
END SUBROUTINE hp_bcast_input

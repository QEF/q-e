!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine bcast_lr_input
  !-----------------------------------------------------------------------
  !
  !     The first processor sends the input to all the other processors
  !
  !
  ! Modified by Osman Baris Malcioglu in 2009 
#ifdef __PARA
#include "f_defs.h"

  use lr_variables
  use realus,              only: real_space, real_space_debug
  use mp,                  only: mp_bcast, mp_barrier
  use io_files,            only: tmp_dir, prefix, wfc_dir
  USE control_flags,       ONLY: tqr
  USE charg_resp,          ONLY: omeg, w_T_prefix, w_T_npol,epsil
  USE io_global,           ONLY: ionode, ionode_id
  USE mp_global,           ONLY: intra_image_comm

  implicit none
  !
  ! 
  !
  call mp_barrier()
  call mp_bcast (itermax, ionode_id )
  call mp_bcast (itermax_int, ionode_id )
  call mp_bcast (charge_response, ionode_id )
  call mp_bcast (project, ionode_id )
  call mp_bcast (restart, ionode_id )
  call mp_bcast (restart_step, ionode_id )
  call mp_bcast (lr_verbosity, ionode_id )
  call mp_bcast (prefix, ionode_id )
  call mp_bcast (tmp_dir, ionode_id )
  call mp_bcast (wfc_dir, ionode_id )
  call mp_bcast (LR_polarization, ionode_id )
  call mp_bcast (ltammd, ionode_id )
  !call mp_bcast (broadening, ionode_id )
  call mp_bcast (real_space, ionode_id )
  call mp_bcast (real_space_debug, ionode_id )
  call mp_bcast (tmp_dir, ionode_id )
  call mp_bcast (tqr, ionode_id )
  call mp_bcast (test_case_no, ionode_id )
  call mp_bcast (omeg, ionode_id )
  call mp_bcast (epsil, ionode_id )
  call mp_bcast (w_T_prefix, ionode_id )
  call mp_bcast (w_T_npol, ionode_id )
  call mp_bcast (n_ipol, ionode_id )
  call mp_bcast (plot_type, ionode_id )
  call mp_bcast (no_hxc, ionode_id )
  call mp_bcast (bgz_suffix, ionode_id )
  call mp_barrier()
  !print *, "bcast lr input finished"
  !print *, "variables"
  !print *, "prefix=", prefix
  !print *, "tmp_dir", tmp_dir
  !print *, "w_T_prefix=", w_T_prefix
  !print *, "charge_response=",charge_response
  !print *, "omeg=",omeg
  !print *, "test_case_no=",test_case_no
#endif
  return
end subroutine bcast_lr_input

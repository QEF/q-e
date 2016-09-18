!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!
!-----------------------------------------------------------------------
subroutine bcast_ph_input ( )
  !-----------------------------------------------------------------------
  !
  !     In this routine the first processor sends the phonon input to all
  !     the other processors
  !
  !
#if defined(__MPI)

  use mp, only: mp_bcast
  use mp_world, only: world_comm
  USE mp_global, only : intra_image_comm
  USE control_ph, ONLY : start_irr, last_irr, start_q, last_q, nmix_ph, &
                         niter_ph, lnoloc, alpha_mix, tr2_ph, recover, &
                         ldisp,  reduce_io, zue, zeu, epsil, trans
  USE gamma_gamma, ONLY : asr
  USE disp, ONLY : nq1, nq2, nq3
  USE partial, ONLY : nat_todo
  USE freq_ph, ONLY : fpol
  USE output, ONLY : fildvscf, fildyn, fildrho
  use io_files, ONLY : tmp_dir, prefix
  USE control_flags, only: iverbosity, modenum
  USE ramanm, ONLY: lraman, elop, dek, eth_rps, eth_ns
  USE input_parameters, ONLY: max_seconds
  USE ions_base,     ONLY : amass
  USE io_global, ONLY : ionode_id
  USE run_info, ONLY : title
  USE wannier_gw,    ONLY : l_head, omega_gauss, n_gauss, grid_type, nsteps_lanczos,&
                             &second_grid_n,second_grid_i,l_scissor,scissor,len_head_block_freq,len_head_block_wfc

  USE control_lr, ONLY : lgamma, lrpa

  implicit none
  !
  ! logicals
  !
  call mp_bcast (lgamma, ionode_id, world_comm )
  call mp_bcast (epsil, ionode_id, world_comm )
  call mp_bcast (trans, ionode_id, world_comm )
  call mp_bcast (zue, ionode_id, world_comm )
  call mp_bcast (zeu, ionode_id, world_comm )
  call mp_bcast (reduce_io, ionode_id, world_comm )
  call mp_bcast (ldisp, ionode_id, world_comm )
  call mp_bcast (lraman, ionode_id, world_comm )
  call mp_bcast (elop, ionode_id, world_comm )
  call mp_bcast (fpol, ionode_id, world_comm )
  call mp_bcast (recover, ionode_id, world_comm )
  call mp_bcast (asr, ionode_id, world_comm )
  call mp_bcast (lrpa, ionode_id, world_comm )
  call mp_bcast (lnoloc, ionode_id, world_comm )
  !
  ! integers
  !
  call mp_bcast (start_irr, ionode_id, world_comm )
  call mp_bcast (last_irr, ionode_id, world_comm )
  call mp_bcast (start_q, ionode_id, world_comm )
  call mp_bcast (last_q, ionode_id, world_comm )
  call mp_bcast (niter_ph, ionode_id, world_comm )
  call mp_bcast (nmix_ph, ionode_id, world_comm )
  call mp_bcast (iverbosity, ionode_id, world_comm )
  call mp_bcast (modenum, ionode_id, world_comm )
  call mp_bcast (nat_todo, ionode_id, world_comm )
  CALL mp_bcast( nq1, ionode_id, world_comm )
  CALL mp_bcast( nq2, ionode_id, world_comm )
  CALL mp_bcast( nq3, ionode_id, world_comm )
  !
  ! real*8
  !
  call mp_bcast (tr2_ph, ionode_id, world_comm )
  call mp_bcast (eth_rps, ionode_id, world_comm )
  call mp_bcast (eth_ns, ionode_id, world_comm )
  call mp_bcast (amass, ionode_id, world_comm )
  call mp_bcast (alpha_mix, ionode_id, world_comm )
  call mp_bcast (max_seconds, ionode_id, world_comm )
  call mp_bcast (dek, ionode_id, world_comm )
  !
  ! characters
  !
  call mp_bcast (title, ionode_id, world_comm )
  call mp_bcast (fildyn, ionode_id, world_comm )
  call mp_bcast (fildvscf, ionode_id, world_comm )
  call mp_bcast (fildrho, ionode_id, world_comm )
  call mp_bcast (tmp_dir, ionode_id, world_comm )
  call mp_bcast (prefix, ionode_id, world_comm )
   !
! head and wings
!
  call mp_bcast(l_head, ionode_id, world_comm)
  call mp_bcast(omega_gauss, ionode_id, world_comm)
  call mp_bcast(n_gauss, ionode_id, world_comm)
  call mp_bcast(grid_type, ionode_id, world_comm)
  call mp_bcast(nsteps_lanczos, ionode_id, world_comm)
  call mp_bcast(second_grid_n, ionode_id, world_comm)
  call mp_bcast(second_grid_i, ionode_id, world_comm)
  call mp_bcast(l_scissor, ionode_id, world_comm)
  call mp_bcast(scissor, ionode_id, world_comm)
  call mp_bcast(len_head_block_freq, ionode_id, world_comm)
  call mp_bcast(len_head_block_wfc, ionode_id,world_comm)

#endif
  return
end subroutine bcast_ph_input

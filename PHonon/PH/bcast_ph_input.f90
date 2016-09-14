!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
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
  USE control_lr, ONLY : lgamma, lrpa
  USE control_ph, ONLY : start_irr, last_irr, start_q, last_q, nmix_ph, &
                         niter_ph, lnoloc, alpha_mix, tr2_ph, recover, &
                         ldisp, reduce_io, zue, zeu, epsil, trans, &
                         ldiag, lqdir, search_sym,  electron_phonon, &
                         qplot, only_init, only_wfc, low_directory_check
  USE gamma_gamma, ONLY : asr
  USE disp, ONLY : nq1, nq2, nq3
  USE partial, ONLY : nat_todo
  USE freq_ph, ONLY : fpol
  USE output, ONLY : fildvscf, fildyn, fildrho
  use io_files, ONLY : tmp_dir, prefix
  USE control_flags, only: iverbosity, modenum
  USE ramanm, ONLY: lraman, elop, dek, eth_rps, eth_ns
  USE input_parameters, ONLY: max_seconds
  USE input_parameters, ONLY : nk1, nk2, nk3, k1, k2, k3
  USE ions_base,     ONLY : amass
  USE io_global,   ONLY : meta_ionode_id
  USE run_info,   ONLY : title
  USE el_phon,    ONLY : elph_nbnd_min,elph_nbnd_max,el_ph_ngauss, el_ph_nsigma, el_ph_sigma
  USE dfile_star, ONLY : drho_star, dvscf_star
  ! YAMBO >
  USE YAMBO,      ONLY : elph_yambo,dvscf_yambo
  ! YAMBO <

  implicit none
  !
  ! logicals
  !
  call mp_bcast (lgamma, meta_ionode_id, world_comm )
  call mp_bcast (epsil, meta_ionode_id, world_comm )
  call mp_bcast (trans, meta_ionode_id, world_comm )
  call mp_bcast (zue, meta_ionode_id, world_comm )
  call mp_bcast (zeu, meta_ionode_id, world_comm )
  call mp_bcast (reduce_io, meta_ionode_id, world_comm )
  call mp_bcast (ldisp, meta_ionode_id, world_comm )
  call mp_bcast (lraman, meta_ionode_id, world_comm )
  call mp_bcast (elop, meta_ionode_id, world_comm )
  call mp_bcast (fpol, meta_ionode_id, world_comm )
  call mp_bcast (recover, meta_ionode_id, world_comm )
  call mp_bcast (asr, meta_ionode_id, world_comm )
  call mp_bcast (lrpa, meta_ionode_id, world_comm )
  call mp_bcast (lnoloc, meta_ionode_id, world_comm )
  call mp_bcast (ldiag, meta_ionode_id, world_comm )
  call mp_bcast (lqdir, meta_ionode_id, world_comm )
  call mp_bcast (qplot, meta_ionode_id, world_comm )
  call mp_bcast (only_wfc, meta_ionode_id, world_comm )
  call mp_bcast (only_init, meta_ionode_id, world_comm )
  call mp_bcast (search_sym, meta_ionode_id, world_comm)
  ! YAMBO >
  call mp_bcast (elph_yambo, meta_ionode_id, world_comm)
  call mp_bcast (dvscf_yambo, meta_ionode_id, world_comm)
  ! YAMBO <
  !
  ! integers
  !
  call mp_bcast (start_irr, meta_ionode_id, world_comm )
  call mp_bcast (last_irr, meta_ionode_id, world_comm )
  call mp_bcast (start_q, meta_ionode_id, world_comm )
  call mp_bcast (last_q, meta_ionode_id, world_comm )
  call mp_bcast (niter_ph, meta_ionode_id, world_comm )
  call mp_bcast (nmix_ph, meta_ionode_id, world_comm )
  call mp_bcast (iverbosity, meta_ionode_id, world_comm )
  call mp_bcast (modenum, meta_ionode_id, world_comm )
  call mp_bcast (nat_todo, meta_ionode_id, world_comm )
  CALL mp_bcast( nq1, meta_ionode_id, world_comm )
  CALL mp_bcast( nq2, meta_ionode_id, world_comm )
  CALL mp_bcast( nq3, meta_ionode_id, world_comm )
  CALL mp_bcast( nk1, meta_ionode_id, world_comm )
  CALL mp_bcast( nk2, meta_ionode_id, world_comm )
  CALL mp_bcast( nk3, meta_ionode_id, world_comm )
  CALL mp_bcast( k1, meta_ionode_id, world_comm )
  CALL mp_bcast( k2, meta_ionode_id, world_comm )
  CALL mp_bcast( k3, meta_ionode_id, world_comm )
  CALL mp_bcast( low_directory_check, meta_ionode_id, world_comm )
  CALL mp_bcast( elph_nbnd_min, meta_ionode_id, world_comm )
  CALL mp_bcast( elph_nbnd_max, meta_ionode_id, world_comm )
  CALL mp_bcast( el_ph_ngauss, meta_ionode_id, world_comm )
  CALL mp_bcast( el_ph_nsigma, meta_ionode_id, world_comm )
  !
  ! real*8
  !
  call mp_bcast (tr2_ph, meta_ionode_id, world_comm )
  call mp_bcast (eth_rps, meta_ionode_id, world_comm )
  call mp_bcast (eth_ns, meta_ionode_id, world_comm )
  call mp_bcast (amass, meta_ionode_id, world_comm )
  call mp_bcast (alpha_mix, meta_ionode_id, world_comm )
  call mp_bcast (max_seconds, meta_ionode_id, world_comm )
  call mp_bcast (dek, meta_ionode_id, world_comm )
  CALL mp_bcast( el_ph_sigma, meta_ionode_id, world_comm )
  !
  ! characters
  !
  call mp_bcast (title, meta_ionode_id, world_comm )
  call mp_bcast (fildyn, meta_ionode_id, world_comm )
  call mp_bcast (fildvscf, meta_ionode_id, world_comm )
  call mp_bcast (fildrho, meta_ionode_id, world_comm )
  call mp_bcast (tmp_dir, meta_ionode_id, world_comm )
  call mp_bcast (prefix, meta_ionode_id, world_comm )
  call mp_bcast (electron_phonon, meta_ionode_id, world_comm )
  !
  ! derived type (one bit at a time)
  !
  call mp_bcast (drho_star%open,  meta_ionode_id, world_comm)
  call mp_bcast (drho_star%pat,   meta_ionode_id, world_comm)
  call mp_bcast (drho_star%dir,   meta_ionode_id, world_comm)
  call mp_bcast (drho_star%ext,   meta_ionode_id, world_comm)
  call mp_bcast (drho_star%basis, meta_ionode_id, world_comm)
  !
  call mp_bcast (dvscf_star%open,  meta_ionode_id, world_comm)
  call mp_bcast (dvscf_star%pat,   meta_ionode_id, world_comm)
  call mp_bcast (dvscf_star%dir,   meta_ionode_id, world_comm)
  call mp_bcast (dvscf_star%ext,   meta_ionode_id, world_comm)
  call mp_bcast (dvscf_star%basis, meta_ionode_id, world_comm)

#endif
  return
end subroutine bcast_ph_input

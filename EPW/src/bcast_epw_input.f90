  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the ionode_id, world_comm directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Code adapted from PH/bcast_ph_input - Quantum-ESPRESSO group
  ! 09/2009 Very little of this subroutine in necessary.  Many 
  ! excess variables
  !
  !-----------------------------------------------------------------------
  SUBROUTINE bcast_ph_input
  !-----------------------------------------------------------------------
  !!
  !!     In this routine the first processor sends the input to all
  !!     the other processors
  !!
#if defined(__MPI)
  USE phcom,         ONLY : zue, trans, tr2_ph, recover, nmix_ph, niter_ph, &
                            lnscf, ldisp, fildvscf, fildrho, epsil, alpha_mix 
  USE epwcom,        ONLY : epexst, epbwrite, ep_coupling, &
                            eliashberg, elecselfen, eig_read, &
                            efermi_read, dvscf_dir, delta_smear, &
                            delta_qsmear, degaussw, degaussq, conv_thr_raxis, &
                            conv_thr_racon, conv_thr_iaxis, broyden_ndim, &
                            broyden_beta, band_plot, a2f, lacon, &
                            kmaps, kerwrite, kerread, imag_read, &
                            gap_edge, fsthick, filukq, filukk, filqf, filkf, &
                            fileig, fila2f, fermi_energy, &
                            etf_mem, epwwrite, epwread, eptemp, &
                            eps_acustic, ephwrite, epbread, nsiter, nqstep, &
                            nqsmear, nqf3, nqf2, nqf1, nkf3, nkf2, nkf1, &
                            ngaussw, nest_fn,  nbndsub, nbndskip, &
                            muc, mp_mesh_q, mp_mesh_k, max_memlt, lunif, &
                            lreal, lpolar, lpade, liso, limag, laniso, &
                            specfun, lifc, asr_typ, &
                            rand_q, rand_nq, rand_nk, rand_k, pwc, phonselfen, &
                            parallel_q, parallel_k, &
                            nw_specfun, nw, nswi, nswfc, nswc, nstemp, nsmear, &
                            wsfc, wscut, write_wfn, wmin_specfun, wmin, &
                            wmax_specfun, wmax, wepexst, wannierize, &
                            vme, longrange, shortrange, system_2d, &
                            tempsmin, tempsmax, temps, delta_approx, title, &
                            scattering, scattering_serta, scattering_0rta, &
                            int_mob, scissor, carrier, ncarrier, iterative_bte, &
                            restart, restart_freq
!  USE epwcom,        ONLY : fildvscf0, tphases
  USE elph2,         ONLY : elph 
  USE mp,            ONLY : mp_bcast
  USE mp_world,      ONLY : world_comm
  USE io_files,      ONLY : prefix, tmp_dir
  USE qpoint,        ONLY : xq
  USE control_lr,    ONLY : lgamma
  USE io_global,     ONLY : meta_ionode_id
  USE control_flags, ONLY : iverbosity
  USE ions_base,     ONLY : amass
  !
  implicit none
  !
  ! logicals
  !
  CALL mp_bcast (lgamma, meta_ionode_id, world_comm)
  CALL mp_bcast (epsil, meta_ionode_id, world_comm)
  CALL mp_bcast (trans, meta_ionode_id, world_comm)
  CALL mp_bcast (zue, meta_ionode_id, world_comm)
  CALL mp_bcast (elph, meta_ionode_id, world_comm)
  CALL mp_bcast (lnscf, meta_ionode_id, world_comm)
  CALL mp_bcast (ldisp, meta_ionode_id, world_comm)
  CALL mp_bcast (elecselfen, meta_ionode_id, world_comm)!
  CALL mp_bcast (phonselfen, meta_ionode_id, world_comm)!
  CALL mp_bcast (ephwrite, meta_ionode_id, world_comm)! RM
  CALL mp_bcast (band_plot, meta_ionode_id, world_comm)! RM
  CALL mp_bcast (vme, meta_ionode_id, world_comm)!
  CALL mp_bcast (recover, meta_ionode_id, world_comm)!
  CALL mp_bcast (epbread, meta_ionode_id, world_comm)   !
  CALL mp_bcast (epbwrite, meta_ionode_id, world_comm)  !
!  CALL mp_bcast (tphases, meta_ionode_id, world_comm)   !
  CALL mp_bcast (fsthick, meta_ionode_id, world_comm)   !
  CALL mp_bcast (wmin, meta_ionode_id, world_comm)      !
  CALL mp_bcast (wmax, meta_ionode_id, world_comm)      !
  CALL mp_bcast (epwread, meta_ionode_id, world_comm)   !
  CALL mp_bcast (epwwrite, meta_ionode_id, world_comm)  !
  CALL mp_bcast (specfun, meta_ionode_id, world_comm)   !
  CALL mp_bcast (wannierize, meta_ionode_id, world_comm)! JN
  CALL mp_bcast (write_wfn, meta_ionode_id, world_comm) ! 
  CALL mp_bcast (kmaps, meta_ionode_id, world_comm) ! 
  CALL mp_bcast (nest_fn, meta_ionode_id, world_comm) ! 
  CALL mp_bcast (eig_read, meta_ionode_id, world_comm) ! 
  CALL mp_bcast (parallel_k, meta_ionode_id, world_comm) 
  CALL mp_bcast (parallel_q, meta_ionode_id, world_comm)
  CALL mp_bcast (a2f, meta_ionode_id, world_comm)
  CALL mp_bcast (etf_mem, meta_ionode_id, world_comm)
  CALL mp_bcast (rand_q, meta_ionode_id, world_comm)
  CALL mp_bcast (rand_k, meta_ionode_id, world_comm)
  CALL mp_bcast (mp_mesh_q, meta_ionode_id, world_comm)
  CALL mp_bcast (mp_mesh_k, meta_ionode_id, world_comm)
  CALL mp_bcast (wepexst, meta_ionode_id, world_comm)
  CALL mp_bcast (epexst, meta_ionode_id, world_comm)
  CALL mp_bcast (lreal, meta_ionode_id, world_comm)     ! RM
  CALL mp_bcast (limag, meta_ionode_id, world_comm)     !
  CALL mp_bcast (lpade, meta_ionode_id, world_comm)     !  
  CALL mp_bcast (lacon, meta_ionode_id, world_comm)     !
  CALL mp_bcast (liso, meta_ionode_id, world_comm)     !
  CALL mp_bcast (laniso, meta_ionode_id, world_comm)     !
  CALL mp_bcast (lpolar, meta_ionode_id, world_comm)     !
  CALL mp_bcast (lifc, meta_ionode_id, world_comm) 
  CALL mp_bcast (lunif, meta_ionode_id, world_comm)     !
  CALL mp_bcast (kerwrite, meta_ionode_id, world_comm)     !
  CALL mp_bcast (kerread, meta_ionode_id, world_comm)     !
  CALL mp_bcast (imag_read, meta_ionode_id, world_comm ) !
  CALL mp_bcast (eliashberg, meta_ionode_id, world_comm ) !
  CALL mp_bcast (ep_coupling, meta_ionode_id, world_comm ) !
  CALL mp_bcast (efermi_read, meta_ionode_id, world_comm)
  CALL mp_bcast (wmin_specfun, meta_ionode_id, world_comm)      !
  CALL mp_bcast (wmax_specfun, meta_ionode_id, world_comm)      !
  CALL mp_bcast (delta_approx, meta_ionode_id, world_comm)      !
  CALL mp_bcast (longrange, meta_ionode_id, world_comm)      !
  CALL mp_bcast (shortrange, meta_ionode_id, world_comm)      !  
  CALL mp_bcast (system_2d, meta_ionode_id, world_comm)
  CALL mp_bcast (scattering, meta_ionode_id, world_comm)
  CALL mp_bcast (scattering_serta, meta_ionode_id, world_comm)
  CALL mp_bcast (scattering_0rta, meta_ionode_id, world_comm)
  CALL mp_bcast (int_mob, meta_ionode_id, world_comm)
  CALL mp_bcast (iterative_bte, meta_ionode_id, world_comm)
  CALL mp_bcast (carrier, meta_ionode_id, world_comm)  
  CALL mp_bcast (restart, meta_ionode_id, world_comm)
  !
  ! integers
  !
  CALL mp_bcast (niter_ph, meta_ionode_id, world_comm)
  CALL mp_bcast (nmix_ph, meta_ionode_id, world_comm)
  CALL mp_bcast (iverbosity, meta_ionode_id, world_comm)
  CALL mp_bcast (ngaussw, meta_ionode_id, world_comm)     ! FG
  CALL mp_bcast (nw, meta_ionode_id, world_comm)          ! 
  CALL mp_bcast (nbndsub, meta_ionode_id, world_comm)     ! 
  CALL mp_bcast (nbndskip, meta_ionode_id, world_comm)    ! 
  CALL mp_bcast (nsmear, meta_ionode_id, world_comm)      ! 
  CALL mp_bcast (rand_nq, meta_ionode_id, world_comm)     ! 
  CALL mp_bcast (rand_nk, meta_ionode_id, world_comm)     ! 
  CALL mp_bcast (nkf1, meta_ionode_id, world_comm)
  CALL mp_bcast (nkf2, meta_ionode_id, world_comm)
  CALL mp_bcast (nkf3, meta_ionode_id, world_comm)
  CALL mp_bcast (nqf1, meta_ionode_id, world_comm)
  CALL mp_bcast (nqf2, meta_ionode_id, world_comm)
  CALL mp_bcast (nqf3, meta_ionode_id, world_comm)
  CALL mp_bcast (nqsmear, meta_ionode_id, world_comm )    ! 
  CALL mp_bcast (nqstep, meta_ionode_id, world_comm)      ! 
  CALL mp_bcast (nswfc, meta_ionode_id, world_comm )      ! 
  CALL mp_bcast (nswc, meta_ionode_id, world_comm )       !
  CALL mp_bcast (nswi, meta_ionode_id, world_comm )       !
  CALL mp_bcast (broyden_ndim, meta_ionode_id, world_comm)!
  CALL mp_bcast (nstemp, meta_ionode_id, world_comm )     !
  CALL mp_bcast (nsiter, meta_ionode_id, world_comm )     !
  CALL mp_bcast (nw_specfun, meta_ionode_id, world_comm)  !
  CALL mp_bcast (restart_freq, meta_ionode_id, world_comm)
  !
  ! real*8
  !
  CALL mp_bcast (tr2_ph, meta_ionode_id, world_comm)
  CALL mp_bcast (amass, meta_ionode_id, world_comm)
  CALL mp_bcast (alpha_mix, meta_ionode_id, world_comm)
  CALL mp_bcast (xq, meta_ionode_id, world_comm)
  CALL mp_bcast (degaussw, meta_ionode_id, world_comm)  ! FG
  CALL mp_bcast (delta_smear, meta_ionode_id, world_comm)    ! 
  CALL mp_bcast (eps_acustic, meta_ionode_id, world_comm)     ! RM
  CALL mp_bcast (degaussq, meta_ionode_id, world_comm)        !
  CALL mp_bcast (delta_qsmear, meta_ionode_id, world_comm)    ! 
  CALL mp_bcast (pwc, meta_ionode_id, world_comm )            !
  CALL mp_bcast (wsfc, meta_ionode_id, world_comm )           !
  CALL mp_bcast (wscut, meta_ionode_id, world_comm )          !
  CALL mp_bcast (broyden_beta, meta_ionode_id, world_comm )   !
  CALL mp_bcast (tempsmin, meta_ionode_id, world_comm )       !
  CALL mp_bcast (tempsmax, meta_ionode_id, world_comm )       !
  CALL mp_bcast (temps, meta_ionode_id, world_comm )       !
  CALL mp_bcast (conv_thr_raxis, meta_ionode_id, world_comm ) !
  CALL mp_bcast (conv_thr_iaxis, meta_ionode_id, world_comm ) !
  CALL mp_bcast (conv_thr_racon, meta_ionode_id, world_comm ) !
  CALL mp_bcast (gap_edge, meta_ionode_id, world_comm ) !
  CALL mp_bcast (muc, meta_ionode_id, world_comm )            !
  CALL mp_bcast (max_memlt, meta_ionode_id, world_comm)       !
  CALL mp_bcast (fermi_energy, meta_ionode_id, world_comm)    !
  CALL mp_bcast (eptemp, meta_ionode_id, world_comm)    !
  CALL mp_bcast (scissor, meta_ionode_id, world_comm)    !
  CALL mp_bcast (ncarrier, meta_ionode_id, world_comm)      
  !
  ! characters
  !
  CALL mp_bcast (title, meta_ionode_id, world_comm)
  CALL mp_bcast (fildvscf, meta_ionode_id, world_comm)
  CALL mp_bcast (fildrho, meta_ionode_id, world_comm)
  CALL mp_bcast (tmp_dir, meta_ionode_id, world_comm)
  CALL mp_bcast (prefix, meta_ionode_id, world_comm)
  !
  CALL mp_bcast (filkf, meta_ionode_id, world_comm)     ! FG
  CALL mp_bcast (filqf, meta_ionode_id, world_comm)     ! FG
  CALL mp_bcast (filukk, meta_ionode_id, world_comm)    ! FG
  CALL mp_bcast (filukq, meta_ionode_id, world_comm)    ! FG
  CALL mp_bcast (fileig, meta_ionode_id, world_comm)    ! FG
!  CALL mp_bcast (fildvscf0, meta_ionode_id, world_comm) !
  CALL mp_bcast (dvscf_dir, meta_ionode_id, world_comm)
  CALL mp_bcast (fila2f, meta_ionode_id, world_comm)     ! RM
  CALL mp_bcast (asr_typ, meta_ionode_id, world_comm)
#endif
  !
END SUBROUTINE bcast_ph_input
!
!-----------------------------------------------------------------------
SUBROUTINE bcast_ph_input1
!-----------------------------------------------------------------------
!
#if defined(__MPI)
  USE pwcom
  USE phcom
  USE mp,         ONLY: mp_bcast
  USE mp_world,   ONLY : world_comm
  USE io_global,  ONLY : meta_ionode_id
  implicit none
  !
  ! integers
  !
  CALL mp_bcast (nat_todo, meta_ionode_id, world_comm)
  IF (nat_todo.gt.0) THEN
     CALL mp_bcast (atomo, meta_ionode_id, world_comm)
  ENDIF
#endif
  !  
END SUBROUTINE bcast_ph_input1

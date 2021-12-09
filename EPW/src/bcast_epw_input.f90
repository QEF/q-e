  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  ! Code adapted from PH/bcast_ph_input - Quantum-ESPRESSO group
  ! 09/2009 Very little of this SUBROUTINE in necessary.  Many
  ! excess variables
  !
  !-----------------------------------------------------------------------
  SUBROUTINE bcast_epw_input
  !-----------------------------------------------------------------------
  !!
  !! In this routine the first processor sends the input to all the other processors
  !!
#if defined(__MPI)
  USE phcom,         ONLY : zue, trans, tr2_ph, nmix_ph, niter_ph, &
                            lnscf, ldisp, fildvscf, fildrho, epsil, alpha_mix
  USE epwcom,        ONLY : epexst, epbwrite, ep_coupling,                    &
                            eliashberg, elecselfen, eig_read, plselfen,       &
                            efermi_read, dvscf_dir, delta_smear, ngaussw,     &
                            delta_qsmear, degaussw, degaussq, conv_thr_raxis, &
                            conv_thr_racon, conv_thr_iaxis, broyden_ndim,     &
                            broyden_beta, band_plot, a2f, lacon, nest_fn,     &
                            kmaps, kerwrite, kerread, imag_read, nkc3,        &
                            gap_edge, fsthick, filqf, filkf, nqc1, nqc2, nqc3,&
                            fileig, fila2f, fermi_energy, nc, nkc1, nkc2,     &
                            etf_mem, epwwrite, epwread, nbndsub, fermi_plot,  &
                            eps_acustic, ephwrite, epbread, nsiter, nqstep,   &
                            nqsmear, nqf3, nqf2, nqf1, nkf3, nkf2, nkf1,      &
                            muc, mp_mesh_q, mp_mesh_k, max_memlt, lunif,      &
                            lreal, lpolar, lpade, liso, limag, laniso, npade, &
                            specfun_el, specfun_ph, lifc, asr_typ,            &
                            lscreen, scr_typ, fermi_diff, smear_rpa,          &
                            rand_q, rand_nq, rand_nk, rand_k, pwc, phonselfen,&
                            specfun_pl, cumulant, bnd_cum, iterative_bte,     &
                            nw_specfun, nw, nswi, nswfc, nswc, nstemp, nsmear,&
                            wsfc, wscut, write_wfn, wmin_specfun, wmin,       &
                            wmax_specfun, wmax, wepexst, wannierize,          &
                            vme, longrange, shortrange, system_2d, lindabs,   &
                            temps, tempsmin, tempsmax, delta_approx, title,   &
                            scattering, scattering_serta, scattering_0rta,    &
                            int_mob, scissor, carrier, ncarrier,              &
                            restart, restart_step, prtgkk, nel, meff, epsiheg,&
                            scatread, restart, restart_step, restart_filq,    &
                            lphase, omegamin, omegamax, omegastep, n_r,       &
                            mob_maxiter, use_ws, epmatkqread, selecqread,     &
                            scdm_proj, scdm_entanglement, scdm_mu, scdm_sigma,&
                            assume_metal, wannier_plot_scale, reduce_unk,     &
                            wannier_plot_supercell, wannier_plot_radius,      &
                            fixsym, epw_no_t_rev, epw_tr, epw_nosym, epw_noinv, &
                            epw_crysym, bfieldx, bfieldy, bfieldz, tc_linear, &
                            tc_linear_solver, mob_maxfreq, mob_nfreq
  USE elph2,         ONLY : elph
  USE mp,            ONLY : mp_bcast
  USE mp_world,      ONLY : world_comm
  USE io_files,      ONLY : prefix, tmp_dir
  USE qpoint,        ONLY : xq
  USE io_global,     ONLY : meta_ionode_id
  USE control_flags, ONLY : iverbosity
  USE ions_base,     ONLY : amass
  ! ---------------------------------------------------------------------------------
  ! Added for polaron calculations. Originally by Danny Sio, modified by Chao Lian.
  ! Shell implementation for future use.
  USE epwcom,        ONLY : wfcelec, model_vertex , polaron_wf, r01, r02, r03,&
                            num_cbands, start_band, start_mode, cb_shift,     &
                            polaron_interpol, polaron_bq, polaron_dos,        &
                            electron_dos, phonon_dos, diag_mode,              &
                            restart_polaron_mode, polaron_type,               &
                            emax_plrn, nDOS_plrn, emin_plrn
  ! -------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  ! logicals
  !
  CALL mp_bcast(epsil           , meta_ionode_id, world_comm)
  CALL mp_bcast(trans           , meta_ionode_id, world_comm)
  CALL mp_bcast(zue             , meta_ionode_id, world_comm)
  CALL mp_bcast(elph            , meta_ionode_id, world_comm)
  CALL mp_bcast(lnscf           , meta_ionode_id, world_comm)
  CALL mp_bcast(ldisp           , meta_ionode_id, world_comm)
  CALL mp_bcast(elecselfen      , meta_ionode_id, world_comm)
  CALL mp_bcast(phonselfen      , meta_ionode_id, world_comm)
  CALL mp_bcast(plselfen        , meta_ionode_id, world_comm)
  CALL mp_bcast(ephwrite        , meta_ionode_id, world_comm)
  CALL mp_bcast(band_plot       , meta_ionode_id, world_comm)
  CALL mp_bcast(fermi_plot      , meta_ionode_id, world_comm)
  CALL mp_bcast(vme             , meta_ionode_id, world_comm)
  CALL mp_bcast(epbread         , meta_ionode_id, world_comm)
  CALL mp_bcast(epbwrite        , meta_ionode_id, world_comm)
  CALL mp_bcast(fsthick         , meta_ionode_id, world_comm)
  CALL mp_bcast(wmin            , meta_ionode_id, world_comm)
  CALL mp_bcast(wmax            , meta_ionode_id, world_comm)
  CALL mp_bcast(epwread         , meta_ionode_id, world_comm)
  CALL mp_bcast(epwwrite        , meta_ionode_id, world_comm)
  CALL mp_bcast(specfun_el      , meta_ionode_id, world_comm)
  CALL mp_bcast(specfun_ph      , meta_ionode_id, world_comm)
  CALL mp_bcast(specfun_pl      , meta_ionode_id, world_comm)
  CALL mp_bcast(wannierize      , meta_ionode_id, world_comm)
  CALL mp_bcast(write_wfn       , meta_ionode_id, world_comm)
  CALL mp_bcast(kmaps           , meta_ionode_id, world_comm)
  CALL mp_bcast(nest_fn         , meta_ionode_id, world_comm)
  CALL mp_bcast(eig_read        , meta_ionode_id, world_comm)
  CALL mp_bcast(a2f             , meta_ionode_id, world_comm)
  CALL mp_bcast(etf_mem         , meta_ionode_id, world_comm)
  CALL mp_bcast(rand_q          , meta_ionode_id, world_comm)
  CALL mp_bcast(rand_k          , meta_ionode_id, world_comm)
  CALL mp_bcast(mp_mesh_q       , meta_ionode_id, world_comm)
  CALL mp_bcast(mp_mesh_k       , meta_ionode_id, world_comm)
  CALL mp_bcast(wepexst         , meta_ionode_id, world_comm)
  CALL mp_bcast(epexst          , meta_ionode_id, world_comm)
  CALL mp_bcast(lreal           , meta_ionode_id, world_comm)
  CALL mp_bcast(limag           , meta_ionode_id, world_comm)
  CALL mp_bcast(lpade           , meta_ionode_id, world_comm)
  CALL mp_bcast(lacon           , meta_ionode_id, world_comm)
  CALL mp_bcast(liso            , meta_ionode_id, world_comm)
  CALL mp_bcast(laniso          , meta_ionode_id, world_comm)
  CALL mp_bcast(tc_linear       , meta_ionode_id, world_comm)
  CALL mp_bcast(tc_linear_solver, meta_ionode_id, world_comm)
  CALL mp_bcast(lpolar          , meta_ionode_id, world_comm)
  CALL mp_bcast(lifc            , meta_ionode_id, world_comm)
  CALL mp_bcast(lscreen         , meta_ionode_id, world_comm)
  CALL mp_bcast(cumulant        , meta_ionode_id, world_comm)
  CALL mp_bcast(lunif           , meta_ionode_id, world_comm)
  CALL mp_bcast(kerwrite        , meta_ionode_id, world_comm)
  CALL mp_bcast(kerread         , meta_ionode_id, world_comm)
  CALL mp_bcast(imag_read       , meta_ionode_id, world_comm)
  CALL mp_bcast(eliashberg      , meta_ionode_id, world_comm)
  CALL mp_bcast(ep_coupling     , meta_ionode_id, world_comm)
  CALL mp_bcast(efermi_read     , meta_ionode_id, world_comm)
  CALL mp_bcast(wmin_specfun    , meta_ionode_id, world_comm)
  CALL mp_bcast(wmax_specfun    , meta_ionode_id, world_comm)
  CALL mp_bcast(delta_approx    , meta_ionode_id, world_comm)
  CALL mp_bcast(longrange       , meta_ionode_id, world_comm)
  CALL mp_bcast(shortrange      , meta_ionode_id, world_comm)
  CALL mp_bcast(system_2d       , meta_ionode_id, world_comm)
  CALL mp_bcast(scattering      , meta_ionode_id, world_comm)
  CALL mp_bcast(scattering_serta, meta_ionode_id, world_comm)
  CALL mp_bcast(scatread        , meta_ionode_id, world_comm)
  CALL mp_bcast(scattering_0rta , meta_ionode_id, world_comm)
  CALL mp_bcast(int_mob         , meta_ionode_id, world_comm)
  CALL mp_bcast(iterative_bte   , meta_ionode_id, world_comm)
  CALL mp_bcast(carrier         , meta_ionode_id, world_comm)
  CALL mp_bcast(restart         , meta_ionode_id, world_comm)
  CALL mp_bcast(prtgkk          , meta_ionode_id, world_comm)
  CALL mp_bcast(lphase          , meta_ionode_id, world_comm)
  CALL mp_bcast(lindabs         , meta_ionode_id, world_comm)
  CALL mp_bcast(use_ws          , meta_ionode_id, world_comm)
  CALL mp_bcast(epmatkqread     , meta_ionode_id, world_comm)
  CALL mp_bcast(selecqread      , meta_ionode_id, world_comm)
  CALL mp_bcast(scdm_proj       , meta_ionode_id, world_comm)
  CALL mp_bcast(assume_metal    , meta_ionode_id, world_comm)
  CALL mp_bcast(reduce_unk      , meta_ionode_id, world_comm)
  CALL mp_bcast(fixsym          , meta_ionode_id, world_comm)
  CALL mp_bcast(epw_no_t_rev    , meta_ionode_id, world_comm)
  CALL mp_bcast(epw_tr          , meta_ionode_id, world_comm)
  CALL mp_bcast(epw_nosym       , meta_ionode_id, world_comm)
  CALL mp_bcast(epw_noinv       , meta_ionode_id, world_comm)
  CALL mp_bcast(epw_crysym      , meta_ionode_id, world_comm)
  !
  ! integers
  !
  CALL mp_bcast(niter_ph    , meta_ionode_id, world_comm)
  CALL mp_bcast(nmix_ph     , meta_ionode_id, world_comm)
  CALL mp_bcast(iverbosity  , meta_ionode_id, world_comm)
  CALL mp_bcast(ngaussw     , meta_ionode_id, world_comm)
  CALL mp_bcast(nw          , meta_ionode_id, world_comm)
  CALL mp_bcast(nbndsub     , meta_ionode_id, world_comm)
  CALL mp_bcast(nsmear      , meta_ionode_id, world_comm)
  CALL mp_bcast(rand_nq     , meta_ionode_id, world_comm)
  CALL mp_bcast(rand_nk     , meta_ionode_id, world_comm)
  CALL mp_bcast(nkc1        , meta_ionode_id, world_comm)
  CALL mp_bcast(nkc2        , meta_ionode_id, world_comm)
  CALL mp_bcast(nkc3        , meta_ionode_id, world_comm)
  CALL mp_bcast(nqc1        , meta_ionode_id, world_comm)
  CALL mp_bcast(nqc2        , meta_ionode_id, world_comm)
  CALL mp_bcast(nqc3        , meta_ionode_id, world_comm)
  CALL mp_bcast(nkf1        , meta_ionode_id, world_comm)
  CALL mp_bcast(nkf2        , meta_ionode_id, world_comm)
  CALL mp_bcast(nkf3        , meta_ionode_id, world_comm)
  CALL mp_bcast(nqf1        , meta_ionode_id, world_comm)
  CALL mp_bcast(nqf2        , meta_ionode_id, world_comm)
  CALL mp_bcast(nqf3        , meta_ionode_id, world_comm)
  CALL mp_bcast(nqsmear     , meta_ionode_id, world_comm)
  CALL mp_bcast(nqstep      , meta_ionode_id, world_comm)
  CALL mp_bcast(nswfc       , meta_ionode_id, world_comm)
  CALL mp_bcast(nswc        , meta_ionode_id, world_comm)
  CALL mp_bcast(nswi        , meta_ionode_id, world_comm)
  CALL mp_bcast(broyden_ndim, meta_ionode_id, world_comm)
  CALL mp_bcast(nstemp      , meta_ionode_id, world_comm)
  CALL mp_bcast(nsiter      , meta_ionode_id, world_comm)
  CALL mp_bcast(npade       , meta_ionode_id, world_comm)
  CALL mp_bcast(nw_specfun  , meta_ionode_id, world_comm)
  CALL mp_bcast(restart_step, meta_ionode_id, world_comm)
  CALL mp_bcast(scr_typ     , meta_ionode_id, world_comm)
  CALL mp_bcast(bnd_cum     , meta_ionode_id, world_comm)
  CALL mp_bcast(mob_maxiter , meta_ionode_id, world_comm)
  CALL mp_bcast(wannier_plot_supercell, meta_ionode_id, world_comm)
  CALL mp_bcast(mob_nfreq  , meta_ionode_id, world_comm)
  !
  ! REAL*8
  !
  CALL mp_bcast(tr2_ph        , meta_ionode_id, world_comm)
  CALL mp_bcast(amass         , meta_ionode_id, world_comm)
  CALL mp_bcast(alpha_mix     , meta_ionode_id, world_comm)
  CALL mp_bcast(xq            , meta_ionode_id, world_comm)
  CALL mp_bcast(degaussw      , meta_ionode_id, world_comm)
  CALL mp_bcast(delta_smear   , meta_ionode_id, world_comm)
  CALL mp_bcast(eps_acustic   , meta_ionode_id, world_comm)
  CALL mp_bcast(degaussq      , meta_ionode_id, world_comm)
  CALL mp_bcast(delta_qsmear  , meta_ionode_id, world_comm)
  CALL mp_bcast(pwc           , meta_ionode_id, world_comm)
  CALL mp_bcast(wsfc          , meta_ionode_id, world_comm)
  CALL mp_bcast(wscut         , meta_ionode_id, world_comm)
  CALL mp_bcast(broyden_beta  , meta_ionode_id, world_comm)
  CALL mp_bcast(tempsmin      , meta_ionode_id, world_comm)
  CALL mp_bcast(tempsmax      , meta_ionode_id, world_comm)
  CALL mp_bcast(temps         , meta_ionode_id, world_comm)
  CALL mp_bcast(conv_thr_raxis, meta_ionode_id, world_comm)
  CALL mp_bcast(conv_thr_iaxis, meta_ionode_id, world_comm)
  CALL mp_bcast(conv_thr_racon, meta_ionode_id, world_comm)
  CALL mp_bcast(gap_edge      , meta_ionode_id, world_comm)
  CALL mp_bcast(muc           , meta_ionode_id, world_comm)
  CALL mp_bcast(max_memlt     , meta_ionode_id, world_comm)
  CALL mp_bcast(fermi_energy  , meta_ionode_id, world_comm)
  CALL mp_bcast(scissor       , meta_ionode_id, world_comm)
  CALL mp_bcast(ncarrier      , meta_ionode_id, world_comm)
  CALL mp_bcast(nel           , meta_ionode_id, world_comm)
  CALL mp_bcast(meff          , meta_ionode_id, world_comm)
  CALL mp_bcast(epsiheg       , meta_ionode_id, world_comm)
  CALL mp_bcast(fermi_diff    , meta_ionode_id, world_comm)
  CALL mp_bcast(smear_rpa     , meta_ionode_id, world_comm)
  CALL mp_bcast(omegamin      , meta_ionode_id, world_comm)
  CALL mp_bcast(omegamax      , meta_ionode_id, world_comm)
  CALL mp_bcast(omegastep     , meta_ionode_id, world_comm)
  CALL mp_bcast(n_r           , meta_ionode_id, world_comm)
  CALL mp_bcast(nc            , meta_ionode_id, world_comm)
  CALL mp_bcast(scdm_mu       , meta_ionode_id, world_comm)
  CALL mp_bcast(scdm_sigma    , meta_ionode_id, world_comm)
  CALL mp_bcast(wannier_plot_radius, meta_ionode_id, world_comm)
  CALL mp_bcast(wannier_plot_scale, meta_ionode_id, world_comm)
  CALL mp_bcast(bfieldx       , meta_ionode_id, world_comm)
  CALL mp_bcast(bfieldy       , meta_ionode_id, world_comm)
  CALL mp_bcast(bfieldz       , meta_ionode_id, world_comm)
  CALL mp_bcast(mob_maxfreq   , meta_ionode_id, world_comm)
  !
  ! characters
  !
  CALL mp_bcast(title            , meta_ionode_id, world_comm)
  CALL mp_bcast(fildvscf         , meta_ionode_id, world_comm)
  CALL mp_bcast(fildrho          , meta_ionode_id, world_comm)
  CALL mp_bcast(tmp_dir          , meta_ionode_id, world_comm)
  CALL mp_bcast(prefix           , meta_ionode_id, world_comm)
  CALL mp_bcast(filkf            , meta_ionode_id, world_comm)
  CALL mp_bcast(filqf            , meta_ionode_id, world_comm)
  CALL mp_bcast(fileig           , meta_ionode_id, world_comm)
  CALL mp_bcast(dvscf_dir        , meta_ionode_id, world_comm)
  CALL mp_bcast(fila2f           , meta_ionode_id, world_comm)
  CALL mp_bcast(restart_filq     , meta_ionode_id, world_comm)
  CALL mp_bcast(asr_typ          , meta_ionode_id, world_comm)
  CALL mp_bcast(scdm_entanglement, meta_ionode_id, world_comm)
  !
  ! ---------------------------------------------------------------------------------
  ! Added for polaron calculations. Originally by Danny Sio, modified by Chao Lian.
  ! Shell implementation for future use.
  CALL mp_bcast (wfcelec         , meta_ionode_id, world_comm)
  CALL mp_bcast (model_vertex    , meta_ionode_id, world_comm)
  CALL mp_bcast (polaron_wf      , meta_ionode_id, world_comm)
  CALL mp_bcast (polaron_interpol, meta_ionode_id, world_comm)
  CALL mp_bcast (polaron_bq      , meta_ionode_id, world_comm)
  CALL mp_bcast (polaron_dos     , meta_ionode_id, world_comm)
  CALL mp_bcast (electron_dos    , meta_ionode_id, world_comm)
  CALL mp_bcast (phonon_dos      , meta_ionode_id, world_comm)
  CALL mp_bcast (num_cbands  , meta_ionode_id, world_comm)
  CALL mp_bcast (start_band  , meta_ionode_id, world_comm)
  CALL mp_bcast (start_mode  , meta_ionode_id, world_comm)
  CALL mp_bcast (cb_shift    , meta_ionode_id, world_comm)
  CALL mp_bcast (diag_mode   , meta_ionode_id, world_comm)
  CALL mp_bcast (restart_polaron_mode, meta_ionode_id, world_comm)
  CALL mp_bcast (polaron_type, meta_ionode_id, world_comm)
  CALL mp_bcast (r01           , meta_ionode_id, world_comm)
  CALL mp_bcast (r02           , meta_ionode_id, world_comm)
  CALL mp_bcast (r03           , meta_ionode_id, world_comm)
  CALL mp_bcast (nDOS_plrn     , meta_ionode_id, world_comm)
  CALL mp_bcast (emax_plrn     , meta_ionode_id, world_comm)
  CALL mp_bcast (emin_plrn     , meta_ionode_id, world_comm)
  ! --------------------------------------------------------------------------------
#endif
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE bcast_epw_input
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE bcast_epw_input1
  !-----------------------------------------------------------------------
  !!
  !! Second batch of broadcasting.
  !!
#if defined(__MPI)
  USE partial,    ONLY : nat_todo, atomo
  USE mp,         ONLY : mp_bcast
  USE mp_world,   ONLY : world_comm
  USE io_global,  ONLY : meta_ionode_id
  !
  IMPLICIT NONE
  !
  ! integers
  !
  CALL mp_bcast (nat_todo, meta_ionode_id, world_comm)
  IF (nat_todo > 0) THEN
     CALL mp_bcast (atomo, meta_ionode_id, world_comm)
  ENDIF
#endif
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE bcast_epw_input1
  !-----------------------------------------------------------------------

  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
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
  SUBROUTINE bcast_input
  !-----------------------------------------------------------------------
  !!
  !! In this routine the first processor sends the input to all the other processors
  !!
#if defined(__MPI)
  USE output,        ONLY : fildvscf, fildrho
  USE input,         ONLY : epexst, epbwrite, ep_coupling,                    &
                            eliashberg, elecselfen, eig_read, plselfen,       &
                            efermi_read, dvscf_dir, delta_smear, ngaussw,     &
                            delta_qsmear, degaussw, degaussq, conv_thr_raxis, &
                            conv_thr_racon, conv_thr_iaxis, broyden_ndim,     &
                            broyden_beta, band_plot, a2f, lacon, nest_fn,     &
                            kmaps, kerwrite, kerread, imag_read, nkc3,        &
                            gap_edge, fsthick, filqf, filkf, nqc1, nqc2, nqc3,&
                            fileig, fila2f, fermi_energy, nc, nkc1, nkc2,     &
                            etf_mem, epwwrite, epwread, nbndsub, fermi_plot,  &
                            eps_acoustic, ephwrite, epbread, nsiter, nqstep,   &
                            nqsmear, nqf3, nqf2, nqf1, nkf3, nkf2, nkf1,      &
                            muc, mp_mesh_q, mp_mesh_k, max_memlt,             &
                            lreal, lpolar, lpade, liso, limag, laniso, npade, &
                            specfun_el, specfun_ph, lifc, asr_typ,            &
                            lscreen, scr_typ, fermi_diff, smear_rpa,          &
                            rand_q, rand_nq, rand_nk, rand_k, phonselfen,     &
                            specfun_pl, lcumulant, bnd_cum, iterative_bte,    &
                            nw_specfun, nw, nswi, nstemp, nsmear,             &
                            wscut, write_wfn, wmin_specfun, wmin, epw_noinv,  &
                            wmax_specfun, wmax, wepexst, wannierize,          &
                            vme, longrange_only, shortrange, system_2d,       &
                            temps, tempsmin, tempsmax, delta_approx, title,   &
                            scattering, scattering_serta, scattering_0rta,    &
                            int_mob, scissor, carrier, ncarrier, lindabs,     &
                            restart, restart_step, prtgkk, nel, meff, epsiheg,&
                            scatread, restart, restart_step, restart_filq,    &
                            lphase, omegamin, omegamax, omegastep, sigma_ref, &
                            mob_maxiter, use_ws, epmatkqread, selecqread,     &
                            scdm_proj, scdm_entanglement, scdm_mu, scdm_sigma,&
                            assume_metal, wannier_plot_scale, reduce_unk,     &
                            wannier_plot_supercell, wannier_plot_radius,      &
                            fixsym, epw_no_t_rev, epw_tr, epw_nosym,          &
                            epw_crysym, bfieldx, bfieldy, bfieldz, tc_linear, &
                            ii_g, ii_charge, ii_n, ii_scattering, ii_only,    &
                            ii_lscreen, ii_eda, ii_partion,                   &
                            gb_scattering, gb_only, gb_size,                  &
                            tc_linear_solver, mob_maxfreq, mob_nfreq,         &
                            fbw, gridsamp, griddens, dos_del, muchem, a_gap0, &
                            filirobj, icoulomb, eps_cut_ir, positive_matsu,   &
                            emax_coulomb, emin_coulomb, filnscf_coul,         &
                            lwfpt, compute_dmat, ahc_nbnd, ahc_nbndskip,      &
                            ahc_win_min, ahc_win_max, plrn,                   &
                            restart_plrn, conv_thr_plrn, end_band_plrn, lrot, &
                            cal_psir_plrn, start_band_plrn,  type_plrn,       &
                            nstate_plrn, interp_Ank_plrn, interp_Bqu_plrn,    &
                            init_sigma_plrn, full_diagon_plrn, mixing_Plrn,   &
                            init_plrn, niter_plrn, nDOS_plrn, edos_max_plrn,  &
                            edos_min_plrn, edos_sigma_plrn, pdos_sigma_plrn,  &
                            pdos_max_plrn, pdos_min_plrn, seed_plrn,          &
                            ethrdg_plrn, time_rev_A_plrn, nhblock_plrn,       &
                            beta_plrn, Mmn_plrn, recal_Mmn_plrn, r0_plrn,     &
                            debug_plrn, time_rev_U_plrn, g_start_band_plrn,   &
                            g_end_band_plrn, g_start_energy_plrn,             &
                            g_end_energy_plrn, model_vertex_plrn,             &
                            model_enband_plrn, model_phfreq_plrn, kappa_plrn, &
                            omega_LO_plrn, m_eff_plrn, step_wf_grid_plrn,     &
                            adapt_ethrdg_plrn, io_lvl_plrn, nethrdg_plrn,     &
                            david_ndim_plrn, init_ntau_plrn, g_tol_plrn,      &
                            init_ethrdg_plrn, init_k0_plrn,                   &
                            scell_mat_plrn, scell_mat, len_mesh, meshnum,     &
                            loptabs, wf_quasi, nq_init, start_mesh,           &
                            DW, mode_res, QD_min, QD_bin, do_CHBB,            &
                            elecselfen_type, calc_nelec_wann, lopt_w2b,       &
                            lfast_kmesh, epw_memdist
  USE global_var,    ONLY : elph
  USE mp,            ONLY : mp_bcast
  USE mp_world,      ONLY : world_comm
  USE io_files,      ONLY : prefix, tmp_dir
  USE qpoint,        ONLY : xq
  USE io_global,     ONLY : meta_ionode_id
  USE control_flags, ONLY : iverbosity
  USE ions_base,     ONLY : amass
  !
  IMPLICIT NONE
  !
  ! logicals
  !
  CALL mp_bcast(elph            , meta_ionode_id, world_comm)
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
  CALL mp_bcast(fbw             , meta_ionode_id, world_comm)
  CALL mp_bcast(gridsamp        , meta_ionode_id, world_comm)
  CALL mp_bcast(griddens        , meta_ionode_id, world_comm)
  CALL mp_bcast(dos_del         , meta_ionode_id, world_comm)
  CALL mp_bcast(muchem          , meta_ionode_id, world_comm)
  CALL mp_bcast(positive_matsu  , meta_ionode_id, world_comm)
  CALL mp_bcast(lpolar          , meta_ionode_id, world_comm)
  CALL mp_bcast(lifc            , meta_ionode_id, world_comm)
  CALL mp_bcast(lscreen         , meta_ionode_id, world_comm)
  CALL mp_bcast(lcumulant       , meta_ionode_id, world_comm)
  CALL mp_bcast(kerwrite        , meta_ionode_id, world_comm)
  CALL mp_bcast(kerread         , meta_ionode_id, world_comm)
  CALL mp_bcast(imag_read       , meta_ionode_id, world_comm)
  CALL mp_bcast(eliashberg      , meta_ionode_id, world_comm)
  CALL mp_bcast(ep_coupling     , meta_ionode_id, world_comm)
  CALL mp_bcast(efermi_read     , meta_ionode_id, world_comm)
  CALL mp_bcast(wmin_specfun    , meta_ionode_id, world_comm)
  CALL mp_bcast(wmax_specfun    , meta_ionode_id, world_comm)
  CALL mp_bcast(delta_approx    , meta_ionode_id, world_comm)
  CALL mp_bcast(longrange_only  , meta_ionode_id, world_comm)
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
  CALL mp_bcast(ii_g            , meta_ionode_id, world_comm)
  CALL mp_bcast(ii_scattering   , meta_ionode_id, world_comm)
  CALL mp_bcast(ii_only         , meta_ionode_id, world_comm)
  CALL mp_bcast(ii_lscreen      , meta_ionode_id, world_comm)
  CALL mp_bcast(ii_partion      , meta_ionode_id, world_comm)
  CALL mp_bcast(gb_scattering   , meta_ionode_id, world_comm)
  CALL mp_bcast(gb_only         , meta_ionode_id, world_comm)
  CALL mp_bcast(adapt_ethrdg_plrn, meta_ionode_id, world_comm)
  CALL mp_bcast(calc_nelec_wann , meta_ionode_id, world_comm)
  CALL mp_bcast(lopt_w2b        , meta_ionode_id, world_comm)
  CALL mp_bcast(lfast_kmesh     , meta_ionode_id, world_comm)
  CALL mp_bcast(epw_memdist     , meta_ionode_id, world_comm)
  !
  ! integers
  !
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
  CALL mp_bcast(mob_nfreq   , meta_ionode_id, world_comm)
  CALL mp_bcast(io_lvl_plrn , meta_ionode_id, world_comm)
  CALL mp_bcast(nethrdg_plrn, meta_ionode_id, world_comm)
  CALL mp_bcast(david_ndim_plrn, meta_ionode_id, world_comm)
  CALL mp_bcast(init_ntau_plrn, meta_ionode_id, world_comm)
  CALL mp_bcast(icoulomb    , meta_ionode_id, world_comm)
  !
  ! REAL*8
  !
  CALL mp_bcast(amass         , meta_ionode_id, world_comm)
  CALL mp_bcast(xq            , meta_ionode_id, world_comm)
  CALL mp_bcast(degaussw      , meta_ionode_id, world_comm)
  CALL mp_bcast(delta_smear   , meta_ionode_id, world_comm)
  CALL mp_bcast(eps_acoustic  , meta_ionode_id, world_comm)
  CALL mp_bcast(degaussq      , meta_ionode_id, world_comm)
  CALL mp_bcast(delta_qsmear  , meta_ionode_id, world_comm)
  CALL mp_bcast(wscut         , meta_ionode_id, world_comm)
  CALL mp_bcast(broyden_beta  , meta_ionode_id, world_comm)
  CALL mp_bcast(tempsmin      , meta_ionode_id, world_comm)
  CALL mp_bcast(tempsmax      , meta_ionode_id, world_comm)
  CALL mp_bcast(temps         , meta_ionode_id, world_comm)
  CALL mp_bcast(conv_thr_raxis, meta_ionode_id, world_comm)
  CALL mp_bcast(conv_thr_iaxis, meta_ionode_id, world_comm)
  CALL mp_bcast(conv_thr_racon, meta_ionode_id, world_comm)
  CALL mp_bcast(gap_edge      , meta_ionode_id, world_comm)
  CALL mp_bcast(a_gap0        , meta_ionode_id, world_comm)
  CALL mp_bcast(muc           , meta_ionode_id, world_comm)
  CALL mp_bcast(emax_coulomb  , meta_ionode_id, world_comm)
  CALL mp_bcast(emin_coulomb  , meta_ionode_id, world_comm)
  CALL mp_bcast(eps_cut_ir    , meta_ionode_id, world_comm)
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
  CALL mp_bcast(sigma_ref     , meta_ionode_id, world_comm)
  CALL mp_bcast(nc            , meta_ionode_id, world_comm)
  CALL mp_bcast(scdm_mu       , meta_ionode_id, world_comm)
  CALL mp_bcast(scdm_sigma    , meta_ionode_id, world_comm)
  CALL mp_bcast(wannier_plot_radius, meta_ionode_id, world_comm)
  CALL mp_bcast(wannier_plot_scale, meta_ionode_id, world_comm)
  CALL mp_bcast(bfieldx       , meta_ionode_id, world_comm)
  CALL mp_bcast(bfieldy       , meta_ionode_id, world_comm)
  CALL mp_bcast(bfieldz       , meta_ionode_id, world_comm)
  CALL mp_bcast(mob_maxfreq   , meta_ionode_id, world_comm)
  CALL mp_bcast(ii_charge     , meta_ionode_id, world_comm)
  CALL mp_bcast(ii_n          , meta_ionode_id, world_comm)
  CALL mp_bcast(ii_eda        , meta_ionode_id, world_comm)
  CALL mp_bcast(gb_size       , meta_ionode_id, world_comm)
  CALL mp_bcast(g_tol_plrn    , meta_ionode_id, world_comm)
  CALL mp_bcast(init_ethrdg_plrn, meta_ionode_id, world_comm)
  CALL mp_bcast(init_k0_plrn  , meta_ionode_id, world_comm)
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
  CALL mp_bcast(lrot             , meta_ionode_id, world_comm)
  CALL mp_bcast(plrn             , meta_ionode_id, world_comm)
  CALL mp_bcast(cal_psir_plrn    , meta_ionode_id, world_comm)
  CALL mp_bcast(interp_Ank_plrn  , meta_ionode_id, world_comm)
  CALL mp_bcast(interp_Bqu_plrn  , meta_ionode_id, world_comm)
  CALL mp_bcast(start_band_plrn  , meta_ionode_id, world_comm)
  CALL mp_bcast(end_band_plrn    , meta_ionode_id, world_comm)
  CALL mp_bcast(type_plrn        , meta_ionode_id, world_comm)
  CALL mp_bcast(nDOS_plrn        , meta_ionode_id, world_comm)
  CALL mp_bcast(edos_max_plrn    , meta_ionode_id, world_comm)
  CALL mp_bcast(edos_min_plrn    , meta_ionode_id, world_comm)
  CALL mp_bcast(nstate_plrn      , meta_ionode_id, world_comm)
  CALL mp_bcast(niter_plrn       , meta_ionode_id, world_comm)
  CALL mp_bcast(conv_thr_plrn    , meta_ionode_id, world_comm)
  CAll mp_bcast(restart_plrn     , meta_ionode_id, world_comm)
  CAll mp_bcast(type_plrn        , meta_ionode_id, world_comm)
  CAll mp_bcast(init_sigma_plrn  , meta_ionode_id, world_comm)
  CAll mp_bcast(ethrdg_plrn      , meta_ionode_id, world_comm)
  CAll mp_bcast(time_rev_A_plrn  , meta_ionode_id, world_comm)
  CAll mp_bcast(time_rev_U_plrn  , meta_ionode_id, world_comm)
  CAll mp_bcast(debug_plrn       , meta_ionode_id, world_comm)
  CAll mp_bcast(full_diagon_plrn , meta_ionode_id, world_comm)
  CAll mp_bcast(mixing_Plrn      , meta_ionode_id, world_comm)
  CAll mp_bcast(init_plrn        , meta_ionode_id, world_comm)
  CAll mp_bcast(Mmn_plrn         , meta_ionode_id, world_comm)
  CAll mp_bcast(recal_Mmn_plrn   , meta_ionode_id, world_comm)
  CAll mp_bcast(r0_plrn          , meta_ionode_id, world_comm)
  CAll mp_bcast(edos_sigma_plrn  , meta_ionode_id, world_comm)
  CAll mp_bcast(pdos_sigma_plrn  , meta_ionode_id, world_comm)
  CAll mp_bcast(pdos_max_plrn    , meta_ionode_id, world_comm)
  CAll mp_bcast(pdos_min_plrn    , meta_ionode_id, world_comm)
  CAll mp_bcast(seed_plrn        , meta_ionode_id, world_comm)
  CAll mp_bcast(nhblock_plrn     , meta_ionode_id, world_comm)
  CAll mp_bcast(beta_plrn        , meta_ionode_id, world_comm)
  CAll mp_bcast(g_start_band_plrn, meta_ionode_id, world_comm)
  CAll mp_bcast(g_end_band_plrn  , meta_ionode_id, world_comm)
  CAll mp_bcast(g_start_energy_plrn, meta_ionode_id, world_comm)
  CAll mp_bcast(g_end_energy_plrn, meta_ionode_id, world_comm)
  CAll mp_bcast(step_wf_grid_plrn, meta_ionode_id, world_comm)
  CAll mp_bcast(model_vertex_plrn, meta_ionode_id, world_comm)
  CAll mp_bcast(model_enband_plrn, meta_ionode_id, world_comm)
  CAll mp_bcast(model_phfreq_plrn, meta_ionode_id, world_comm)
  CAll mp_bcast(kappa_plrn       , meta_ionode_id, world_comm)
  CAll mp_bcast(omega_LO_plrn    , meta_ionode_id, world_comm)
  CAll mp_bcast(m_eff_plrn       , meta_ionode_id, world_comm)
  CALL mp_bcast(scell_mat_plrn   , meta_ionode_id, world_comm)
  CALL mp_bcast(scell_mat        , meta_ionode_id, world_comm)
  CALL mp_bcast(loptabs          , meta_ionode_id, world_comm)
  CALL mp_bcast(meshnum          , meta_ionode_id, world_comm)
  CALL mp_bcast(len_mesh         , meta_ionode_id, world_comm)
  CALL mp_bcast(wf_quasi         , meta_ionode_id, world_comm)
  CALL mp_bcast(nq_init          , meta_ionode_id, world_comm)
  CALL mp_bcast(start_mesh       , meta_ionode_id, world_comm)
  CALL mp_bcast(DW               , meta_ionode_id, world_comm)
  CALL mp_bcast(mode_res         , meta_ionode_id, world_comm)
  CALL mp_bcast(QD_min           , meta_ionode_id, world_comm)
  CALL mp_bcast(QD_bin           , meta_ionode_id, world_comm)
  CALL mp_bcast(do_CHBB          , meta_ionode_id, world_comm)
  CALL mp_bcast(lwfpt            , meta_ionode_id, world_comm)
  CALL mp_bcast(compute_dmat     , meta_ionode_id, world_comm)
  CALL mp_bcast(ahc_nbnd         , meta_ionode_id, world_comm)
  CALL mp_bcast(ahc_nbndskip     , meta_ionode_id, world_comm)
  CALL mp_bcast(ahc_win_min      , meta_ionode_id, world_comm)
  CALL mp_bcast(ahc_win_max      , meta_ionode_id, world_comm)
  CALL mp_bcast(elecselfen_type  , meta_ionode_id, world_comm)
  CALL mp_bcast(filirobj         , meta_ionode_id, world_comm)
  CALL mp_bcast(filnscf_coul     , meta_ionode_id, world_comm)
  !
#endif
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE bcast_input
  !-----------------------------------------------------------------------

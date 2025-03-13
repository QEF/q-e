  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.ORg/copyleft.gpl.txt .
  !
  ! Adapted from the code PH/phq_readin - Quantum-ESPRESSO group
  !-----------------------------------------------------------------------
  SUBROUTINE readin()
  !-----------------------------------------------------------------------
  !!
  !! This routine reads the control variables for the program epw.
  !! A second routine, read_file, reads the variables saved on a file
  !! by the self-consistent program.
  !!
  !! @Note:
  !!   SP: Image parallelization added
  !!
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ntyp => nsp
  USE cell_base,     ONLY : at
  USE mp,            ONLY : mp_bcast
  USE wvfct,         ONLY : nbnd, et
  USE klist,         ONLY : nks, xk, nkstot
  USE lsda_mod,      ONLY : lsda, isk
  USE fixed_occ,     ONLY : tfixed_occ
  USE qpoint,        ONLY : xq
  USE output,        ONLY : fildvscf, fildrho
  USE start_k,       ONLY : nk1, nk2, nk3
  USE disp,          ONLY : nq1, nq2, nq3
  USE input,         ONLY : delta_smear, nsmear, dis_win_min, dis_win_max, wannierize, &
                            ngaussw, dvscf_dir, bands_skipped, wdata, kmaps, ntempxx,  &
                            num_iter, dis_froz_max, fsthick, dis_froz_min, eig_read,   &
                            vme, degaussw, epexst, epwwrite, epbread, phonselfen, nqc2,&
                            elecselfen, a2f, plselfen, specfun_pl, nest_fn, filukk,    &
                            rand_nk, rand_k, rand_nq, rand_q, nkc1, nkc2, nkc3, nqc1,  &
                            nqc3, nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, eps_acoustic, nw, &
                            wmax, wmin, mp_mesh_q, mp_mesh_k, filqf, filkf, nswi, nc,  &
                            delta_qsmear, degaussq, band_plot, ephwrite, nstemp,       &
                            broyden_beta, conv_thr_raxis, temps, tempsmin, tempsmax,   &
                            broyden_ndim, wscut, nqstep, limag, lreal, muc,            &
                            gap_edge, conv_thr_iaxis, nqsmear, iprint, wepexst,        &
                            epwread, eliashberg, imag_read, kerread, kerwrite,         &
                            fermi_energy, efermi_read, max_memlt, fila2f, sigma_ref,   &
                            ep_coupling, nw_specfun, wmax_specfun, wmin_specfun,       &
                            laniso, lpolar, lifc, asr_typ, lscreen, scr_typ, nbndsub,  &
                            fermi_diff, smear_rpa, lcumulant, bnd_cum, proj, write_wfn,&
                            iswitch, liso, lacon, lpade, etf_mem, epbwrite,            &
                            nsiter, conv_thr_racon, npade, specfun_el, specfun_ph,     &
                            system_2d, delta_approx, title, int_mob, scissor,          &
                            iterative_bte, scattering, selecqread, epmatkqread,        &
                            ncarrier, carrier, scattering_serta, restart, restart_step,&
                            scattering_0rta, longrange_only, shortrange, scatread,     &
                            restart_filq, prtgkk, nel, meff, epsiheg, lphase, use_ws,  &
                            omegamin, omegamax, omegastep, lindabs, mob_maxiter,       &
                            auto_projections, scdm_proj, scdm_entanglement, scdm_mu,   &
                            scdm_sigma, assume_metal, wannier_plot, wannier_plot_list, &
                            wannier_plot_supercell, wannier_plot_scale, reduce_unk,    &
                            wannier_plot_radius, fermi_plot, fixsym, epw_no_t_rev,     &
                            epw_tr, epw_nosym, epw_noinv, epw_crysym, mob_maxfreq,     &
                            bfieldx, bfieldy, bfieldz, ii_eda, ii_partion,             &
                            ii_g, ii_charge, ii_n, ii_scattering, ii_only, ii_lscreen, &
                            gb_scattering, gb_only, gb_size,                           &
                            lwfpt, compute_dmat, ahc_nbnd, ahc_nbndskip,               &
                            ahc_win_min, ahc_win_max, mob_nfreq, plrn, restart_plrn,   &
                            conv_thr_plrn, end_band_plrn, init_sigma_plrn,             &
                            cal_psir_plrn, start_band_plrn,  type_plrn, nstate_plrn,   &
                            interp_Ank_plrn, interp_Bqu_plrn, init_k0_plrn,            &
                            full_diagon_plrn, mixing_Plrn, init_plrn, niter_plrn,      &
                            nDOS_plrn, edos_max_plrn, edos_min_plrn, edos_sigma_plrn,  &
                            pdos_sigma_plrn, pdos_max_plrn, pdos_min_plrn,             &
                            seed_plrn, ethrdg_plrn, time_rev_A_plrn, nhblock_plrn,     &
                            beta_plrn, Mmn_plrn, recal_Mmn_plrn, r0_plrn, debug_plrn,  &
                            time_rev_U_plrn,  g_start_band_plrn, g_end_band_plrn,      &
                            g_start_energy_plrn, g_end_energy_plrn, lrot,              &
                            model_vertex_plrn, model_enband_plrn, model_phfreq_plrn,   &
                            kappa_plrn, omega_LO_plrn, m_eff_plrn, step_wf_grid_plrn,  &
                            g_power_order_plrn, g_tol_plrn, io_lvl_plrn,               &
                            scell_mat_plrn, scell_mat, init_ntau_plrn, nethrdg_plrn,   &
                            adapt_ethrdg_plrn, init_ethrdg_plrn, gridsamp, griddens,   &
                            tc_linear, tc_linear_solver, fbw, dos_del, muchem,         &
                            filirobj, eps_cut_ir, icoulomb, emax_coulomb, emin_coulomb,&
                            filnscf_coul, positive_matsu,                              &
                            len_mesh, meshnum, wf_quasi, do_CHBB, nq_init, a_gap0,     &
                            start_mesh, DW, loptabs, mode_res, QD_min, QD_bin,         &
                            david_ndim_plrn, elecselfen_type, calc_nelec_wann,         &
                            lopt_w2b, epw_memdist, lfast_kmesh
  USE input,         ONLY : xk_all, xk_loc, xk_cryst, isk_all, isk_loc, et_all, et_loc
  USE global_var,    ONLY : elph, num_wannier_plot, wanplotlist, gtemp, qrpl
  USE ep_constants,  ONLY : ryd2mev, ryd2ev, ev2cmm1, kelvin2eV, zero, eps20, ang2m, one, bohr2nm
  USE constants,     ONLY : electron_si, AMU_RY, eps16
  USE io_files,      ONLY : tmp_dir, prefix
  USE control_flags, ONLY : iverbosity, modenum, gamma_only
  USE ions_base,     ONLY : amass
  USE mp_world,      ONLY : world_comm, nproc
  USE partial,       ONLY : atomo, nat_todo
  USE mp_pools,      ONLY : my_pool_id, me_pool, npool
  USE mp_images,     ONLY : nimage
  USE io_global,     ONLY : meta_ionode, meta_ionode_id, qestdin, stdout
  USE io_var,        ONLY : iunkf, iunqf
  USE paw_variables, ONLY : okpaw
  USE io,            ONLY : param_get_range_vector
  USE noncollin_module, ONLY : noncolin
  USE open_close_input_file, ONLY : open_input_file, close_input_file
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst
  !! Find if a file exists.
  LOGICAL, EXTERNAL :: imatches
  !! Does the title match
  CHARACTER(LEN = 256) :: outdir
  !! Output directory
  CHARACTER(LEN = 512) :: line
  !! Line in input file
  CHARACTER(LEN = 256), EXTERNAL :: trimcheck
  !! Trim the file name
  INTEGER :: ios
  !! INTEGER variable for I/O control
  INTEGER :: ios2
  !! INTEGER variable for I/O control
  INTEGER :: na
  !! counter on polarizations
  INTEGER :: it
  !! counter on iterations
  INTEGER :: modenum_aux
  !! auxilary variable for saving the modenum
  INTEGER :: i
  !! Counter for loops
  INTEGER :: ik
  !! Counter on k-points
  INTEGER :: itemp
  !! counter on temperatures
  INTEGER :: nstemp_hold = 0
  !! placeholder for nstemp
  INTEGER :: nk1tmp
  !! temp vars for saving kgrid info
  INTEGER :: nk2tmp
  !! temp vars for saving kgrid info
  INTEGER :: nk3tmp
  !! temp vars for saving kgrid info
  INTEGER :: ierr
  !! Error status
  REAL(kind = DP) :: b_abs
  !! Absolute magnetic field
  !
  NAMELIST / inputepw / &
       amass, outdir, prefix, iverbosity, fildvscf, rand_q, rand_nq, rand_k,   &
       elph, nq1, nq2, nq3, nk1, nk2, nk3, nbndsub, rand_nk, specfun_pl,       &
       filukk, epbread, epbwrite, epwread, epwwrite, etf_mem,                  &
       eig_read, wepexst, epexst, vme, elecselfen, phonselfen, use_ws, nc,     &
       degaussw, fsthick, nsmear, delta_smear, nqf1, nqf2, nqf3, nkf1, nkf2,   &
       dvscf_dir, ngaussw, epmatkqread, selecqread, nkf3, mp_mesh_k, mp_mesh_q,&
       wannierize, dis_win_max, dis_win_min, dis_froz_min, dis_froz_max, nswi, &
       num_iter, proj, bands_skipped, wdata, iprint, write_wfn, ephwrite,      &
       wmin, wmax, nw, eps_acoustic, a2f, nest_fn, plselfen, filqf, filkf,     &
       band_plot, fermi_plot, degaussq, delta_qsmear, nqsmear, nqstep,         &
       broyden_beta, broyden_ndim, nstemp, temps, bfieldx, bfieldy, bfieldz,   &
       conv_thr_raxis, conv_thr_iaxis, conv_thr_racon, wscut, system_2d,       &
       gap_edge, nsiter, muc, lreal, limag, lpade, lacon, liso, laniso, lpolar,&
       npade, lscreen, scr_typ, fermi_diff, smear_rpa, lcumulant, bnd_cum,     &
       lifc, asr_typ, kerwrite, kerread, imag_read, eliashberg,                &
       ep_coupling, fila2f, max_memlt, efermi_read, fermi_energy,              &
       specfun_el, specfun_ph, wmin_specfun, wmax_specfun, nw_specfun,         &
       delta_approx, scattering, int_mob, scissor, ncarrier, carrier,          &
       iterative_bte, scattering_serta, scattering_0rta, longrange_only,       &
       scatread, restart, restart_step, restart_filq, prtgkk, nel, meff,       &
       epsiheg, lphase, omegamin, omegamax, omegastep, lindabs, sigma_ref,     &
       mob_maxiter, auto_projections, scdm_proj, scdm_entanglement, scdm_mu,   &
       scdm_sigma, assume_metal, wannier_plot, wannier_plot_list, reduce_unk,  &
       wannier_plot_supercell, wannier_plot_scale, wannier_plot_radius,        &
       fixsym, epw_no_t_rev, epw_tr, epw_nosym, epw_noinv, epw_crysym,         &
       mob_maxfreq, mob_nfreq, lrot, lwfpt, shortrange,                        &
       ii_g, ii_charge, ii_n, ii_scattering, ii_only, ii_lscreen, ii_eda,      &
       ii_partion, plrn, restart_plrn, conv_thr_plrn, end_band_plrn,           &
       gb_scattering, gb_only, gb_size,                                        &
       cal_psir_plrn, start_band_plrn,  type_plrn, nstate_plrn,                &
       interp_Ank_plrn, interp_Bqu_plrn, init_sigma_plrn, init_k0_plrn,        &
       full_diagon_plrn, mixing_Plrn, init_plrn, niter_plrn,                   &
       nDOS_plrn, edos_max_plrn, edos_min_plrn, edos_sigma_plrn,               &
       pdos_sigma_plrn, pdos_max_plrn, pdos_min_plrn,                          &
       seed_plrn, ethrdg_plrn, time_rev_A_plrn, nhblock_plrn,                  &
       beta_plrn, Mmn_plrn, recal_Mmn_plrn, r0_plrn, debug_plrn,               &
       time_rev_U_plrn,  g_start_band_plrn, g_end_band_plrn,                   &
       g_start_energy_plrn, g_end_energy_plrn,                                 &
       model_vertex_plrn, model_enband_plrn, model_phfreq_plrn,                &
       kappa_plrn, omega_LO_plrn, m_eff_plrn, step_wf_grid_plrn,               &
       g_power_order_plrn, g_tol_plrn, io_lvl_plrn,                            &
       scell_mat_plrn, scell_mat, init_ntau_plrn,                              &
       adapt_ethrdg_plrn, init_ethrdg_plrn, nethrdg_plrn, david_ndim_plrn,     &
       tc_linear, tc_linear_solver, gridsamp, griddens, fbw, dos_del, muchem,  &
       filirobj, eps_cut_ir, icoulomb, emax_coulomb, emin_coulomb,             &
       filnscf_coul, positive_matsu,                                           &
       loptabs, len_mesh, meshnum, wf_quasi, nq_init, start_mesh, DW,          &
       mode_res, QD_min, QD_bin, do_CHBB, lwfpt, compute_dmat, ahc_nbnd,       &
       ahc_nbndskip, ahc_win_min, ahc_win_max, a_gap0, elecselfen_type,        &
       calc_nelec_wann, lopt_w2b, epw_memdist
  ! --------------------------------------------------------------------------------
  !
  ! amass    : atomic masses
  ! iverbosity : verbosity control
  ! outdir   : directory where input,     output, temporary files reside
  ! elph     : if true calculate electron-phonon coefficients
  ! prefix   : the prefix of files produced by pwscf
  ! fildvscf : output file containing deltavsc
  ! fildrho  : output file containing deltarho
  !
  ! Added by Feliciano Giustino
  ! ngaussw    : smearing type after wann interp
  !              (n >= 0) : Methfessel-Paxton case. See PRB 40, 3616 (1989)
  !              (n=-1)   : Cold smearing See PRL 82, 3296 (1999)
  !              (n=-99)  : Fermi-Dirac case: 1.0/(1.0+exp(-x)).
  !                         If n = -99 you probably want assume_metal == .true. as well.
  ! degaussw   : corresponding width (units of eV)
  ! filqf      : file with fine q kmesh for interpolation
  ! filkf      : file with fine kmesh for interpolation
  ! filukk     : file with rotation matrix U(k) for interpolation
  ! tphases    : if true set absolute unitary gauge for eigenvectors
  ! epstrict   : if true use strict selection rule for phonon linewidht calculation
  ! fsthick    : the thickness of the Fermi shell for averaging the e-ph matrix elements (units of eV)
  ! eptemp     : temperature for the electronic Fermi occupations in the e-p calculation (units of Kelvin)
  ! fildvscf0  : file containing deltavscf to be used as fake perturbation to set phases
  ! nw         : nr. of bins for frequency scan in \delta( e_k - e_k+q - w ) (units of eV)
  ! wmin       : min frequency for frequency scan in \delta( e_k - e_k+q - w ) (units of eV)
  ! wmax       : max    "  "  "                                    (units of eV)
  ! nbndsub    : number of bands in the optimal subspace (when disentanglement is used)
  ! elecselfen : if .TRUE. calculate imaginary part of electron selfenergy due to e-p interaction
  ! phonselfen : if .TRUE. calculate imaginary part of phonon selfenergy due to e-p interaction
  ! dvscf_dir  : the dir containing all the .dvscf and .dyn files
  ! epbread    : read epmatq array from .epb files
  ! epbwrite   : write epmatq array to .epb files
  ! epwread    : read all quantities in Wannier representation from file epwdata.fmt
  ! epwwrite   : write all quantities in Wannier representation to file epwdata.fmt
  !
  ! Added by jn
  ! wannierize   : if .TRUE. run the wannier90 code to maximally localize the WFs
  ! dis_win_min  : lower bound on wannier90 disentanglement window
  ! dis_win_max  : upper bound on wannier90 disentanglement window
  ! dis_froz_min : lower bound on frozen wannier90 disentanglement window
  ! dis_froz_max : upper bound on frozen wannier90 disentanglement window
  ! num_iter     : number of iterations used in the wannier90 minimisation
  ! proj         : initial projections (states) of the wannier functions before minimization
  ! auto_prjections: if .TRUE. automatically generate initial projections for W90
  ! scdm_proj    : if .TRUE. calculate MLWFs without an initial guess via the SCDM algorithm
  ! scdm_entanglement : disentanglement type in the SCDM algorithm
  ! scdm_mu      : parameter for Wannier functions via SCDM algorithm
  ! scdm_sigma   : parameter for Wannier functions via SCDM algorithm
  ! bands_skipped: k-point independent list of bands excluded from the calculation of overlap and projection matrices in W90
  ! wdata        : Empty array that can be used to pass extra info to prefix.win file, for things not explicitly declared here
  ! iprint       : verbosity of the wannier90 code
  ! write_wfn    : writes out UNK files from pwscf run for plotting of XSF files
  ! kmaps        : if true, read kmap and kgmap from disk (prior run)
  ! eig_read     : if .TRUE. then readin a set of electronic eigenvalues in eV to replace the calcualted ones
  ! wepexst      : if .TRUE. prefix.epmatwe files are already on disk. don't recalculate. debugging param
  ! epexst       : if .TRUE. prefix.epmatwp files are already on disk. don't recalculate  debugging param
  ! nest_fn      : if .TRUE., calculate the nesting function for a given set of q's
  ! nsmear       : number of smearing values to use for the selfen_phon call
  ! delta_smear  : change in energy for each additional nsmear ( units of eV)
  !
  ! Added by Roxana Margine
  ! ephwrite    : if .TRUE. write el-phonon matrix elements on the fine mesh to file
  ! eps_acoustic: min phonon frequency for e-p and a2f calculations (units of cm-1)
  ! band_plot   : if .TRUE. write files to plot band structure and phonon dispersion
  ! fermi_plot  : if .TRUE. write files to plot Fermi surface
  ! degaussq    : smearing for sum over q in e-ph coupling (units of meV)
  ! delta_qsmear: change in energy for each additional smearing in the a2f (units of meV)
  ! nqsmear     : number of smearings used to calculate a2f
  ! nqstep      : number of bins for frequency used to calculate a2f
  ! wscut       : upper limit for frequency integration in Eliashberg equations (at least 5 times wsphmax) (units of eV)
  ! broyden_beta : mixing factor for broyden mixing
  ! broyden_ndim : number of iterations used in mixing scheme
  ! nstemp      : number of temperature points for which the Eliashberg equations are solved
  ! tempsmin    : minimum temperature for which the Eliashberg equations are solved
  ! tempsmax    : maximum temperature " "
  ! conv_thr_raxis : convergence threshold for iterative solution of real-axis Eliashberg equations
  ! conv_thr_iaxis : convergence threshold for iterative solution of imag-axis Eliashberg equations
  ! conv_thr_racon : convergence threshold for iterative solution of analytic continuation of
  !                  Eliashberg equations from imag- to real-axis
  ! gap_edge     : initial guess of the superconducting gap (in eV)
  ! nsiter       : nr of iterations for self-consitency cycle
  ! npade        : percentage of Matsubara points used in Pade continuation
  ! muc          : effective Coulomb potential
  ! lreal        : if .TRUE. solve the real-axis Eliashberg eqautions
  ! limag        : if .TRUE. solve the imag-axis Eliashberg eqautions
  ! lpade        : if .TRUE. use pade approximants to continue imag-axis
  !                Eliashberg equtions to real-axis
  ! lacon        : if .TRUE. use analytic continuation to continue imag-axis
  !                Eliashberg equtions to real-axis
  ! liso         : if .TRUE. solve isotropic case
  ! laniso       : if .TRUE. solve anisotropic case
  ! kerwrite     : if .TRUE. write kp and km to files .ker for real-axis calculations
  ! kerread      : if .TRUE. read kp and km from files .ker for real-axis calculations
  ! imag_read    : if .TRUE. read from files Delta and Znorm on the imaginary-axis
  ! eliashberg   : if .TRUE. solve the Eliashberg equations
  ! ep_coupling  : if .TRUE. run e-p coupling calculation
  ! fila2f       : input file with eliashberg spectral function
  ! max_memlt    : maximum memory that can be allocated per pool
  ! efermi_read  : if. true. read from input file
  ! fermi_energy : fermi eneergy read from input file (units of eV)
  ! wmin_specfun : min frequency in electron spectral function due to e-p interaction (units of eV)
  ! wmax_specfun : max frequency in electron spectral function due to e-p interaction (units of eV)
  ! nw_specfun   : nr. of bins for frequency in electron spectral function due to e-p interaction
  ! delta_approx : if .TRUE. the double delta approximation is used to compute the phonon self-energy
  !
  ! Added by Samad Hajinazar
  ! tc_linear        : if .TRUE. linearized Eliashberg eqn. for Tc will be solved
  ! tc_linear_solver : Algorithm to solve eigenvalue problem for Tc (default='power', 'lapack')
  ! gridsamp         : Type of the Matsubara freq. sampling 
  !                    (-1= read from file; 0= uniform; 1= sparse; 2= sparse-ir, 3= uniform (FFT))
  ! griddens         : Measure of sparsity of the grid (default=1.d0, larger values give denser mesh)
  ! fbw              : if .TRUE. full-bandwidth calculations will be performed
  ! dos_del          : Delta_E in electronic dos for Fermi window (in eV)
  ! muchem           : if .TRUE. chem. pot. is updated in fbw calculations
  !
  ! Added by Hitoshi Mori
  ! a_gap0           : This determines the shape of initial guess of gap function
  !                    a_gap0 = negative, step function will be used.
  !                    a_gap0 = 0.0, we will use an initial guess with no frequency-dependence.
  !                    a_gap0 > eps8, we will use the Lorentzian: f(iw) = gap0 / (1 + a_gap0 * (iw / wsphmax)**2).
  ! positive_matsu   : If true, the domain of Matsubara frequency is limited to positive. (default=true) 
  ! eps_cut_ir       : This is a threshold to ignore negligibly small IR coefficients (default=1d-5).
  !                    If eps_cut_ir is zero or negative, all coefficients are used in iterative calculations.             
  ! filirobj         : input file with the objects of IR-basis
  ! icoulomb         : this tag specifies how to calculate the Coulomb contribution to the Eliashberg eqs.
  ! emax_coulomb     : Only the bands lower than "emax_coulomb + efermi" are considered
  ! emin_coulomb     : Only the bands upper than "emin_coulomb + efermi" are considered
  ! filnscf_coul     : input file with electronic eigenvalues for Coulomb contribution
  !
  ! Added by Carla Verdi & Samuel Pon\'e
  ! lpolar     : if .TRUE. enable the correct Wannier interpolation in the case of polar material.
  ! lifc       : if .TRUE. reads interatomic force constants produced by q2r.x for phonon interpolation
  ! asr_typ    : select type of ASR if lifc=.TRUE. (as in matdyn); otherwise it is the usual simple sum rule
  ! lscreen    : if .TRUE. the e-ph matrix elements are screened by the RPA or TF dielectric function
  ! scr_typ    : if 0 calculates the Lindhard screening, if 1 the Thomas-Fermi screening
  ! fermi_diff : difference between Fermi energy and band edge (in eV)
  ! smear_rpa  : smearing for the calculation of the Lindhard function (in eV)
  ! lcumulant  : if .TRUE. calculates the electron spectral function using the cumulant expansion method
  !              (can be used as independent postprocessing by setting ep_coupling=.FALSE.)
  ! bnd_cum    : band index for which the cumulant calculation is done
  !              (for more than one band, perform multiple calculations and add the results together)
  !
  ! Added by Samuel Ponc\'e
  ! specfun_el      : if .TRUE. calculate electron spectral function due to e-p interaction
  ! specfun_ph      : if .TRUE. calculate phonon spectral function due to e-p interaction
  ! specfun_pl      : if .TRUE. calculate plason spectral function
  ! restart         : if .TRUE. a run can be restarted from the interpolation level
  ! restart_step    : Create a restart point every restart_step q/k-points
  ! restart_filq    : Use to merge different q-grid scattering rates (name of the file)
  ! scattering      : if .TRUE. scattering rates are calculated
  ! scattering_serta: if .TRUE. scattering rates are calculated using self-energy relaxation-time-approx
  ! scatread        : if .TRUE. the current scattering rate file is read from file.
  ! scattering_0rta : if .TRUE. scattering rates are calculated using 0th order relaxation-time-approx
  ! int_mob         : if .TRUE. computes the intrinsic mobilities. This means that the
  !                   electron and hole carrier density is equal.
  ! iterative_bte   : if .TRUE. computes the iterative solution to the BTE. Need a
  !                   prior run with ERTA.
  ! scissor         : Value of the scissor shitf in eV. This only affects the CBM of etf. Do you use in
  !                   metals obviously.
  ! carrier         : if .TRUE. computes the doped carrier mobilities.
  ! ncarrier        : Set the Fermi level so that the carrier concentration is
  !                   " ncarrier". If ncarrier > 0, electron doping, hole doping otherwise
  ! longrange_only  : if .TRUE. computes the long-range part of the el-ph (can
  !                   only be used with lpolar = .TRUE. )
  ! shortrange      : if .TRUE. computes the short-range part of the el-ph (can
  !                   only be used with lpolar = .TRUE. )
  ! prtgkk          : Print the vertex |g| [meV]. This generates huge outputs.
  ! etf_mem         : if 0 no optimization, if 1 less memory is used for the fine grid interpolation
  !                   When etf_mem == 2, an additional loop is done on mode for the fine grid interpolation
  !                   part. This reduces the memory further by a factor "nmodes". [This is now depreciated]
  ! plselfen        : Calculate the electron-plasmon self-energy.
  ! nel             : Fractional number of electrons in the unit cell
  ! meff            : Density of state effective mass (in unit of the electron mass)
  ! epsiheg         : Dielectric constant at zero doping
  ! lphase          : If .TRUE., fix the gauge on the phonon eigenvectors and electronic eigenvectors - DS
  ! mob_maxiter     : Maximum number of iteration for the IBTE.
  ! use_ws          : If .TRUE., use the Wannier-center to create the Wigner-Seitz cell.
  ! epmatkqread     : If .TRUE., restart an IBTE calculation from scattering written to files.
  ! selecqread      : If .TRUE., restart from the selecq.fmt file
  ! nc              : Number of carrier for the Ziman resistivity formula (can be fractional)
  ! bfieldx, y, z   : Value of the magnetic field in Tesla along x, y, z direction.
  ! system_2d       : if 'no' then 3D bulk materials [default]
  !                   if 'gaussian' then the long-range terms include dipoles only and the range separation
  !                   function is approximated by a Gaussian following Ref. Phys. Rev. B 94, 085415 (2016)
  !                   if 'dipole_sp' then the long-range terms include dipoles following PRB 107, 155424 (2023)
  !                   if 'quadrupole' then the long-range terms include dipoles and quadrupoles terms
  !                   following PRL 130, 166301 (2023) and requires the presence of a "quadrupole.fmt" file.
  !                   if 'dipole_sh' then the long-range terms include dipoles following [PRB 105, 115414 (2022)]
  !
  ! Addad by Zhe Liu
  ! calc_nelec_wann : compute number of electrons in wannierized band
  ! epw_memdist     : distributed storage of epmatwp in MPI processes, only works with etf_mem = 0
  !
  ! Added by Manos Kioupakis
  ! omegamin        : Photon energy minimum
  ! omegamax        : Photon energy maximum
  ! omegastep       : Photon energy step in evaluating phonon-assisted absorption spectra (in eV)
  ! sigma_ref       : Reference conductivity for resistive contribution of FCA
  ! lindabs         : If .TRUE., do phonon-assisted absorption
  !
  ! Added by Felix Goudreault
  ! assume_metal     : If .TRUE. => we are dealing with a metal
  !
  ! Added by Hyungjun Lee
  ! wannier_plot           : If .TRUE., plot Wannier functions
  ! wannier_plot_list      : Field read for parsing Wannier function list
  ! wannier_plot_supercell : Size of supercell for plotting Wannier functions
  ! wannier_plot_scale     : Scaling parameter for cube files
  ! wannier_plot_radius    : Cut-off radius for plotting Wannier functions
  ! reduce_unk             : If .TRUE., plot Wannier functions on reduced grids
  !
  ! Added by Samuel Ponc\'e, Hyungjun Lee and Roxana Margine
  ! vme : if 'dipole' then computes the velocity as dipole+commutator = <\psi_mk|p/m+i[V_NL,r]|\psi_nk>
  !     : if 'wannier' then computes the velocity as dH_nmk/dk - i(e_nk-e_mk)A_nmk where A is the Berry connection
  !     : Note: Before v5.4, vme = .FALSE. was the velocity in the local approximation as <\psi_mk|p|\psi_nk>
  !             Before v5.4, vme = .TRUE. was = to 'wannier'
  !
  ! Added by Jae-Mo Lihm for Wannier function perturbation theory
  ! lwfpt : enable Wannier function perturbation theory calculations
  ! compute_dmat: compute dmat, the overlap matrix between psi(Sk) and S*psi(k)
  ! ahc_nbnd: Number of bands included in ph.x Allen-Heine-Cardona calculation
  ! ahc_nbndskip: Number of low-lying bands excluded in ph.x AHC calculation
  ! ahc_win_min : Lower bound of AHC window for the lower Fan term.
  ! ahc_win_max : Upper bound of AHC window for the lower Fan term.
  ! Added by Samuel Ponc\'e for adiabatic support of AHC
  ! elecselfen_type : if 'adiabatic' then computes the adiabatic electron self-energy.
  !     Only to be used in non-IR-active materials. See J. Chem. Phys. 143, 102813 (2015)
  !     for more information.
  !     If 'nonadiabatic' then computes the non-adiabatic electron self-energy (default).
  ! lopt_w2b : if .TRUE. use optimized version of Wannier-to-Bloch Fourier transformation
  !
  IF ( npool * nimage /= nproc ) THEN
    CALL errore("readin", "Number of processes must be equal to product "//&
                "of number of pools and"//new_line('n')//"     number of "//&
                "images"//new_line('n')//"     Image parallelization can be used "//&
                "only in calculations on coarse grid.", 1)
  END IF
  !
  nk1tmp = 0
  nk2tmp = 0
  nk3tmp = 0
  !
  IF (meta_ionode) THEN
    !
    ! ... Input from file (ios=0) or standard input (ios=-1) on unit "qestdin"
    !
    ios = open_input_file (  )
    !
    ! ... Read the first line of the input file
    !
    IF ( ios <= 0 ) READ(qestdin, '(A)', IOSTAT = ios) title
    !
  ENDIF
  !
  CALL mp_bcast(ios, meta_ionode_id, world_comm)
  CALL errore('readin', 'reading input file ', ABS(ios))
  CALL mp_bcast(title, meta_ionode_id, world_comm)
  !
  ! Rewind the input if the title is actually the beginning of inputph namelist
  !
  IF (imatches("&inputepw", title)) THEN
    WRITE(stdout, '(6x,a)') "Title line not specified: using 'default'."
    title = 'default'
    IF (meta_ionode) REWIND(qestdin, IOSTAT = ios)
    CALL mp_bcast(ios, meta_ionode_id, world_comm  )
    CALL errore('readin', 'Title line missing from input.', ABS(ios))
  ENDIF
  !
  ! Set default values for variables in namelist
  amass(:)     = 0.d0
  iverbosity   = 0
  elph         = .FALSE.
  elecselfen   = .FALSE.
  phonselfen   = .FALSE.
  plselfen     = .FALSE.
  specfun_el   = .FALSE.
  specfun_ph   = .FALSE.
  specfun_pl   = .FALSE.
  epbread      = .FALSE.
  epbwrite     = .FALSE.
  epwread      = .FALSE.
  epwwrite     = .TRUE.
  restart      = .FALSE.
  restart_step = 100
  wannierize   = .FALSE.
  write_wfn    = .FALSE.
  kmaps        = .FALSE.
  nest_fn      = .FALSE.
  wepexst      = .FALSE.
  epexst       = .FALSE.
  eig_read     = .FALSE.
  dis_win_max  = 9999.d0
  dis_win_min  = -9999.d0
  dis_froz_max = 9999.d0
  dis_froz_min = -9999.d0
  num_iter     = 200
  proj(:)      = ''
  auto_projections = .FALSE.
  scdm_proj    = .FALSE.
  scdm_entanglement = 'isolated'
  scdm_mu      = 0.d0
  scdm_sigma   = 1.d0
  bands_skipped= ''
  wdata(:)     = ''
  iprint       = 2
  wannier_plot = .FALSE.
  wannier_plot_scale = 1.0d0
  wannier_plot_radius = 3.5d0
  wannier_plot_supercell = (/5,5,5/)
  wannier_plot_list = ''
  reduce_unk   = .FALSE.
  wmin         = 0.d0
  wmax         = 0.3d0
  eps_acoustic  = 5.d0 ! cm-1
  nw           = 10
  fsthick      = 1.d10 ! eV
  degaussw     = 0.025d0 ! eV
  a2f          = .FALSE.
  etf_mem      = 1
  ngaussw      = 1
  outdir       = '.'
  dvscf_dir    = '.'
  prefix       = 'pwscf'
  filqf        = ' '
  restart_filq = ' '
  filkf        = ' '
  fildrho      = ' '
  fildvscf     = ' '
  filukk       = ' '
  rand_q       = .FALSE.
  delta_approx = .FALSE.
  rand_nq      = 1
  rand_k       = .FALSE.
  rand_nk      = 1
  nq1          = 0
  nq2          = 0
  nq3          = 0
  nk1          = 0
  nk2          = 0
  nk3          = 0
  nqf1         = 0
  nqf2         = 0
  nqf3         = 0
  nkf1         = 0
  nkf2         = 0
  nkf3         = 0
  mp_mesh_k    = .FALSE.
  mp_mesh_q    = .FALSE.
  nbndsub      = 0
  nsmear       = 1
  delta_smear  = 0.01d0 ! eV
  modenum      = 0 ! Was -1 previously and read from Modules/input_parameters.f90
                   ! INTEGER :: modenum = 0. In QE 5, modenum variable does not exist
                   ! anymore. Change the default EPW value to match the previous QE one.
  vme          = 'wannier' ! Note: Was .FALSE. by default until EPW v5.1 and then .TRUE. until EPW v5.4
  elecselfen_type = 'nonadiabatic'
  ephwrite     = .FALSE.
  band_plot    = .FALSE.
  fermi_plot   = .FALSE.
  nqsmear      = 10
  nqstep       = 500
  delta_qsmear = 0.05d0 ! meV
  degaussq     = 0.05d0 ! meV
  lreal        = .FALSE.
  limag        = .FALSE.
  lpade        = .FALSE.
  lacon        = .FALSE.
  liso         = .FALSE.
  laniso       = .FALSE.
  lpolar       = .FALSE.
  lifc         = .FALSE.
  asr_typ      = 'simple'
  lscreen      = .FALSE.
  scr_typ      = 0
  fermi_diff   = 1.d0
  smear_rpa    = 0.05d0
  lcumulant    = .FALSE.
  bnd_cum      = 1
  kerwrite     = .FALSE.
  kerread      = .FALSE.
  imag_read    = .FALSE.
  eliashberg   = .FALSE.
  ep_coupling  = .TRUE.
  tc_linear    = .FALSE.
  tc_linear_solver = 'power'
  gridsamp     = 0
  griddens     = 1.d0
  fbw          = .FALSE.
  dos_del      = 1.d-02
  muchem       = .FALSE.
  a_gap0       = 1.0d0
  nswi         = 0
  wscut        = 0.d0
  broyden_beta = 0.7d0
  broyden_ndim = 8
  conv_thr_raxis = 5.d-04
  conv_thr_iaxis = 1.d-05
  conv_thr_racon = 5.d-04
  gap_edge     = 0.d0
  positive_matsu = .TRUE.
  icoulomb     = 0
  emax_coulomb = 1.d5
  emin_coulomb = -1.d5
  filirobj     = ' '
  filnscf_coul = ' '
  eps_cut_ir   = 1.d-5
  nstemp       = 0
  temps(:)     = -1.d0
  nsiter       = 40
  npade        = 90
  muc          = 0.d0
  fila2f       = ' '
  max_memlt    = 2.85d0
  efermi_read  = .FALSE.
  fermi_energy = 0.d0
  wmin_specfun = 0.d0 ! eV
  wmax_specfun = 0.3d0 ! eV
  nw_specfun   = 100
  system_2d    = 'no'   ! Previously was .FALSE.
  scattering   = .FALSE.
  scattering_serta = .FALSE.
  scatread     = .FALSE.
  scattering_0rta = .FALSE.
  int_mob      = .FALSE.
  iterative_bte = .FALSE.
  scissor      = 0.d0 ! eV
  carrier      = .FALSE.
  ncarrier     = 0.d0 ! cm^-3
  longrange_only = .FALSE.
  shortrange   = .FALSE.
  prtgkk       = .FALSE.
  nel          = 0.0d0
  meff         = 1.d0
  epsiheg      = 1.d0
  lphase       = .FALSE.
  lrot         = .FALSE.
  omegamin     = 0.d0  ! eV
  omegamax     = 10.d0 ! eV
  omegastep    = 1.d0  ! eV
  sigma_ref    = 1.d0  ! 1/(Ohm m)
  lindabs      = .FALSE.
  mob_maxiter  = 50
  use_ws       = .FALSE.
  epmatkqread  = .FALSE.
  selecqread   = .FALSE.
  nc           = 4.0d0
  assume_metal = .FALSE.  ! default is we deal with an insulator
  fixsym       = .FALSE.
  epw_no_t_rev = .TRUE.
  epw_tr       = .TRUE.
  epw_nosym    = .FALSE.
  epw_noinv    = .FALSE.
  epw_crysym   = .FALSE.
  bfieldx      = 0.d0  ! Tesla
  bfieldy      = 0.d0  ! Tesla
  bfieldz      = 0.d0  ! Tesla
  mob_maxfreq  = 100 ! Maximum frequency for spectral decomposition in meV
  mob_nfreq    = 100 ! Number of frequency for the spectral decomposition
  ii_g         = .FALSE.
  ii_charge    = 1.0d0
  ii_n         = 0.0d0
  ii_scattering = .FALSE.
  ii_only      = .FALSE.
  ii_lscreen   = .TRUE.
  ii_eda       = 0.0d0
  ii_partion   = .FALSE.
  gb_scattering = .FALSE.
  gb_only      = .FALSE.
  gb_size      = 0.0d0  ! nm
  !
  ! Added for polaron calculations by Chao Lian
  nstate_plrn  = 1
  niter_plrn   = 50
  plrn         = .FALSE.
  restart_plrn = .FALSE.
  model_vertex_plrn = .false.
  model_enband_plrn = .false.
  model_phfreq_plrn = .false.
  kappa_plrn    = 0.0
  omega_LO_plrn = 0.0
  m_eff_plrn    = 0.0
  conv_thr_plrn = 1E-5
  g_power_order_plrn = 1
  step_wf_grid_plrn  = 1
  cal_psir_plrn      = .FALSE.
  interp_Ank_plrn    = .FALSE.
  interp_Bqu_plrn    = .FALSE.

  start_band_plrn   = 0
  end_band_plrn     = 0
  g_start_band_plrn = 0
  g_end_band_plrn   = 0
  g_start_energy_plrn = -10.0
  g_end_energy_plrn   = 10.0

  full_diagon_plrn = .FALSE.
  mixing_Plrn      = 1.0
  init_plrn        = 1
  Mmn_plrn         = .FALSE.
  recal_Mmn_plrn   = .FALSE.
  debug_plrn       = .FALSE.
  r0_plrn          = 0.0
  nDOS_plrn        = 1000
  edos_min_plrn    = 0.0 ! eV
  pdos_min_plrn    = 0.0 ! meV
  edos_max_plrn    = 0.0 ! eV
  pdos_max_plrn    = 0.0 ! meV
  edos_sigma_plrn  = 0.01d0 ! eV
  pdos_sigma_plrn  = 0.1 ! meV
  type_plrn        = -1
  init_sigma_plrn  = 4.6
  init_k0_plrn     = (/1000.d0, 1000.d0, 1000.d0/)
  ethrdg_plrn      = 1E-6
  time_rev_A_plrn  = .FALSE.
  time_rev_U_plrn  = .FALSE.
  nhblock_plrn     = 1
  beta_plrn        = 0.0
  g_tol_plrn       = -0.01
  io_lvl_plrn      = 0
  scell_mat_plrn   = .FALSE.
  scell_mat(1, 1:3) = (/1, 0, 0/)
  scell_mat(2, 1:3) = (/0, 1, 0/)
  scell_mat(3, 1:3) = (/0, 0, 1/)
  init_ntau_plrn    = 1
  adapt_ethrdg_plrn = .FALSE.
  init_ethrdg_plrn  = 1.d-2
  nethrdg_plrn      = 11
  david_ndim_plrn   = 4
  loptabs     = .FALSE.
  do_CHBB     = .FALSE.
  len_mesh    = 1
  meshnum     = 1
  wf_quasi    = -1
  nq_init     = -1
  start_mesh  = 0
  DW          = 0
  mode_res    = 0
  QD_bin      = 0.0
  QD_min      = 0.005
  lwfpt       = .FALSE.
  compute_dmat = .FALSE.
  calc_nelec_wann = .FALSE.
  ahc_nbnd       = -1
  ahc_nbndskip   = 0
  ahc_win_min    = -9999.d0
  ahc_win_max    = -9999.d0
  lopt_w2b    = .FALSE.
  lfast_kmesh = .FALSE.
  epw_memdist = .FALSE.
  !
  ! Reading the namelist inputepw and check
  IF (meta_ionode) THEN
    READ(qestdin, inputepw, IOSTAT = ios)
    ios2 = 0
    IF (ios /= 0) THEN
      BACKSPACE(qestdin)
      READ(qestdin, '(A512)', IOSTAT = ios2) line
    ENDIF
    IF (ios2 /= 0) CALL errore('readin', 'Could not find namelist &inputepw', 2)
    IF (ios /= 0) THEN
      CALL errore('readin', 'Bad line in namelist &inputepw: "' &
                 & //TRIM(line)//'" (error could be in the previous line)', 1)
    ENDIF
    ios = close_input_file ( )
  ENDIF ! meta_ionode
  !
  IF (meta_ionode) THEN
    IF (wannier_plot) THEN
      IF (wannier_plot_radius < 0.0d0) &
        CALL errore('readin', 'Error: wannier_plot_radius must be positive', 1)
      IF (wannier_plot_scale < 0.0d0) &
        CALL errore('readin', 'Error: wannier_plot_scale must be positive', 1)
      IF (ANY(wannier_plot_supercell <= 0)) &
        CALL errore('readin', &
          'Error: Three positive integers must be explicitly provided &
                  for wannier_plot_supercell', 1)
      CALL param_get_range_vector(wannier_plot_list, num_wannier_plot, .TRUE.)
      IF (num_wannier_plot == 0) THEN
        num_wannier_plot = nbndsub
        ALLOCATE(wanplotlist(num_wannier_plot), STAT = ierr)
        IF (ierr /= 0) CALL errore('readin', 'Error allocating wanplotlist', 1)
        DO i = 1, num_wannier_plot
          wanplotlist(i) = i
        ENDDO
      ELSE
        ALLOCATE(wanplotlist(num_wannier_plot), STAT = ierr)
        IF (ierr /= 0) CALL errore('readin', 'Error allocating wanplotlist', 1)
        CALL param_get_range_vector(wannier_plot_list, num_wannier_plot, .FALSE., wanplotlist)
        IF (ANY(wanplotlist < 1) .OR. ANY(wanplotlist > nbndsub)) &
          CALL errore('readin', &
            'Error: wannier_plot_list asks for a non-valid wannier function to be plotted', 1)
      ENDIF
    ENDIF
  ENDIF ! meta_ionode
  CALL mp_bcast(wannier_plot, meta_ionode_id, world_comm)
  IF (wannier_plot) CALL mp_bcast(num_wannier_plot, meta_ionode_id, world_comm)
  IF ((wannier_plot) .AND. (.NOT. meta_ionode)) THEN
    ALLOCATE(wanplotlist(num_wannier_plot), STAT = ierr)
    IF (ierr /= 0) CALL errore('readin', 'Error allocating wanplotlist', 1)
  ENDIF
  IF (wannier_plot) CALL mp_bcast(wanplotlist, meta_ionode_id, world_comm)
  !
  nk1tmp = nk1
  nk2tmp = nk2
  nk3tmp = nk3
  !
  ! Explaination: nk? and nq? are used by QE modules and therefore needs to be define
  !               We define a EPW coarse grid nkc? and nqc? which is the same as nk? and nq?
  !               but internal to EPW.
  nkc1 = nk1
  nkc2 = nk2
  nkc3 = nk3
  nqc1 = nq1
  nqc2 = nq2
  nqc3 = nq3
  !
  ! Check all namelist variables
  !
  ! file with rotation matrix U(k) for interpolation
  filukk = TRIM(prefix) // '.ukk'
  IF (nsmear < 1) CALL errore('readin', 'Wrong number of nsmears', 1)
  IF (iverbosity < 0 .OR. iverbosity > 4) CALL errore('readin', 'Wrong iverbosity', 1)
  IF (epbread .AND. epbwrite) CALL errore('readin', 'epbread cannot be used with epbwrite', 1)
  IF (epbread .AND. epwread) CALL errore('readin', 'epbread cannot be used with epwread', 1)
  IF (degaussw * 4.d0 > fsthick) CALL errore('readin', ' degaussw too close to fsthick', 1)
  IF ((nw < 1) .OR. (nw > 1000)) CALL errore ('readin', 'unreasonable nw', 1)
  IF (elecselfen .AND. plselfen) CALL errore('readin', &
      'Electron-plasmon self-energy cannot be computed with electron-phonon', 1)
  IF (phonselfen .AND. plselfen) CALL errore('readin', &
      'Electron-plasmon self-energy cannot be computed with electron-phonon', 1)
  IF (specfun_el .AND. plselfen) CALL errore('readin', &
      'Electron-plasmon self-energy cannot be computed with el-ph spectral function', 1)
  IF (specfun_ph .AND. plselfen) CALL errore('readin', &
      'Electron-plasmon self-energy cannot be computed with el-ph spectral function', 1)
  IF (elecselfen .AND. specfun_pl ) CALL errore('readin', &
      'Electron-plasmon spectral function cannot be computed with electron-phonon', 1)
  IF (phonselfen .AND. specfun_pl) CALL errore('readin', &
      'Electron-plasmon spectral function cannot be computed with electron-phonon', 1)
  IF (specfun_el .AND. specfun_pl) CALL errore('readin', &
      'Electron-plasmon spectral function cannot be computed with el-ph spectral function', 1)
  IF (specfun_ph .AND. specfun_pl) CALL errore('readin', &
      'Electron-plasmon spectral function cannot be computed with el-ph spectral function', 1)
  IF (ABS(degaussw) < eps16 .AND. specfun_el) CALL errore('readin', &
      'adapt_smearing cannot be used with spectral functions. Set degaussw > 0.', 1)
  IF (a2f .AND. .NOT. phonselfen) CALL errore('readin', 'a2f requires phonoselfen', 1)
  IF (elph .AND. .NOT. ep_coupling) CALL errore('readin', 'elph requires ep_coupling=.true.', 1)
  IF ((elph .AND. wannierize) .AND. (epwread)) CALL errore('readin', &
      'must use same w90 rotation matrix for entire run', 1)
  IF (wannierize .AND. .NOT. ep_coupling) CALL errore('readin', &
      'wannierize requires ep_coupling=.true.', 1)
  IF ((wmin > wmax)) CALL errore('readin', ' check wmin, wmax ', 1)
  IF ((wmin_specfun > wmax_specfun)) CALL errore('readin', 'check wmin_specfun, wmax_specfun', 1)
  IF ((nw_specfun < 2)) CALL errore('readin', 'nw_specfun must be at least 2', 1)
  IF ((nqstep < 2)) CALL errore('readin', 'nqstep must be at least 2', 1)
  IF ((nbndsub > 200)) CALL errore('readin', 'too many wannier functions increase size of projx', 1)
  IF ((phonselfen .OR. specfun_ph) .AND. mp_mesh_k) &
    CALL errore('readin', 'phonon self-energy only works with full uniform k-mesh', 1)
  IF ((elecselfen .OR. specfun_el) .AND. mp_mesh_q) &
    CALL errore('readin', 'electron self-energy only works with full uniform q-mesh', 1)
  IF (ephwrite) THEN
    IF (.NOT. ep_coupling .AND. .NOT. elph) CALL errore('readin', &
      'ephwrite requires ep_coupling=.TRUE., elph=.TRUE.', 1)
    IF (rand_k .OR. rand_q) CALL errore('readin', 'ephwrite requires a uniform grid', 1)
    IF (MOD(nkf1,nqf1) /= 0 .OR. MOD(nkf2,nqf2) /= 0 .OR. MOD(nkf3,nqf3) /= 0) &
    CALL errore('readin', 'ephwrite requires nkf1,nkf2,nkf3 to be multiple of nqf1,nqf2,nqf3', 1)
  ENDIF
  IF (band_plot .AND. filkf == ' ' .AND. filqf == ' ') CALL errore('readin', &
      'plot band structure and phonon dispersion requires k- and q-points read from filkf and filqf files', 1)
  IF (band_plot .AND. filkf /= ' ' .AND. (nkf1 > 0 .OR. nkf2 > 0 .OR. nkf3 > 0)) CALL errore('readin', &
      'You should define either filkf or nkf when band_plot = .true.', 1)
  IF (band_plot .AND. filqf /= ' ' .AND. (nqf1 > 0 .OR. nqf2 > 0 .OR. nqf3 > 0)) CALL errore('readin', &
     'You should define either filqf or nqf when band_plot = .true.', 1)
  IF (filkf /= ' ' .AND. .NOT. efermi_read) CALL errore('readin', &
     'WARNING: if k-points are along a line, then efermi_read=.true. and fermi_energy must be given in the input file', -1)
  IF (MAXVAL(temps(:)) < 0.d0 .AND. nstemp > 0) &
    CALL errore('readin', 'temps(:) must be specified if nstemp > 0', 1)
  IF (nstemp > ntempxx) &
    CALL errore('readin', 'Maximum value of nstemp that can be used is 50', 1)
  IF ((ABS(ncarrier) > 1E+5) .AND. .NOT. carrier) CALL errore('readin', &
      'carrier must be .TRUE. if you specify ncarrier.', 1)
  IF (carrier .AND. (ABS(ncarrier) < 1E+5))  CALL errore('readin', &
      'The absolute value of the doping carrier concentration must be larger than 1E5 cm^-3', 1)
  IF ((longrange_only .OR. shortrange) .AND. (.NOT. lpolar)) CALL errore('readin', &
      'Error: longrange_only or shortrange can only be true if lpolar is true as well.', 1)
  IF (longrange_only .AND. shortrange) CALL errore('readin',&
      'Error: longrange_only and shortrange cannot be both true.', 1)
  IF (.NOT. epwread .AND. .NOT. epwwrite) CALL errore('readin', &
      'Error: Either epwread or epwwrite needs to be true. ', 1)
  !
  IF (etf_mem == 2) THEN
    WRITE(stdout, '(/,5x,a)') 'Warning: the option etf_mem=2 &
        &is deprecated and switched to etf_mem = 1 automatically.'
    etf_mem = 1
  ENDIF
  IF (etf_mem /= 0 .AND. epw_memdist) CALL errore('epw_readin', 'Error: epw_memdist only work with etf_mem = 0.', 1)
  IF (etf_mem == 0 .AND. .NOT. epw_memdist) &
     WRITE(stdout, '(/,5x,a)') 'Suggestion: epw_memdist == .true. can reduce the memory usage when etf_mem == 0.'
  IF (etf_mem == 3) THEN
    etf_mem = 1
    lfast_kmesh = .TRUE.
     WRITE(stdout, '(/5x, a)') 'etf_mem = 3 results in the internal parameter &
         &lfast_kmesh switched on for efficient k mesh generation.'
  ENDIF
  !
  IF (lfast_kmesh) THEN
    IF (.NOT. mp_mesh_k) CALL errore('readin', 'When lfast_kmesh (etf_mem = 3), you have to use mp_mesh_k == .true.', 1)
    IF (.NOT. efermi_read) CALL errore('readin', 'When lfast_kmesh (etf_mem = 3), you have to use efermi_read == .true.', 1)
    IF (int_mob) CALL errore('readin', 'When lfast_kmesh (etf_mem = 3), you have to use int_mob == .false.', 1)
    IF (.NOT. carrier) CALL errore('readin', 'When lfast_kmesh (etf_mem = 3), you have to use carrier == .true.', 1)
    IF (phonselfen) CALL errore('readin', 'phonselfen is not implemented with lfast_kmesh (etf_mem = 3)', 1)
    IF (filkf /= ' ' .OR. filqf /= ' ' .OR. rand_k .OR. rand_q) THEN
      CALL errore('readin', 'lfast_kmesh (etf_mem = 3) requires homogeneous grids', 1)
    ENDIF
  ENDIF
  !
  IF (etf_mem > 1 .OR. etf_mem < 0) CALL errore('readin', 'etf_mem can only be 0 or 1.', 1)
  !
  IF (lwfpt .AND. (ahc_nbnd <= 0)) CALL errore('readin', &
      'Error: ahc_nbnd must be set if wfpt is used.', 1)
  IF (lwfpt .AND. (ahc_win_min < -9990.d0)) CALL errore('readin', &
      'Error: ahc_win_min must be set if WFPT is used.', 1)
  IF (lwfpt .AND. (ahc_win_max < -9990.d0)) CALL errore('readin', &
      'Error: ahc_win_max must be set if WFPT is used.', 1)
  IF (lwfpt .AND. (fsthick < 1.d8)) CALL errore('readin', &
      'Error: fsthick cannot be used with WFPT.', 1)
  IF (lwfpt) compute_dmat = .TRUE.
  ahc_win_min = ahc_win_min / ryd2ev
  ahc_win_max = ahc_win_max / ryd2ev
  !
  IF ((asr_typ /= 'simple') .AND. (asr_typ /= 'crystal') .AND. (asr_typ /= 'one-dim') .AND. &
      (asr_typ /= 'zero-dim')) THEN
    CALL errore('set_asr','invalid Acoustic Sum Rule:' // asr_typ, 1)
  ENDIF
  !
  IF (scdm_proj) THEN
    IF ((TRIM(scdm_entanglement) /= 'isolated') .AND. &
        (TRIM(scdm_entanglement) /= 'erfc')     .AND. &
        (TRIM(scdm_entanglement) /= 'gaussian')) THEN
      CALL errore('pw2wan90epw', 'Can not recognize the choice for scdm_entanglement. ' &
                                  // 'Valid options are: isolated, erfc and gaussian', 1)
    ENDIF
  ENDIF
  !
  ! Check consistency with lopt_w2b
  !
  IF (lopt_w2b) THEN
    IF (etf_mem /= 0) CALL errore('readin', &
      'Only etf_mem == 0 is supported with lopt_w2b = .true.', 1)
    IF (filqf /= ' ') CALL errore('readin', &
      'lopt_w2b = .true. works only for homogeneous grids', 1)
  ENDIF
  !
  !
  ! Make sure the files exists
  !
  IF (meta_ionode) THEN
    IF (filkf /= ' ') THEN
      OPEN(UNIT = iunkf, FILE = filkf, STATUS = 'old', FORM = 'formatted', ERR = 100, IOSTAT = ios)
100   CALL errore('readin', 'opening file ' // filkf, ABS(ios))
      CLOSE(iunkf)
    ENDIF
    IF (filqf /= ' ') THEN
      OPEN(UNIT = iunqf, FILE = filqf, STATUS = 'old', FORM = 'formatted', ERR = 101, IOSTAT = ios)
101   CALL errore('readin', 'opening file ' // filqf, ABS(ios))
      CLOSE(iunqf)
    ENDIF
  ENDIF ! meta_ionode
  IF (iterative_bte) THEN
    ! The fine grids have to be homogeneous and the same. Otherwise the populations can oscillate.
    IF (nkf1 /= nqf1 .OR. nkf2 /= nqf2 .OR. nkf3 /= nqf3) THEN
      CALL errore('readin', 'Error: the fine k-points and q-points grids have to be the same when doing IBTE.', 1)
    ENDIF
  ENDIF
  IF (auto_projections .AND. proj(1) /= ' ') CALL errore('readin', &
      'Cannot specify both auto_projections and projections block', 1)
  IF ((auto_projections .AND. .NOT. scdm_proj) .OR. (.NOT. auto_projections .AND. scdm_proj)) &
    CALL errore('readin', 'auto_projections require both scdm_proj=.true. and auto_projections=.true.', 1)
  IF (dis_win_min > -9999.d0 + eps16) THEN
    dis_win_min  = -9999.d0
    WRITE(stdout, '(/,5x,a)') 'WARNING: The specified dis_win_min is ignored.'
    WRITE(stdout, '(5x,a)') "         You should instead use bands_skipped = 'exclude_bands = ...'"
    WRITE(stdout, '(5x,a)') "         to control the lower bound of band manifold."
  ENDIF
  !
  ! 2D interpolation
  IF ( ALL((/system_2d /= 'no', system_2d /= 'gaussian', system_2d /= 'dipole_sp', &
            system_2d /= 'quadrupole', system_2d /= 'dipole_sh'/)) ) THEN
    CALL errore('readin', 'invalid value system_2d = "'//TRIM(system_2d)//'"', 1)
  ENDIF
  CALL mp_bcast(system_2d, meta_ionode_id, world_comm)
  IF (system_2d == 'quadrupole' .OR. system_2d =='no') THEN
    ! If quadrupole file exist, read it
    IF (meta_ionode) THEN
      INQUIRE(FILE = 'quadrupole.fmt', EXIST = exst)
    ENDIF
    CALL mp_bcast(exst, meta_ionode_id, world_comm)
    !
    qrpl = .FALSE.
    IF (exst) qrpl = .TRUE.
    IF (system_2d == 'quadrupole' .AND. .NOT. exst) CALL errore('readin', 'Error: the file "quadrupole.fmt" was not found.', 1)
  ENDIF ! system_2d
  !
  IF (lpolar .AND. elecselfen_type == 'adiabatic') THEN
    WRITE(stdout, '(5x,a)') "WARNING: Adiabatic AHC should not be used if lpolar == .true. because"
    WRITE(stdout, '(5x,a)') "         it yields a non-convergent integral. Use 'nonadiabatic' instead."
    WRITE(stdout, '(5x,a)') "         If you used this option intentionally, use it at your own risk."
    WRITE(stdout, '(5x,a)') "         See J. Chem. Phys. 143, 102813 (2015) for more information."
  ENDIF
  IF (elecselfen_type /= 'adiabatic' .AND. elecselfen_type /= 'nonadiabatic') THEN
    CALL errore('readin', "Error: elecselfen_type must be 'adiabatic' or 'nonadiabatic'", 1)
  ENDIF
  !
  IF (lfast_kmesh) THEN
    WRITE(stdout,'(5x,a)') 'WARNING: The use of lfast_kmesh has been tested and validated for cubic and hexagonal materials.'
    WRITE(stdout,'(5x,a)') '         For other materials, use with care and possibly use lfast_kmesh = .FALSE.'
  ENDIF
  !
  b_abs = ABS(bfieldx) + ABS(bfieldy) + ABS(bfieldz)
  IF (b_abs > eps20 .AND. (.NOT. mp_mesh_k)) THEN
    WRITE(stdout,'(5x,a)') 'WARNING: Finite magnetic field is much more stable with k-point symmetry: mp_mesh_k == .true.'
    CALL errore('readin', 'Error: Finite magnetic field only implemented with k-point symmetry: mp_mesh_k == .true.', 1)
  ENDIF
  IF (b_abs > eps20 .AND. (filkf /= ' ')) THEN
    CALL errore('readin', 'Error: Finite magnetic field only implemented with homogeneous k-grids', 1)
  ENDIF
  IF (b_abs > eps20 .AND. (filqf /= ' ')) THEN
    CALL errore('readin', 'Error: Finite magnetic field only implemented with homogeneous k-grids', 1)
  ENDIF
  IF (vme /= 'dipole' .AND. vme /= 'wannier') THEN
    CALL errore('readin', "Error: vme must be 'dipole' or 'wannier'", 1)
  ENDIF
  !
  ! setup temperature array
  DO itemp = 1, ntempxx
    IF (temps(itemp) >= 0.d0) THEN
      nstemp_hold = itemp
    ENDIF
  ENDDO
  !
  ! Case of nstemp > 0 but temps(:) = 0 is caught above
  IF (nstemp_hold == 0 .AND. nstemp == 0) THEN !default mode (nstemp_hold == 0 if temps(:) = 0)
    nstemp = 1
    temps(1) = 300
    WRITE(stdout, '(/,5x,a)') 'No temperature supplied. Setting temps(:) to 300 K.'
  ELSE IF (nstemp == 0 .OR. nstemp_hold == nstemp) THEN !list mode
    nstemp = nstemp_hold !catches if nstemp not supplied, no effect if it is
    WRITE(stdout, '(/,5x,a)') 'Reading supplied temperature list.'
  ELSE IF (nstemp_hold < nstemp .AND. nstemp_hold == 2) THEN !even spacing mode
    tempsmin = temps(1)
    tempsmax = temps(2)
    IF (tempsmin >= tempsmax) THEN !bad start and end points
      CALL errore('readin', 'Error generating temperatures: need temps(1) < temps(2)', 1)
    ELSE
      DO itemp = 1, nstemp
        temps(itemp) = tempsmin + DBLE(itemp - 1) * (tempsmax - tempsmin) / DBLE(nstemp - 1)
      ENDDO
    ENDIF
    WRITE(stdout, '(/,5x,a)') 'Generating evenly spaced temperature list.'
  ELSE IF (nstemp_hold .NE. nstemp) THEN !temps and nstemp not match
    ! Ignore nstemp setting, print warning
    WRITE(stdout, '(/,5x,a)') 'WARNING: Mismatch between temps(:) and nstemp'
    WRITE(stdout, '(/,5x,a)') 'WARNING: Using supplied temperature list and ignoring nstemp'
    nstemp = nstemp_hold
  ELSE
    CALL errore('readin', 'Error generating temperatures: unknown error', 1)
  END IF
  ! go from K to Ry
  temps(:) = temps(:) * kelvin2eV / ryd2ev
  !
  ! Check whether calculations are compatible with zero temperature
  !
  IF (MAXVAL(temps(:)) == 0.0) THEN
    IF (.NOT. (elecselfen .OR. phonselfen)) CALL errore('readin', &
        'temps = 0.0 can be used only for elecselfen and phonselfen', 1)
    IF (.NOT. efermi_read) CALL errore('readin', 'If temps = 0.0, fermi energy &
      must be explicitly given by efermi_read = .true. and fermi_energy = value', 1)
  ENDIF
  !
  CALL mp_bcast(nstemp, meta_ionode_id, world_comm)
  CALL mp_bcast(temps, meta_ionode_id, world_comm)
  ALLOCATE(gtemp(nstemp), STAT = ierr)
  IF (ierr /= 0) CALL errore('readin', 'Error allocating gtemp', 1)
  gtemp(:) = temps(1:nstemp)
  !
  ! In the case of Fermi-Dirac distribution one should probably etemp instead of degauss.
  ! This is achieved with assume_metal == .true.
  IF (ngaussw == -99 .AND. .NOT. assume_metal) THEN
    WRITE(stdout, '(/,5x,a)') 'WARNING: You are using ngaussw == -99 (Fermi-Dirac).'
    WRITE(stdout, '(/,5x,a)') '         You probably need assume_metal == .true '
  ENDIF
  !!!!!
  ! Some controls on impurity scattering input
  IF (assume_metal .AND. ii_g) THEN
    CALL errore('readin', 'Error: impurity matrix elements not compatable with metals', 1)
  ENDIF
  IF (ii_scattering .AND. .NOT. ii_g) THEN
    CALL errore('readin', 'Error: ii_g must = .true. if ii_scattering = .true.', 1)
  ENDIF
  IF (ii_partion .AND. ii_eda < eps20) THEN
    WRITE(stdout, '(/,5x,a)') '--------------------------------------------------------------------------------------'
    WRITE(stdout, '(/,5x,a)') 'WARNING: ii_partion == .TRUE. but dopant ionization energy ii_eda == 0.0 eV.'
    WRITE(stdout, '(/,5x,a)') 'Results for partial ionizaton may not be physical.'
    WRITE(stdout, '(/,5x,a)') 'if ii_partion == .true., please set ii_eda to a reasonable physical ionization energy in eV.'
    WRITE(stdout, '(/,5x,a)') '--------------------------------------------------------------------------------------'
  ENDIF
  IF (ii_partion .AND. ii_eda > one) THEN
    WRITE(stdout, '(/,5x,a)') '--------------------------------------------------------------------------------------'
    WRITE(stdout, '(/,5x,a)') 'WARNING: dopant ionization energy ii_eda > 1.0 eV.'
    WRITE(stdout, '(/,5x,a)') 'Please check if correct, results may not be physical.'
    WRITE(stdout, '(/,5x,a)') '--------------------------------------------------------------------------------------'
  ENDIF
  IF (ii_n < zero) CALL errore('readin', 'Error: Ionized impurity density ii_n must be > 0', 1)
  !
  IF (gb_scattering) THEN
    IF (gb_size < eps20) THEN
      CALL errore('readin', 'Error: Grain size gb_size must be > 0', 1)
    ELSE
      ! from nm to Bohr
      gb_size = gb_size / bohr2nm
    ENDIF
    IF (ii_only) CALL errore('readin', 'Error: Either gb_scattering or ii_only can be .TRUE.', 1 )
  ELSE
    IF (gb_only) CALL errore('readin', 'Error: gb_scattering must be .TRUE.', 1)
  ENDIF
  IF (ii_only .AND. gb_only) CALL errore('readin', 'Error: Either ii_only or gb_only can be .TRUE.', 1)
  IF (ii_scattering .AND. gb_only) CALL errore('readin', 'Error: Either ii_scattering or gb_only can be .TRUE.', 1)
  !
  ! thickness and smearing width of the Fermi surface
  ! from eV to Ryd
  fsthick     = fsthick / ryd2ev
  degaussw    = degaussw / ryd2ev
  ii_eda      = ii_eda / ryd2ev
  delta_smear = delta_smear / ryd2ev
  !
  ! smearing of phonon in a2f
  ! from meV to Ryd
  degaussq = degaussq / ryd2mev
  delta_qsmear = delta_qsmear / ryd2mev
  !
  ! Max frequency for the spectral decomposition of mobility from meV to Ry.
  mob_maxfreq = mob_maxfreq / ryd2mev
  !
  ! fermi_energy read from the input file from eV to Ryd
  IF (efermi_read) THEN
    fermi_energy = fermi_energy / ryd2ev
  ENDIF
  ! bfield: input in Tesla=[V/(m^2/s)] , convert in [eV/(Ang^2/s)]
  bfieldx = bfieldx * electron_si * ang2m**(-2)
  bfieldy = bfieldy * electron_si * ang2m**(-2)
  bfieldz = bfieldz * electron_si * ang2m**(-2)
  !
  ! eptemp : temperature for the electronic Fermi occupations in the e-p calculation (units of Kelvin)
  ! 1 K in eV = 8.6173423e-5
  ! Out-of bound issue with GCC compiler. Multiple Fermi temp is not used anyway.
  !
  ! from cm-1 to Ryd
  eps_acoustic = eps_acoustic / ev2cmm1 / ryd2ev
  !
  ! wmin and wmax from eV to Ryd
  wmin = wmin / ryd2ev
  wmax = wmax / ryd2ev
  !
  ! wmin_specfun and wmax_specfun from eV to Ryd
  wmin_specfun = wmin_specfun / ryd2ev
  wmax_specfun = wmax_specfun / ryd2ev
  !
  ! Scissor going from eV to Ryd
  scissor = scissor / ryd2ev
  !
  ! Photon energies for indirect absorption from eV to Ryd
  omegamin = omegamin / ryd2ev
  omegamax = omegamax / ryd2ev
  omegastep = omegastep / ryd2ev
  !
  xq(:) = zero
  !
  tmp_dir = trimcheck(outdir)
  dvscf_dir = trimcheck(dvscf_dir)
  !
  CALL bcast_input()
  !
  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !
  modenum_aux = modenum
  !
  ! SP: This initialized xk, nspin and nspin_mag
  IF (epwread .AND. .NOT. epbread) THEN
    CONTINUE
  ELSE
    ! In read_file, the call to allocate_wfc allocate evc with dimension ALLOCATE(evc(npwx*npol, nbnd))
    CALL read_file()
    !
    ! We define the global list of coarse grid k-points (cart and cryst)
    ALLOCATE(xk_all(3, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('readin', 'Error allocating xk_all', 1)
    ALLOCATE(et_all(nbnd, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('readin', 'Error allocating et_all', 1)
    ALLOCATE(isk_all(nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('readin', 'Error allocating isk_all', 1)
    ALLOCATE(xk_cryst(3, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('readin', 'Error allocating xk_cryst', 1)
    xk_all(:, :)   = zero
    et_all(:, :)   = zero
    isk_all(:)    = 0
    xk_cryst(:, :) = zero
    DO ik = 1, nkstot
      xk_all(:, ik)   = xk(:, ik)
      isk_all(ik)     = isk(ik)
      et_all(:, ik)   = et(:, ik)
      xk_cryst(:, ik) = xk(:, ik)
    ENDDO
    !  bring k-points from cartesian to crystal coordinates
    CALL cryst_to_cart(nkstot, xk_cryst, at, -1)
    ! Only master has the correct full list of kpt. Therefore bcast to all cores
    CALL mp_bcast(xk_all, meta_ionode_id, world_comm)
    CALL mp_bcast(et_all, meta_ionode_id, world_comm)
    CALL mp_bcast(isk_all, meta_ionode_id, world_comm)
    CALL mp_bcast(xk_cryst, meta_ionode_id, world_comm)
    !
    ! We define the local list of kpt
    ALLOCATE(xk_loc(3, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('readin', 'Error allocating xk_loc', 1)
    ALLOCATE(et_loc(nbnd, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('readin', 'Error allocating et_loc', 1)
    ALLOCATE(isk_loc(nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('readin', 'Error allocating isk_loc', 1)
    xk_loc(:, :) = zero
    et_loc(:, :) = zero
    isk_loc(:) = 0
    DO ik = 1, nks
      xk_loc(:, ik) = xk(:, ik)
      et_loc(:, ik) = et(:, ik)
      isk_loc(ik)   = isk(ik)
    ENDDO
    !
    ! 04-2019 - SP
    ! isk_loc and isk_all are spin index (LSDA only) on the local or all k-points.
    ! Those variable are introduced here for potential use but are not currently used further in EPW
    ! One would need to interpolate isk on the fine grids in ephwann_shuffle.
    !
  ENDIF
  !
  ! nbnd comes out of readfile
  IF (nbndsub == 0) nbndsub = nbnd
  !
#if defined(__MPI)
  IF (.NOT. (me_pool /= 0 .OR. my_pool_id /= 0)) THEN
    nk1 = nk1tmp
    nk2 = nk2tmp
    nk3 = nk3tmp
  ENDIF
#else
  nk1 = nk1tmp
  nk2 = nk2tmp
  nk3 = nk3tmp
#endif
  !
  IF (gamma_only) CALL errore('readin',&
     'cannot start from pw.x data file using Gamma-point tricks',1)
  !
  IF (modenum_aux /= -1) THEN
    modenum = modenum_aux
    iswitch = -4
  ELSEIF (modenum == 0) THEN
    iswitch = -2
  ELSE
    iswitch = -4
  ENDIF
  !
  CALL mp_bcast(iswitch, meta_ionode_id, world_comm)
  CALL mp_bcast(modenum, meta_ionode_id, world_comm)
  !
  IF (tfixed_occ) CALL errore('readin', 'phonon with arbitrary occupations not tested', 1)
  !
  IF (elph .AND. lsda) CALL errore('readin', 'El-ph and spin not implemented', 1)
  !
  IF (noncolin .AND. okpaw) THEN
    WRITE(stdout, '(a)') &
       'WARNING: readin: Some features are still experimental in ph.x with PAW and noncolin=.true.'
    WRITE(stdout, '(21x,a)') 'In this case, use the EPW results at your own risk.'
  ENDIF
  !
  !   There might be other variables in the input file which describe
  !   partial computation of the dynamical matrix. Read them here
  !
  CALL allocate_part(nat)
  IF (me_pool /= 0 .OR. my_pool_id /=0) GOTO 800
  IF (nat_todo < 0 .OR. nat_todo > nat) &
    CALL errore('readin', 'nat_todo is wrong', 1)
  IF (nat_todo /= 0) THEN
    IF (meta_ionode) READ(5, *, IOSTAT = ios) (atomo(na), na = 1, nat_todo)
    CALL mp_bcast(ios, meta_ionode_id, world_comm)
    CALL mp_bcast(atomo, meta_ionode_id, world_comm)
  ENDIF
800 CONTINUE
  !
#if defined(__MPI)
  CALL mp_bcast(nat_todo, meta_ionode_id, world_comm)
  IF (nat_todo > 0) THEN
    CALL mp_bcast(atomo, meta_ionode_id, world_comm)
  ENDIF
#endif
  !
  DO it = 1, ntyp
    IF (amass(it) <= 0.d0) CALL errore('readin', 'Wrong masses', it)
  ENDDO
  !
  ! broadcast the values of nq1, nq2, nq3
  !
  CALL mp_bcast(nq1, meta_ionode_id, world_comm)
  CALL mp_bcast(nq2, meta_ionode_id, world_comm)
  CALL mp_bcast(nq3, meta_ionode_id, world_comm)
  CALL mp_bcast(nk1, meta_ionode_id, world_comm)
  CALL mp_bcast(nk2, meta_ionode_id, world_comm)
  CALL mp_bcast(nk3, meta_ionode_id, world_comm)
  amass = AMU_RY * amass
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE readin
  !-----------------------------------------------------------------------

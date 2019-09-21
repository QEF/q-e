  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.ORg/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/phq_readin - Quantum-ESPRESSO group               
  !-----------------------------------------------------------------------
  SUBROUTINE epw_readin()
  !-----------------------------------------------------------------------
  !!
  !! This routine reads the control variables for the program epw.
  !! from standard input (unit 5).
  !! A second routine readfile reads the variables saved on a file
  !! by the self-consistent program.
  !!
  !! @Note:
  !!   SP: Image parallelization added
  !!
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
  USe disp,          ONLY : nq1, nq2, nq3
  USE epwcom,        ONLY : delta_smear, nsmear, dis_win_min, dis_win_max, wannierize, &
                            ngaussw, dvscf_dir, eptemp, bands_skipped, wdata, kmaps,   &
                            num_iter, dis_froz_max, fsthick, dis_froz_min, eig_read,   &
                            vme, degaussw, epexst, epwwrite, epbread, phonselfen, nqc2,&
                            elecselfen, a2f, plselfen, specfun_pl, nest_fn, filukk,    &
                            rand_nk, rand_k, rand_nq, rand_q, nkc1, nkc2, nkc3, nqc1,  &
                            nqc3, nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, eps_acustic, nw, &
                            wmax, wmin, mp_mesh_q, mp_mesh_k, filqf, filkf, nswi, nc,  &
                            delta_qsmear, degaussq, band_plot, ephwrite, nstemp,       &
                            broyden_beta, conv_thr_raxis, tempsmax, tempsmin, temps,   &
                            broyden_ndim, wscut, wsfc, nqstep, limag, lreal, muc,      &
                            gap_edge, conv_thr_iaxis, nqsmear, iprint, wepexst, nswfc, &
                            epwread, eliashberg, imag_read, kerread, kerwrite, lunif,  &
                            fermi_energy, efermi_read, max_memlt, fila2f, pwc, nswc,   &
                            ep_coupling, nw_specfun, wmax_specfun, wmin_specfun,       &
                            laniso, lpolar, lifc, asr_typ, lscreen, scr_typ, nbndsub,  &
                            fermi_diff, smear_rpa, cumulant, bnd_cum, proj, write_wfn, &
                            iswitch, ntempxx, liso, lacon, lpade, etf_mem, epbwrite,   &
                            nsiter, conv_thr_racon, specfun_el, specfun_ph, nbndskip,  &
                            system_2d, delta_approx, title, int_mob, scissor,          &
                            iterative_bte, scattering, selecqread, epmatkqread,        &  
                            ncarrier, carrier, scattering_serta, restart, restart_freq,&
                            scattering_0rta, longrange, shortrange, scatread, use_ws,  &
                            restart_filq, prtgkk, nel, meff, epsiHEG, lphase,          &
                            omegamin, omegamax, omegastep, n_r, lindabs, mob_maxiter,  & 
                            auto_projections, scdm_proj, scdm_entanglement, scdm_mu,   & 
                            scdm_sigma 
  USE klist_epw,     ONLY : xk_all, xk_loc, xk_cryst, isk_all, isk_loc, et_all, et_loc
  USE elph2,         ONLY : elph
  USE constants_epw, ONLY : ryd2mev, ryd2ev, ev2cmm1, kelvin2eV, zero, eps20, electron_SI, ang2m
  USE io_files,      ONLY : tmp_dir, prefix
  USE control_flags, ONLY : iverbosity, modenum, gamma_only
  USE ions_base,     ONLY : amass
  USE mp_world,      ONLY : world_comm, mpime
  USE partial,       ONLY : atomo, nat_todo
  USE constants,     ONLY : AMU_RY, eps16
  USE mp_global,     ONLY : my_pool_id, me_pool
  USE io_global,     ONLY : meta_ionode, meta_ionode_id, ionode_id
  USE io_var,        ONLY : iunkf, iunqf
  USE noncollin_module, ONLY : npol
  USE wvfct,         ONLY : npwx
#if defined(__NAG)
  USE F90_UNIX_ENV,  ONLY : iargc, getarg
#endif
  !
  IMPLICIT NONE
  !
  LOGICAL, EXTERNAL :: imatches
  !! Does the title match 
  CHARACTER(LEN = 256) :: outdir
  !! Output directory
#if ! defined(__NAG)
  INTEGER :: iargc
#endif
  INTEGER :: ios
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
  INTEGER :: nk1tmp
  !! temp vars for saving kgrid info
  INTEGER :: nk2tmp
  !! temp vars for saving kgrid info
  INTEGER :: nk3tmp  
  !! temp vars for saving kgrid info
  INTEGER :: ierr
  !! Error status
  ! 
  NAMELIST / inputepw / &
       amass, outdir, prefix, iverbosity, fildvscf,                            &
       elph, nq1, nq2, nq3, nk1, nk2, nk3, nbndskip,  nbndsub,                 &
       filukk, epbread, epbwrite, epwread, epwwrite, etf_mem, kmaps,           &
       eig_read, wepexst, epexst, vme,                                         &
       degaussw, fsthick, eptemp,  nsmear, delta_smear,                        &
       dvscf_dir, ngaussw, epmatkqread, selecqread,                            &
       wannierize, dis_win_max, dis_win_min, dis_froz_min, dis_froz_max,       &
       num_iter, proj, bands_skipped, wdata, iprint, write_wfn,                &
       wmin, wmax, nw, eps_acustic, a2f, nest_fn, plselfen,                    & 
       elecselfen, phonselfen, use_ws, nc,                                     &
       rand_q, rand_nq, rand_k, rand_nk, specfun_pl,                           &
       nqf1, nqf2, nqf3, nkf1, nkf2, nkf3,                                     &
       mp_mesh_k, mp_mesh_q, filqf, filkf, ephwrite,                           & 
       band_plot, degaussq, delta_qsmear, nqsmear, nqstep,                     &
       nswfc, nswc, nswi, pwc, wsfc, wscut, system_2d,                         &
       broyden_beta, broyden_ndim, nstemp, tempsmin, tempsmax, temps,          &
       conv_thr_raxis, conv_thr_iaxis, conv_thr_racon,                         &
       gap_edge, nsiter, muc, lreal, limag, lpade, lacon, liso, laniso, lpolar,& 
       lscreen, scr_typ, fermi_diff, smear_rpa, cumulant, bnd_cum,             &
       lifc, asr_typ, lunif, kerwrite, kerread, imag_read, eliashberg,         & 
       ep_coupling, fila2f, max_memlt, efermi_read, fermi_energy,              &
       specfun_el, specfun_ph, wmin_specfun, wmax_specfun, nw_specfun,         & 
       delta_approx, scattering, int_mob, scissor, ncarrier, carrier,          &
       iterative_bte, scattering_serta, scattering_0rta, longrange, shortrange,&
       scatread, restart, restart_freq, restart_filq, prtgkk, nel, meff,       &
       epsiHEG, lphase, omegamin, omegamax, omegastep, n_r, lindabs,           & 
       mob_maxiter, auto_projections, scdm_proj, scdm_entanglement, scdm_mu,   & 
       scdm_sigma
  ! tphases, fildvscf0
  !
  ! amass    : atomic masses
  ! iverbosity   : verbosity control
  ! outdir   : directory where input, output, temporary files reside
  ! elph     : if true calculate electron-phonon coefficients
  ! prefix   : the prefix of files produced by pwscf
  ! fildvscf : output file containing deltavsc
  ! fildrho  : output file containing deltarho
  !
  ! added by @ FG
  !
  ! ngaussw  : smearing type for FS average after wann interp
  ! degaussw : corresponding width (units of eV)
  ! filqf    : file with fine q kmesh for interpolation
  ! filkf    : file with fine kmesh for interpolation
  ! filukk   : file with rotation matrix U(k) for interpolation
  ! tphases  : if true set absolute unitary gauge for eigenvectors
  ! epstrict : if true use strict selection rule for phonon linewidht calculation
  ! fsthick  : the thickness of the Fermi shell for averaging the e-ph matrix elements (units of eV)
  ! eptemp   : temperature for the electronic Fermi occupations in the e-p calculation (units of Kelvin)
  ! fildvscf0: file containing deltavscf to be used as fake perturbation to set phases
  ! nw       : nr. of bins for frequency scan in \delta( e_k - e_k+q - w ) (units of eV)
  ! wmin     : min frequency for frequency scan in \delta( e_k - e_k+q - w ) (units of eV)
  ! wmax     : max    "  "  "                                    (units of eV)                   
  ! nbndsub  : number of bands in the optimal subspace (when disentanglement is used)
  ! elecselfen: if .TRUE. calculate imaginary part of electron selfenergy due to e-p interaction
  ! phonselfen: if .TRUE. calculate imaginary part of phonon selfenergy due to e-p interaction
  ! dvscf_dir: the dir containing all the .dvscf and .dyn files
  ! epbread  : read epmatq array from .epb files
  ! epbwrite : write epmatq array to .epb files
  ! nbndskip : number of bands to be skipped from the original Hamitonian (nfirstwin-1 in Marzari's notation)
  ! epwread  : read all quantities in Wannier representation from file epwdata.fmt
  ! epwwrite : write all quantities in Wannier representation to file epwdata.fmt
  !
  !  added by @jn
  !
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
  ! vme          : if .TRUE. calculate velocity matrix elements
  ! nest_fn      : if true, calculate the nesting function for a given set of q's
  ! nsmear       : number of smearing values to use for the selfen_phon call
  ! delta_smear  : change in energy for each additional nsmear ( units of eV)
  !
  ! added by @ RM
  !
  ! ephwrite    : if true write el-phonon matrix elements on the fine mesh to file
  ! eps_acustic : min phonon frequency for e-p and a2f calculations (units of cm-1)
  ! band_plot   : if true write files to plot band structure and phonon dispersion
  ! degaussq    : smearing for sum over q in e-ph coupling (units of meV)
  ! delta_qsmear: change in energy for each additional smearing in the a2f (units of meV)
  ! nqsmear     : number of smearings used to calculate a2f
  ! nqstep   : number of bins for frequency used to calculate a2f
  ! nswfc    : nr. of grid points between (0,wsfc) in Eliashberg equations 
  ! nswc     : nr. of grid points between (wsfc,wscut)
  ! pwc      : power used to define nswc for non-uniform grid real-axis calculations
  ! wsfc     : intermediate freqeuncy used for integration in Eliashberg equations (at least 2-3 times wsphmax) 
  ! wscut    : upper limit for frequency integration in Eliashberg equations (at least 5 times wsphmax) (units of eV)
  ! broyden_beta : mixing factor for broyden mixing
  ! broyden_ndim : number of iterations used in mixing scheme
  ! nstemp   : number of temperature points for which the Eliashberg equations are solved
  ! tempsmin : minimum temperature for which the Eliashberg equations are solved
  ! tempsmax : maximum temperature " "
  ! conv_thr_raxis : convergence threshold for iterative solution of real-axis Eliashberg equations
  ! conv_thr_iaxis : convergence threshold for iterative solution of imag-axis Eliashberg equations 
  ! conv_thr_racon : convergence threshold for iterative solution of analytic continuation of  
  !                  Eliashberg equations from imag- to real-axis
  ! gap_edge : initial guess of the superconducting gap (in eV)
  ! nsiter   : nr of iterations for self-consitency cycle
  ! muc     : effective Coulomb potential
  ! lreal   : if .TRUE. solve the real-axis Eliashberg eqautions
  ! limag   : if .TRUE. solve the imag-axis Eliashberg eqautions
  ! lpade   : if .TRUE. use pade approximants to continue imag-axis 
  !           Eliashberg equtions to real-axis
  ! lacon   : if .TRUE. use analytic continuation to continue imag-axis 
  !           Eliashberg equtions to real-axis
  ! liso    : if .TRUE. solve isotropic case
  ! laniso  : if .TRUE. solve anisotropic case
  ! lunif   : if .TRUE. a uniform grid is defined between wsfc and wscut for real-axis calculations
  ! kerwrite: if .TRUE. write Kp and Km to files .ker for real-axis calculations
  ! kerread : if .TRUE. read Kp and Km from files .ker for real-axis calculations
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
  ! system_2d    : if .TRUE. two-dimensional system (vaccum is in z-direction)  
  ! delta_approx : if .TRUE. the double delta approximation is used to compute the phonon self-energy 
  !
  ! added by CV & SP
  ! lpolar  : if .TRUE. enable the correct Wannier interpolation in the case of polar material.  
  ! lifc    : if .TRUE. reads interatomic force constants produced by q2r.x for phonon interpolation
  ! asr_typ : select type of ASR if lifc=.TRUE. (as in matdyn); otherwise it is the usual simple sum rule
  ! lscreen : if .TRUE. the e-ph matrix elements are screened by the RPA or TF dielectric function
  ! scr_typ : if 0 calculates the Lindhard screening, if 1 the Thomas-Fermi screening
  ! fermi_diff : difference between Fermi energy and band edge (in eV)
  ! smear_rpa  : smearing for the calculation of the Lindhard function (in eV)
  ! cumulant   : if .TRUE. calculates the electron spectral function using the cumulant expansion method
  !              (can be used as independent postprocessing by setting ep_coupling=.FALSE.)
  ! bnd_cum    : band index for which the cumulant calculation is done 
  !              (for more than one band, perform multiple calculations and add the results together)
  !
  ! 
  ! Added by SP
  !
  ! specfun_el      : if .TRUE. calculate electron spectral function due to e-p interaction
  ! specfun_ph      : if .TRUE. calculate phonon spectral function due to e-p interaction
  ! specfun_pl      : if .TRUE. calculate plason spectral function 
  ! restart         : if .TRUE. a run can be restarted from the interpolation level
  ! restart_freq    : Create a restart point every restart_freq q/k-points
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
  ! longrange       : if .TRUE. computes the long-range part of the el-ph (can
  !                   only be used with lpolar = .TRUE. )
  ! shortrange      : if .TRUE. computes the short-range part of the el-ph (can
  !                   only be used with lpolar = .TRUE. )
  ! prtgkk          : Print the vertex |g| [meV]. This generates huge outputs.   
  ! etf_mem         : if 0 no optimization, if 1 less memory is used for the fine grid interpolation
  !                   When etf_mem == 2, an additional loop is done on mode for the fine grid interpolation
  !                   part. This reduces the memory further by a factor "nmodes".    
  ! plselfen        : Calculate the electron-plasmon self-energy.
  ! nel             : Fractional number of electrons in the unit cell
  ! meff            : Density of state effective mass (in unit of the electron mass)
  ! epsiHEG         : Dielectric constant at zero doping
  ! lphase          : If .TRUE., fix the gauge on the phonon eigenvectors and electronic eigenvectors - DS 
  ! mob_maxiter     : Maximum number of iteration for the IBTE. 
  ! use_ws          : If .TRUE., use the Wannier-center to create the Wigner-Seitz cell. 
  ! epmatkqread     : If .TRUE., restart an IBTE calculation from scattering written to files. 
  ! selecqread      : If .TRUE., restart from the selecq.fmt file
  ! nc              : Number of carrier for the Ziman resistivity formula (can be fractional)
  !  
  ! Added by Manos Kioupakis
  ! omegamin  : Photon energy minimum
  ! omegamax  : Photon energy maximum
  ! omegastep : Photon energy step in evaluating phonon-assisted absorption spectra (in eV)
  ! n_r       :  constant refractive index
  ! lindabs   : do phonon-assisted absorption
  ! 
  nk1tmp = 0
  nk2tmp = 0
  nk3tmp = 0
  !
  IF (meta_ionode) THEN
    !
    ! ... Input from file ?
    CALL input_from_file( )
    !
    ! ... Read the first line of the input file
    !
    READ(5, '(A)', IOSTAT = ios) title
    !
  ENDIF
  ! 
  CALL mp_bcast(ios, meta_ionode_id, world_comm)
  CALL errore('epw_readin', 'reading title ', ABS(ios))
  CALL mp_bcast(title, meta_ionode_id, world_comm)
  !
  ! Rewind the input if the title is actually the beginning of inputph namelist
  !
  IF (imatches("&inputepw", title)) THEN
    WRITE(*, '(6x,a)') "Title line not specified: using 'default'."
    title = 'default'
    IF (meta_ionode) REWIND(5, IOSTAT = ios)
    CALL mp_bcast(ios, meta_ionode_id, world_comm  )
    CALL errore('epw_readin', 'Title line missing from input.', ABS(ios))
  ENDIF
  !
  IF (.NOT. meta_ionode) GOTO 400
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
  restart_freq = 100
  wannierize   = .FALSE.
  write_wfn    = .FALSE.
  kmaps        = .FALSE.
  nest_fn      = .FALSE.
  wepexst      = .FALSE.
  epexst       = .FALSE.
  eig_read     = .FALSE.
  dis_win_max  = 1d3
  dis_win_min  = -1d3
  dis_froz_max =  1d3
  dis_froz_min = -1d3
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
  wmin         = 0.d0
  wmax         = 0.3d0
  eps_acustic  = 5.d0 ! cm-1
  nw           = 10
  fsthick      = 1.d10 ! eV
  eptemp       = 300.0d0
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
  nbndskip     = 0
  nsmear       = 1
  delta_smear  = 0.01d0 ! eV
  modenum = 0 ! Was -1 previously and read from Modules/input_parameters.f90
              ! INTEGER :: modenum = 0. In QE 5, modenum variable does not exist
              ! anymore. Change the default EPW value to match the previous QE
              ! one. 
  vme = .TRUE. ! Was false by default until EPW 5.1
  ephwrite = .FALSE.
  band_plot = .FALSE.
  nqsmear = 10
  nqstep = 500
  delta_qsmear = 0.05d0 ! meV 
  degaussq = 0.05d0 ! meV
  lreal   = .FALSE.
  lunif   = .TRUE.
  limag   = .FALSE.
  lpade   = .FALSE.
  lacon   = .FALSE.
  liso    = .FALSE.
  laniso  = .FALSE.
  lpolar  = .FALSE.
  lifc    = .FALSE.
  asr_typ = 'simple'
  lscreen = .FALSE. 
  scr_typ = 0
  fermi_diff  = 1.d0
  smear_rpa   = 0.05d0
  cumulant    = .FALSE.
  bnd_cum     = 1
  kerwrite    = .FALSE.
  kerread     = .FALSE.
  imag_read   = .FALSE.
  eliashberg  = .FALSE.
  ep_coupling = .TRUE.
  nswfc    = 0
  nswc     = 0
  nswi     = 0 
  pwc      = 1.d0
  wscut    = 0.d0
  wsfc     = 0.5d0 * wscut 
  broyden_beta = 0.7d0
  broyden_ndim = 8
  conv_thr_raxis = 5.d-04
  conv_thr_iaxis = 1.d-05
  conv_thr_racon = 5.d-04
  gap_edge = 0.d0
  nstemp   = 1
  tempsmin = 0.d0
  tempsmax = 0.d0
  temps(:) = 0.d0
  nsiter   = 40
  muc     = 0.d0
  fila2f  = ' '
  max_memlt = 2.85d0
  efermi_read = .FALSE.
  fermi_energy = 0.d0
  wmin_specfun = 0.d0 ! eV
  wmax_specfun = 0.3d0 ! eV
  nw_specfun = 100
  system_2d = .FALSE.
  scattering = .FALSE.
  scattering_serta = .FALSE.
  scatread = .FALSE.
  scattering_0rta = .FALSE.
  int_mob    = .FALSE.
  iterative_bte = .FALSE.
  scissor    = 0.d0 ! eV
  carrier    = .FALSE.
  ncarrier   = 0.d0 ! cm^-3
  longrange  = .FALSE.
  shortrange = .FALSE.  
  prtgkk     = .FALSE.
  nel        = 0.0d0
  meff       = 1.d0
  epsiHEG    = 1.d0 
  lphase     = .FALSE. 
  omegamin   = 0.d0  ! eV
  omegamax   = 10.d0 ! eV
  omegastep  = 1.d0  ! eV
  n_r        = 1.d0
  lindabs    = .FALSE.
  mob_maxiter= 50
  use_ws     = .FALSE.
  epmatkqread = .FALSE.
  selecqread = .FALSE.
  nc         = 4.0d0
  !
  !     reading the namelist inputepw
  !
#if defined(__CRAYY) 
  !   The Cray does not accept "err" and "iostat" together with a namelist
  READ(5, inputepw)
  ios = 0
#else
  !
  IF (meta_ionode) READ(5, inputepw, ERR = 200, IOSTAT = ios)
#endif
200 CALL errore('epw_readin', 'reading input_epw namelist', ABS(ios))
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
  IF (nsmear < 1) CALL errore('epw_readin', 'Wrong number of nsmears', 1)
  IF (iverbosity < 0 .OR. iverbosity > 4) CALL errore('epw_readin', 'Wrong iverbosity', 1)
  IF (epbread .AND. epbwrite) CALL errore('epw_readin', 'epbread cannot be used with epbwrite', 1)
  IF (epbread .AND. epwread) CALL errore('epw_readin', 'epbread cannot be used with epwread', 1)
  IF (degaussw * 4.d0 > fsthick) CALL errore('epw_readin', ' degaussw too close to fsthick', 1)
  IF (nbndskip < 0) CALL errore('epw_readin', ' nbndskip must not be less than 0', 1)
  IF ((nw < 1) .OR. (nw > 1000)) CALL errore ('epw_readin', 'unreasonable nw', 1)
  IF (elecselfen .AND. plselfen) CALL errore('epw_readin', &
      'Electron-plasmon self-energy cannot be computed with electron-phonon', 1)
  IF (phonselfen .AND. plselfen) CALL errore('epw_readin', &
      'Electron-plasmon self-energy cannot be computed with electron-phonon', 1)
  IF (specfun_el .AND. plselfen) CALL errore('epw_readin', &
      'Electron-plasmon self-energy cannot be computed with el-ph spectral function', 1)
  IF (specfun_ph .AND. plselfen) CALL errore('epw_readin', &
      'Electron-plasmon self-energy cannot be computed with el-ph spectral function', 1)
  IF (elecselfen .AND. specfun_pl ) CALL errore('epw_readin', &
      'Electron-plasmon spectral function cannot be computed with electron-phonon', 1)
  IF (phonselfen .AND. specfun_pl) CALL errore('epw_readin', &
      'Electron-plasmon spectral function cannot be computed with electron-phonon', 1)
  IF (specfun_el .AND. specfun_pl) CALL errore('epw_readin', &
      'Electron-plasmon spectral function cannot be computed with el-ph spectral function', 1)
  IF (specfun_ph .AND. specfun_pl) CALL errore('epw_readin', &
      'Electron-plasmon spectral function cannot be computed with el-ph spectral function', 1)
  IF (a2f .AND. .NOT. phonselfen) CALL errore('epw_readin', 'a2f requires phonoselfen', 1)
  IF (elph .AND. .NOT. ep_coupling) CALL errore('epw_readin', 'elph requires ep_coupling=.true.', 1)
  IF ((elph .AND. wannierize) .AND. (epwread)) CALL errore('epw_readin', &
      'must use same w90 rotation matrix for entire run', 1)
  IF (wannierize .AND. .NOT. ep_coupling) CALL errore('epw_readin', &
      'wannierize requires ep_coupling=.true.', 1)
  IF ((wmin > wmax)) CALL errore('epw_readin', ' check wmin, wmax ', 1)
  IF ((wmin_specfun > wmax_specfun)) CALL errore('epw_readin', 'check wmin_specfun, wmax_specfun', 1)
  IF ((nw_specfun < 2)) CALL errore('epw_readin', 'nw_specfun must be at least 2', 1)
  IF ((nqstep < 2)) CALL errore('epw_readin', 'nqstep must be at least 2', 1)
  IF ((nbndsub > 200)) CALL errore('epw_readin', 'too many wannier functions increase size of projx', 1)
  IF ((phonselfen .OR. elecselfen .OR. specfun_el .OR. specfun_ph) .AND. (mp_mesh_k .OR. mp_mesh_q)) & 
    CALL errore('epw_readin', 'can only work with full uniform mesh', 1)
  IF (ephwrite .AND. .NOT. ep_coupling .AND. .NOT. elph) CALL errore('epw_readin', &
      'ephwrite requires ep_coupling=.TRUE., elph=.TRUE.', 1)
  IF (ephwrite .AND. (rand_k .OR. rand_q)) &
    CALL errore('epw_readin', 'ephwrite requires a uniform grid', 1) 
  IF (ephwrite .AND. (MOD(nkf1, nqf1) /= 0 .OR. MOD(nkf2, nqf2) /= 0 .OR. MOD(nkf3, nqf3) /= 0)) &
    CALL errore('epw_readin', 'ephwrite requires nkf1,nkf2,nkf3 to be multiple of nqf1,nqf2,nqf3', 1)
  IF (band_plot .AND. filkf == ' ' .AND. filqf == ' ') CALL errore('epw_readin', &
      'plot band structure and phonon dispersion requires k- and q-points read from filkf and filqf files', 1)
  IF (band_plot .AND. filkf /= ' ' .AND. (nkf1 > 0 .OR. nkf2 > 0 .OR. nkf3 > 0)) CALL errore('epw_readin', &
      'You should define either filkf or nkf when band_plot = .true.', 1)
  IF (band_plot .AND. filqf /= ' ' .AND. (nqf1 > 0 .OR. nqf2 > 0 .OR. nqf3 > 0)) CALL errore('epw_readin', &
     'You should define either filqf or nqf when band_plot = .true.', 1)
  IF (filkf /= ' ' .AND. .NOT. efermi_read) CALL errore('epw_readin', &
     'WARNING: if k-points are along a line, then efermi_read=.true. and fermi_energy must be given in the input file', -1)
  IF (scattering .AND. nstemp < 1) CALL errore('epw_readin', 'wrong number of nstemp', 1)
  IF (scattering .AND. MAXVAL(temps(:)) > 0.d0 .AND. tempsmin > 0.d0 .AND. tempsmax > 0.d0) &
    CALL errore('epw_readin', 'define either (tempsmin and tempsmax) or temps(:)', 1)
  IF (scattering .AND. tempsmax < tempsmin) CALL errore('epw_readin', 'tempsmax should be greater than tempsmin', 1)
  IF ((ABS(ncarrier) > 1E+5) .AND. .NOT. carrier) CALL errore('epw_readin', &
      'carrier must be .TRUE. if you specify ncarrier.', 1)
  IF (carrier .AND. (ABS(ncarrier) < 1E+5))  CALL errore('epw_readin', &
      'The absolute value of the doping carrier concentration must be larger than 1E5 cm^-3', 1)
  IF ((longrange .OR. shortrange) .AND. (.NOT. lpolar)) CALL errore('epw_readin', &
      'Error: longrange or shortrange can only be true if lpolar is true as well.', 1)
  IF (longrange .AND. shortrange) CALL errore('epw_readin',&
      'Error: longrange and shortrange cannot be both true.', 1)
  IF (epwread .AND. .NOT. kmaps .AND. .NOT. epbread) CALL errore('epw_readin', &
      'Error: kmaps has to be true for a restart run. ', 1)
  IF (.NOT. epwread .AND. .NOT. epwwrite) CALL errore('epw_readin', &
      'Error: Either epwread or epwwrite needs to be true. ', 1)
  IF (lscreen .AND. etf_mem == 2) CALL errore('epw_readin', 'Error: lscreen not implemented with etf_mem=2', 1)
  IF (ABS(degaussw) < eps16 .AND. etf_mem == 2) CALL errore('epw_readin', &
      'Error: adapt_smearing not implemented with etf_mem=2', 1)
  ! Make sure the files exists
  IF (filkf /= ' ') THEN
    OPEN(UNIT = iunkf, FILE = filkf, STATUS = 'old', FORM = 'formatted', ERR = 100, IOSTAT = ios)
100 CALL errore('epw_readin', 'opening file ' // filkf, ABS(ios))
    CLOSE(iunkf)
  ENDIF
  IF (filqf /= ' ') THEN
    OPEN(UNIT = iunqf, FILE = filqf, STATUS = 'old', FORM = 'formatted', ERR = 101, IOSTAT = ios)
101 CALL errore('epw_readin', 'opening file ' // filqf, ABS(ios))
    CLOSE(iunqf)
  ENDIF  
  IF (iterative_bte) THEN
    ! The fine grids have to be homogeneous and the same. Otherwise the populations can oscillate. 
    IF (nkf1 /= nqf1 .OR. nkf2 /= nqf2 .OR. nkf3 /= nqf3) THEN
      CALL errore('epw_readin', 'Error: the fine k-points and q-points grids have to be the same when doing IBTE.', 1) 
    ENDIF
  ENDIF
  IF (auto_projections .AND. proj(1) /= ' ') CALL errore('epw_readin', & 
      'Cannot specify both auto_projections and projections block', 1)
  IF ((auto_projections .AND. .NOT. scdm_proj) .OR. (.NOT. auto_projections .AND. scdm_proj)) & 
    CALL errore('epw_readin', 'auto_projections require both scdm_proj=.true. and auto_projections=.true.', 1)
  !
  ! thickness and smearing width of the Fermi surface  
  ! from eV to Ryd
  fsthick     = fsthick / ryd2ev 
  degaussw    = degaussw / ryd2ev 
  delta_smear = delta_smear / ryd2ev 
  !
  ! smearing of phonon in a2f
  ! from meV to Ryd 
  degaussq = degaussq / ryd2mev
  delta_qsmear = delta_qsmear / ryd2mev
  !
  ! fermi_energy read from the input file from eV to Ryd
  IF (efermi_read) THEN
    fermi_energy = fermi_energy / ryd2ev
  ENDIF
  ! eptemp : temperature for the electronic Fermi occupations in the e-p calculation (units of Kelvin)
  ! 1 K in eV = 8.6173423e-5
  ! from K to Ryd
  ! Out-of bound issue with GCC compiler. Multiple Fermi temp is not used anyway.
  eptemp = eptemp * kelvin2eV / ryd2ev
  ! 
  ! from cm-1 to Ryd
  eps_acustic = eps_acustic / ev2cmm1 / ryd2ev 
  !
  !    reads the q point (just if ldisp = .FALSE.)
  !
  ! wmin and wmax from eV to Ryd
  wmin = wmin / ryd2ev
  wmax = wmax / ryd2ev
  !
  ! wmin_specfun and wmax_specfun from eV to Ryd
  wmin_specfun = wmin_specfun / ryd2ev
  wmax_specfun = wmax_specfun / ryd2ev
  !
  ! scissor going from eV to Ryd
  scissor = scissor / ryd2ev
  ! 
  ! Photon energies for indirect absorption from eV to Ryd
  omegamin = omegamin / ryd2ev
  omegamax = omegamax / ryd2ev
  omegastep = omegastep / ryd2ev
  IF (scattering) THEN
    DO i = 1, ntempxx
      IF (temps(i) > 0.d0) THEN
        nstemp = i
      ENDIF
    ENDDO
    !
    ! go from K to Ry
    temps(:) = temps(:) * kelvin2eV / ryd2ev
    tempsmin = tempsmin * kelvin2eV / ryd2ev
    tempsmax = tempsmax * kelvin2eV / ryd2ev
    !
  ENDIF
  ! 
  xq(:) = zero
  !
  tmp_dir = TRIM(outdir)
  dvscf_dir = TRIM(dvscf_dir) // '/'
  !
400 CONTINUE
  CALL bcast_epw_input()
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
    IF (ierr /= 0) CALL errore('epw_readin', 'Error allocating xk_all', 1)
    ALLOCATE(et_all(nbnd, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_readin', 'Error allocating et_all', 1)
    ALLOCATE(isk_all(nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_readin', 'Error allocating isk_all', 1)
    ALLOCATE(xk_cryst(3, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_readin', 'Error allocating xk_cryst', 1)
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
    CALL mp_bcast(xk_all, ionode_id, world_comm)
    CALL mp_bcast(et_all, ionode_id, world_comm)
    CALL mp_bcast(isk_all, ionode_id, world_comm)
    CALL mp_bcast(xk_cryst, ionode_id, world_comm)
    ! 
    ! We define the local list of kpt
    ALLOCATE(xk_loc(3, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_readin', 'Error allocating xk_loc', 1)
    ALLOCATE(et_loc(nbnd, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_readin', 'Error allocating et_loc', 1)
    ALLOCATE(isk_loc(nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_readin', 'Error allocating isk_loc', 1)
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
  IF (.NOT. (me_pool /=0 .OR. my_pool_id /=0)) THEN
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
  IF (gamma_only) CALL errore('epw_readin',&
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
  IF (tfixed_occ) CALL errore('epw_readin', 'phonon with arbitrary occupations not tested', 1)
  !
  IF (elph .AND. lsda) CALL errore('epw_readin', 'El-ph and spin not implemented', 1)
  !
  !   There might be other variables in the input file which describe
  !   partial computation of the dynamical matrix. Read them here
  !
  CALL allocate_part(nat)
  IF (me_pool /= 0 .OR. my_pool_id /=0) GOTO 800 
  IF (nat_todo < 0 .OR. nat_todo > nat) & 
    CALL errore('epw_readin', 'nat_todo is wrong', 1)
  IF (nat_todo /= 0) THEN
    IF (meta_ionode) READ(5, *, IOSTAT = ios) (atomo(na), na = 1, nat_todo)
    CALL mp_bcast(ios, meta_ionode_id, world_comm)
    CALL mp_bcast(atomo, meta_ionode_id, world_comm)
  ENDIF
800 CONTINUE
  CALL bcast_epw_input1()
  !
  DO it = 1, ntyp
    IF (amass(it) <= 0.d0) CALL errore('epw_readin', 'Wrong masses', it)
  ENDDO
  ! 
  !  broadcast the values of nq1, nq2, nq3
  !
  CALL mp_bcast(nq1, meta_ionode_id, world_comm)
  CALL mp_bcast(nq2, meta_ionode_id, world_comm)
  CALL mp_bcast(nq3, meta_ionode_id, world_comm)
  CALL mp_bcast(nk1, meta_ionode_id, world_comm)
  CALL mp_bcast(nk2, meta_ionode_id, world_comm)
  CALL mp_bcast(nk3, meta_ionode_id, world_comm)
  !
  amass = AMU_RY * amass
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE epw_readin
  !-----------------------------------------------------------------------

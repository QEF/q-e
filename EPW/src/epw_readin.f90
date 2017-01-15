  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/phq_readin - Quantum-ESPRESSO group               
  !-----------------------------------------------------------------------
  SUBROUTINE epw_readin
  !-----------------------------------------------------------------------
  !!
  !!    This routine reads the control variables for the program epw.
  !!    from standard input (unit 5).
  !!    A second routine readfile reads the variables saved on a file
  !!    by the self-consistent program.
  !!
  !!   @Note:
  !!     SP: Image parallelization added
  !!
  USE ions_base,     ONLY : nat, ntyp => nsp
  USE mp,            ONLY : mp_bcast 
  USE pwcom,         ONLY : xqq
  USE wvfct,         ONLY : nbnd
  USE klist,         ONLY : nks 
  USE lsda_mod,      ONLY : lsda
  USE fixed_occ,     ONLY : tfixed_occ
  USE check_stop,    ONLY : time_max => max_seconds
  USE qpoint,        ONLY : xq
  USE disp,          ONLY : nq1, nq2, nq3
  USE output,        ONLY : fildvscf, fildrho
  USE epwcom,        ONLY : delta_smear, nsmear, dis_win_min, dis_win_max, wannierize, &
                            ngaussw, dvscf_dir, eptemp, wdata, &
                            num_iter, dis_froz_max, fsthick, dis_froz_min, &
                            vme, degaussw, epexst, eig_read, kmaps, &
                            epwwrite, epbread, phonselfen, elecselfen, &
                            a2f, rand_k, rand_nq, rand_q, &
                            parallel_q, parallel_k, nkf1, &
                            nkf2, nkf3, nqf1, nqf2, nqf3, rand_nk, &
                            nest_fn, eps_acustic, nw, wmax, wmin, &
                            mp_mesh_q, filqf, filkf, delta_qsmear, degaussq, &
                            band_plot, ephwrite, mp_mesh_k, nstemp, broyden_beta, &
                            conv_thr_raxis, tempsmax, tempsmin, temps, broyden_ndim, &
                            wscut, wsfc, nqstep, limag, lreal, muc, gap_edge, &
                            conv_thr_iaxis, nqsmear, iprint, wepexst, epwread, & 
                            eliashberg, imag_read, kerread, kerwrite, lunif, specfun, &
                            fermi_energy, efermi_read, max_memlt, fila2f, &
                            ep_coupling, nw_specfun, wmax_specfun, &
                            wmin_specfun, laniso, lpolar, lifc, asr_typ, &
                            proj, write_wfn, iswitch, ntempxx, &
                            liso, lacon, lpade, etf_mem, epbwrite, &
                            nsiter, conv_thr_racon, &
                            pwc, nswc, nswfc, nswi, filukq, filukk, &
                            nbndsub, nbndskip, system_2d, delta_approx, &
                            title, int_mob, scissor, iterative_bte, scattering, &
                            ncarrier, carrier, scattering_serta, &
                            scattering_0rta, longrange, shortrange,restart, &
                            restart_freq
  USE elph2,         ONLY : elph
  USE start_k,       ONLY : nk1, nk2, nk3
  USE constants_epw, ONLY : ryd2mev, ryd2ev, ev2cmm1, kelvin2eV
  USE io_files,      ONLY : tmp_dir, prefix
  USE control_flags, ONLY : iverbosity, modenum, gamma_only
  USE ions_base,     ONLY : amass
  USE mp_world,      ONLY : world_comm
  USE partial,       ONLY : atomo, nat_todo
  USE constants,     ONLY : AMU_RY
  USE control_lr,    ONLY : lgamma
  USE mp_global,     ONLY : my_pool_id, me_pool
  USE io_global,     ONLY : meta_ionode, meta_ionode_id, ionode, ionode_id, stdout
#if defined(__NAG)
  USE F90_UNIX_ENV,  ONLY : iargc, getarg
#endif
  !
  implicit none
  !
#if ! defined(__NAG)
  integer :: iargc
#endif
  INTEGER :: ios
  !! integer variable for I/O control
  INTEGER :: na
  !! counter on polarizations
  INTEGER :: it
  !! counter on iterations
  INTEGER :: modenum_aux
  !! auxilary variable for saving the modenum
  INTEGER :: i
  !! Counter for loops
  INTEGER :: nk1tmp
  !! temp vars for saving kgrid info
  INTEGER :: nk2tmp
  !! temp vars for saving kgrid info
  INTEGER :: nk3tmp  
  !! temp vars for saving kgrid info
  LOGICAL, EXTERNAL  :: imatches
  character(len=256) :: outdir
  namelist / inputepw / &
       amass, outdir, prefix, iverbosity, time_max, fildvscf,                  &
       elph, nq1, nq2, nq3, nk1, nk2, nk3, nbndskip,  nbndsub,                 &
       filukk, filukq, epbread, epbwrite, epwread, epwwrite, etf_mem, kmaps,   &
       eig_read, wepexst, epexst, vme,                                         &
       degaussw, fsthick, eptemp,  nsmear, delta_smear,                        &
       dvscf_dir, ngaussw,                                                     &
       wannierize, dis_win_max, dis_win_min, dis_froz_min, dis_froz_max,       &
       num_iter, proj,  wdata, iprint, write_wfn, wmin, wmax, nw,              &
       eps_acustic, a2f, nest_fn,                                              & 
       elecselfen, phonselfen, parallel_k, parallel_q,                         &
       rand_q, rand_nq, rand_k, rand_nk,                                       &
       nqf1, nqf2, nqf3, nkf1, nkf2, nkf3,                                     &
       mp_mesh_k, mp_mesh_q, filqf, filkf, ephwrite,                           & 
       band_plot, degaussq, delta_qsmear, nqsmear, nqstep,                     &
       nswfc, nswc, nswi, pwc, wsfc, wscut,                                    &
       broyden_beta, broyden_ndim, nstemp, tempsmin, tempsmax, temps,          &
       conv_thr_raxis, conv_thr_iaxis, conv_thr_racon,                         &
       gap_edge, nsiter, muc, lreal, limag, lpade, lacon, liso, laniso, lpolar,& 
       lifc, asr_typ, lunif, kerwrite, kerread, imag_read, eliashberg,         & 
       ep_coupling, fila2f, max_memlt, efermi_read, fermi_energy,              &
       specfun, wmin_specfun, wmax_specfun, nw_specfun, system_2d,             & 
       delta_approx, scattering, int_mob, scissor, ncarrier, carrier,          &
       iterative_bte, scattering_serta, scattering_0rta, longrange, shortrange,&
       restart, restart_freq

  ! tphases, fildvscf0
  !
  ! amass    : atomic masses
  ! iverbosity   : verbosity control
  ! outdir   : directory where input, output, temporary files reside
  ! elph     : if true calculate electron-phonon coefficients
  ! time_max : maximum cputime for this run
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
  ! filukq   : file with rotation matrix U(k+q) for interpolation
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
  ! wannierize : if .TRUE. run the wannier90 code to maximally localize the WFs
  ! dis_win_min : lower bound on wannier90 disentanglement window
  ! dis_win_max : upper bound on wannier90 disentanglement window
  ! dis_froz_min : lower bound on frozen wannier90 disentanglement window
  ! dis_froz_max : upper bound on frozen wannier90 disentanglement window
  ! num_iter     : number of iterations used in the wannier90 minimisation
  ! proj         : initial projections (states) of the wannier functions before minimization
  ! wdata        : Empty array that can be used to pass extra info to prefix.win file, for things not explicitly declared here 
  ! iprint       : verbosity of the wannier90 code
  ! write_wfn    : writes out UNK files from pwscf run for plotting of XSF files
  ! etf_mem      : if .true., the fine Bloch-space e-ph matrix elements are stored in memory
  ! kmaps        : if true, read kmap and kgmap from disk (prior run)
  ! eig_read     : if .true. then readin a set of electronic eigenvalues in eV to replace the calcualted ones
  ! wepexst      : if .TRUE. prefix.epmatwe files are already on disk. don't recalculate. debugging param
  ! epexst       : if .TRUE. prefix.epmatwp files are already on disk. don't recalculate  debugging param
  ! vme          : if .TRUE. calculate velocity matrix elements
  ! nest_fn      : if true, calculate the nesting function for a given set of q's
  ! nsmear       : number of smearing values to use for the selfen_phon call
  ! delta_smear  : change in energy for each additional nsmear ( units of eV)
  !
  ! added by @ RM
  !
  ! ephwrite : if true write el-phonon matrix elements on the fine mesh to file
  ! eps_acustic : min phonon frequency for e-p and a2f calculations (units of cm-1)
  ! band_plot : if true write files to plot band structure and phonon dispersion
  ! degaussq : smearing for sum over q in e-ph coupling (units of meV)
  ! delta_qsmear : change in energy for each additional smearing in the a2f (units of meV)
  ! nqsmear : number of smearings used to calculate a2f
  ! nqstep : number of bins for frequency used to calculate a2f
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
  ! lreal   : if .true. solve the real-axis Eliashberg eqautions
  ! limag   : if .true. solve the imag-axis Eliashberg eqautions
  ! lpade   : if .true. use pade approximants to continue imag-axis 
  !           Eliashberg equtions to real-axis
  ! lacon   : if .true. use analytic continuation to continue imag-axis 
  !           Eliashberg equtions to real-axis
  ! liso    : if .true. solve isotropic case
  ! laniso  : if .true. solve anisotropic case
  ! lunif   : if .true. a uniform grid is defined between wsfc and wscut for real-axis calculations
  ! kerwrite: if .true. write Kp and Km to files .ker for real-axis calculations
  ! kerread : if .true. read Kp and Km from files .ker for real-axis calculations
  ! imag_read : if .true. read from files Delta and Znorm on the imaginary-axis
  ! eliashberg : if .true. solve the Eliashberg equations
  ! ep_coupling : if .true. run e-p coupling calculation
  ! fila2f  : input file with eliashberg spectral function 
  ! max_memlt : maximum memory that can be allocated per pool 
  ! efermi_read : if. true. read from input file
  ! fermi_energy : fermi eneergy read from input file (units of eV)
  ! specfun  : if .TRUE. calculate electron spectral function due to e-p interaction
  ! wmin_specfun : min frequency in electron spectral function due to e-p interaction (units of eV)
  ! wmax_specfun : max frequency in electron spectral function due to e-p interaction (units of eV)
  ! nw_specfun   : nr. of bins for frequency in electron spectral function due to e-p interaction 
  ! system_2d : if .true. two-dimensional system (vaccum is in z-direction)  
  ! delta_approx: if .true. the double delta approximation is used to compute the phonon self-energy 
  !
  ! added by CV & SP
  ! lpolar : if .true. enable the correct Wannier interpolation in the case of polar material.  
  ! lifc : if .true. reads interatomic force constants produced by q2r.x for phonon interpolation
  ! asr_typ : select type of ASR if lifc=.true. (as in matdyn); otherwise it is the usual simple sum rule
  ! 
  ! Added by SP
  !
  ! restart         : if .true. a run can be restarted from the interpolation level
  ! restart_freq    : Create a restart point every restart_freq q/k-points
  ! scattering      : if .true. scattering rates are calculated
  ! scattering_serta: if .true. scattering rates are calculated using self-energy relaxation-time-approx
  ! scattering_0rta : if .true. scattering rates are calculated using 0th order relaxation-time-approx
  ! int_mob         : if .true. computes the intrinsic mobilities. This means that the
  !                   electron and hole carrier density is equal. 
  ! iterative_bte   : if .true. computes the iterative solution to the BTE. Need a
  !                   prior run with ERTA. 
  ! scissor         : Value of the scissor shitf in eV. This only affects the CBM of etf. Do you use in
  !                   metals obviously.
  ! carrier         : if .true. computes the doped carrier mobilities. 
  ! ncarrier        : Set the Fermi level so that the carrier concentration is
  !            "      ncarrier". If ncarrier > 0, electron doping, hole doping otherwise
  ! longrange       : if .true. computes the long-range part of the el-ph (can
  !                   only be used with lpolar = .true. )
  ! shortrange      : if .true. computes the short-range part of the el-ph (can
  !                   only be used with lpolar = .true. )
  !  
  CHARACTER (LEN=80)  :: input_file
  INTEGER             :: nargs, iiarg, ierr
  !
  nk1tmp = 0
  nk2tmp = 0
  nk3tmp = 0
  !
  IF (meta_ionode) THEN
  !
  ! ... Input from file ?
     CALL input_from_file ( )
  !
  ! ... Read the first line of the input file
  !
     READ( 5, '(A)', IOSTAT = ios ) title
  !
  ENDIF
  ! 
  CALL mp_bcast(ios, meta_ionode_id, world_comm )
  CALL errore( 'epw_readin', 'reading title ', ABS( ios ) )
  CALL mp_bcast(title, meta_ionode_id, world_comm  )
  !
  ! Rewind the input if the title is actually the beginning of inputph namelist
  !
  IF( imatches("&inputepw", title) ) THEN
    WRITE(*, '(6x,a)') "Title line not specified: using 'default'."
    title='default'
    IF (meta_ionode) REWIND(5, iostat=ios)
    CALL mp_bcast(ios, meta_ionode_id, world_comm  )
    CALL errore('epw_readin', 'Title line missing from input.', abs(ios))
  ENDIF
  !
  IF (.NOT. meta_ionode) goto 400
  !
  !   set default values for variables in namelist
  !
  amass(:)     = 0.d0
  iverbosity   = 0
  elph         = .false.
  elecselfen   = .false.
  phonselfen   = .false.
  specfun      = .false.
  epbread      = .false.
  epbwrite     = .false.
  epwread      = .false.
  epwwrite     = .false.
  restart      = .false.
  restart_freq = 100
  wannierize   = .false.
  write_wfn    = .false.
  kmaps        = .false.
  nest_fn      = .false.
  wepexst      = .false.
  epexst       = .false.
  eig_read     = .false.
  dis_win_max  = 1d3
  dis_win_min  = -1d3
  dis_froz_max = -0.9d3
  dis_froz_min = -1d3
  num_iter     = 200
  proj(:)      = ''
  wdata(:)     = ''
  iprint       = 2
  wmin         = 0.d0
  wmax         = 0.3d0
  eps_acustic  = 5.d0 ! cm-1
  nw           = 10
  fsthick      = 1.d10 ! eV
  eptemp       = 300.0d0
  degaussw     = 0.025d0 ! eV
!  tphases      = .false.
  parallel_k   = .true.
  parallel_q   = .false.
  a2f          = .false.
  etf_mem      = .true.
!  fildvscf0    = ' '
  ngaussw      = 1
  time_max     = 10000000.d0
  outdir       = '.'
  dvscf_dir    = '.'
  prefix       = 'pwscf'
  filqf        = ' '
  filkf        = ' '
  fildrho      = ' '
  fildvscf     = ' '
  filukk       = ' '
  rand_q       = .false.
  delta_approx = .false.
  rand_nq      = 1
  rand_k       = .false.
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
  mp_mesh_k    = .false.
  mp_mesh_q    = .false.
  nbndsub      = 0
  nbndskip     = 0
  nsmear       = 1
  delta_smear  = 0.01d0 ! eV
  modenum = 0 ! Was -1 previously and read from Modules/input_parameters.f90
              ! INTEGER :: modenum = 0. In QE 5, modenum variable does not exist
              ! anymore. Change the default EPW value to match the previous QE
              ! one. 
  vme = .false.
  ephwrite = .false.
  band_plot = .false.
  nqsmear = 10
  nqstep = 500
  delta_qsmear = 0.05d0 ! meV 
  degaussq = 0.05d0 ! meV
  lreal   = .false.
  lunif   = .true.
  limag   = .false.
  lpade   = .false.
  lacon   = .false.
  liso    = .false.
  laniso  = .false.
  lpolar  = .false.
  lifc    = .false.
  asr_typ = 'simple'
  kerwrite= .false.
  kerread = .false.
  imag_read = .false.
  eliashberg = .false.
  ep_coupling = .true.
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
  efermi_read = .false.
  fermi_energy = 0.d0
  wmin_specfun = 0.d0 ! eV
  wmax_specfun = 0.3d0 ! eV
  nw_specfun = 100
  system_2d = .false.
  scattering = .false.
  scattering_serta = .false.
  scattering_0rta = .false.
  int_mob = .false.
  iterative_bte = .false.
  scissor = 0.d0 ! eV
  carrier = .false.
  ncarrier = 0.d0 ! cm^-3
  longrange = .false.
  shortrange = .false.  
  !
  !     reading the namelist inputepw
  !
#if defined(__CRAYY) 
  !   The Cray does not accept "err" and "iostat" together with a namelist
  READ (5, inputepw)
  ios = 0
#else
  !
  IF (meta_ionode) READ (5, inputepw, err = 200, iostat = ios)
#endif
200 CALL errore ('epw_readin', 'reading input_epw namelist', abs (ios) )
  !
  nk1tmp = nk1
  nk2tmp = nk2
  nk3tmp = nk3
  !
  !     Check all namelist variables
  !
  IF (filukk.eq.' ') filukk=trim(prefix)//'.ukk'
  IF (nsmear .lt. 1) CALL errore ('epw_readin', &
       & 'Wrong number of nsmears',1)
  IF (iverbosity.ne.0.and.iverbosity.ne.1.and.iverbosity.ne.2.and.iverbosity.ne.3) CALL errore ('epw_readin', &
       &' Wrong  iverbosity ', 1)
  IF (time_max.lt.1.d0) CALL errore ('epw_readin', ' Wrong time_max', 1)
!  IF (tphases.and.fildvscf0.eq.' ') CALL errore ('epw_readin', &
!       &' tphases requires fildvscf0', 1)
  IF (epbread.and.epbwrite) CALL errore ('epw_readin', &
       &' epbread cannot be used with epbwrite', 1)
  IF (degaussw*4.d0.gt.fsthick) CALL errore ('epw_readin', &
       &' degaussw too close to fsthick', 1)
  IF ( nbndskip .lt. 0) CALL errore('epw_readin', &
       &' nbndskip must not be less than 0', 1)
  IF ((nw.lt.1).or.(nw.gt.1000)) CALL errore ('epw_readin', &
       &' unreasonable nw', 1)
  IF (parallel_k .and. parallel_q) CALL errore('epw_readin', &
       &'can only parallelize over k OR q',1)
  IF (.not.(parallel_k .or. parallel_q)) CALL errore('epw_readin', &
       &'must parallelize over k OR q',1)
  IF (parallel_k .and. elecselfen) CALL errore('epw_readin', &
       &'Electron selfenergy is more efficient with k_para',-1)
  IF (a2f .and. .not.phonselfen) CALL errore('epw_readin', &
       &'a2f requires phonoselfen',1)
  IF (parallel_q .and. phonselfen) CALL errore('epw_readin', &
       &'Phonon selfenergy is more efficient with q_para',-1)
  IF (parallel_q) CALL errore('epw_readin', &
       &'WARNING: Parallel q not tested!', -1)
  IF (elph .and. .not.ep_coupling ) CALL errore('epw_readin', &
      &'elph requires ep_coupling=.true.',1)
  IF ( (elph.and.wannierize) .and. (epwread) ) CALL errore('epw_readin', &
       & 'must use same w90 rotation matrix for entire run', 1)
  IF (wannierize .and. .not.ep_coupling ) CALL errore('epw_readin', &
      &'wannierize requires ep_coupling=.true.',1)
  IF ((wmin.gt.wmax)) &
       CALL errore ('epw_readin', ' check wmin, wmax ', 1)
  IF ((wmin_specfun.gt.wmax_specfun)) &
       CALL errore ('epw_readin', ' check wmin_specfun, wmax_specfun ', 1)
  IF ((nw_specfun.lt.2)) CALL errore ('epw_readin', &
       &' nw_specfun must be at least 2', 1)
  IF ((nqstep.lt.2)) CALL errore ('epw_readin', &
       &' nqstep must be at least 2', 1)
  IF ((nbndsub.gt.200)) CALL errore ('epw_readin', & 
       ' too many wannier functions increase size of projx', 1)
  IF (( phonselfen .or. elecselfen .or. specfun ) .and. ( mp_mesh_k .or. mp_mesh_q )) & 
     CALL errore('epw_readin', 'can only work with full uniform mesh',1)
  IF (ephwrite .and. .not.ep_coupling .and. .not.elph ) CALL errore('epw_readin', &
      &'ephwrite requires ep_coupling=.true., elph=.true.',1)
  IF (ephwrite .and. (rand_k .OR. rand_q ) ) &
     CALL errore('epw_readin', 'ephwrite requires a uniform grid',1) 
  IF (ephwrite .and. (mod(nkf1,nqf1) .ne. 0 .OR. mod(nkf2,nqf2) .ne. 0 .OR. mod(nkf3,nqf3) .ne. 0 ) ) &
     CALL errore('epw_readin', 'ephwrite requires nkf1,nkf2,nkf3 to be multiple of nqf1,nqf2,nqf3',1)
  IF (band_plot .and. filkf .eq. ' ' .and. filqf .eq. ' ') CALL errore('epw_readin', &
      &'plot band structure and phonon dispersion requires k- and q-points read from filkf and filqf files',1)
  IF (band_plot .and. parallel_q ) CALL errore('epw_readin', &
      &'band_plot can only be used with parallel_k',1)    
  IF (band_plot .and. filkf .ne. ' ' .and. (nkf1 > 0 .or. nkf2 > 0 .or. nkf3 > 0) ) CALL errore('epw_readin', &
                &'You should define either filkf or nkf when band_plot = .true.',1)
  IF (band_plot .and. filqf .ne. ' ' .and. (nqf1 > 0 .or. nqf2 > 0 .or. nqf3 > 0) ) CALL errore('epw_readin', &
                &'You should define either filqf or nqf when band_plot = .true.',1)
  IF ( filkf .ne. ' ' .and. .not.efermi_read ) CALL errore('epw_readin', &
      &'WARNING: if k-points are along a line, then efermi_read=.true. and fermi_energy must be given in the input file',-1)
  IF ( scattering .AND. nstemp < 1 ) CALL errore('epw_init', &
       'wrong number of nstemp',1)
  IF ( scattering .AND. maxval(temps(:)) > 0.d0 .AND. tempsmin > 0.d0 .AND. tempsmax > 0.d0 ) &
       CALL errore('epw_init', 'define either (tempsmin and tempsmax) or temps(:)',1)
  IF ( scattering .AND. tempsmax < tempsmin ) &
       CALL errore('epw_init', 'tempsmax should be greater than tempsmin',1)
  IF ( int_mob .AND. efermi_read)  CALL errore('epw_init', &
       'Fermi level can not be set (efermi_read) when computing intrinsic mobilities',1)
  IF ( int_mob .AND. (ABS(ncarrier) > 1E+5) )  CALL errore('epw_init', &
       'You cannot compute intrinsic mobilities and doped mobilities at the same time',1)
  IF ( carrier .AND. (ABS(ncarrier) < 1E+5) )  CALL errore('epw_init', &
       'The absolute value of the doping carrier concentration must be larger than 1E5 cm^-3',1)
  IF ( (iterative_bte .AND. scattering_serta) .OR. (iterative_bte .AND.scattering_0rta)  ) CALL errore('epw_init', &
       'You should first do a run in the SRTA to get the initial scattering_rate files.',1)
  IF ( (longrange .OR. shortrange) .AND. (.not. lpolar)) CALL errore('epw_init',&
       &'Error: longrange or shortrange can only be true if lpolar is true as well.',1)
  IF ( longrange .AND. shortrange) CALL errore('epw_init',&
       &'Error: longrange and shortrange cannot be both true.',1)
  IF ( epwread .AND. .not. kmaps .AND. .not. epbread) CALL errore('epw_init',&
       &'Error: kmaps has to be true for a restart run. ',1)
  !
  ! thickness and smearing width of the Fermi surface  
  ! from eV to Ryd
  fsthick = fsthick / ryd2ev 
  degaussw = degaussw / ryd2ev 
  delta_smear = delta_smear / ryd2ev 
  !
  ! smearing of phonon in a2f
  ! from meV to Ryd 
  degaussq = degaussq / ryd2mev
  delta_qsmear = delta_qsmear / ryd2mev
  !
  ! fermi_energy read from the input file
  ! from eV to Ryd
  IF ( efermi_read ) THEN
     fermi_energy = fermi_energy / ryd2ev
  ENDIF
  ! eptemp : temperature for the electronic Fermi occupations in the e-p calculation (units of Kelvin)
  ! 1 K in eV = 8.6173423e-5
  ! from K to Ryd
  ! Out-of bound issue with GCC compiler. Multiple Fermi temp is not used anyway.
  eptemp = eptemp * kelvin2eV / ryd2ev
  !DO i = 1, ntempxx
  !   IF (eptemp(i) .gt. 0.d0) THEN
  !      ! 1 K in eV = 8.6173423e-5
  !      ! from K to Ryd
  !      eptemp(i) = eptemp(i) * kelvin2eV / ryd2ev
  !   ENDIF
  !ENDDO
  !
  ! from cm-1 to Ryd
  eps_acustic = eps_acustic / ev2cmm1 / ryd2ev 
  !
  !    reads the q point (just if ldisp = .false.)
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
  IF ( scattering ) THEN
    DO i = 1, ntempxx
      IF (temps(i) .gt. 0.d0) THEN
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
  xq(:) = 0.d0
  !
  lgamma = .false.
  tmp_dir = trim(outdir)
  dvscf_dir = trim(dvscf_dir)//'/'
  !
400 continue
  CALL bcast_ph_input
  xqq(:) = xq(:) 
  !
  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !
  modenum_aux = modenum
  !
  ! SP: This initialized nspin and nspin_mag
  IF ( epwread .and. .not. epbread ) then
    continue
  ELSE
    CALL read_file
  ENDIF
  !
  ! nbnd comes out of readfile
  IF (nbndsub.eq.0) nbndsub = nbnd
  !
#if defined(__MPI)
  IF (.not.(me_pool /=0 .or. my_pool_id /=0)) THEN
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
  IF (modenum_aux .ne. -1) THEN
     modenum = modenum_aux
     iswitch = -4
  ELSE IF (modenum .eq. 0) THEN
     iswitch = -2
  ELSE
     iswitch = -4
  END IF
  !
  CALL mp_bcast ( iswitch, meta_ionode_id, world_comm )
  CALL mp_bcast ( modenum, meta_ionode_id, world_comm)
  !
  IF (tfixed_occ) &
     CALL errore('epw_readin','phonon with arbitrary occupations not tested',1)
  !
  IF (elph.and.lsda) CALL errore ('epw_readin', 'El-ph and spin not implemented', 1)
  !
  !   There might be other variables in the input file which describe
  !   partial computation of the dynamical matrix. Read them here
  !
  CALL allocate_part ( nat )
  IF (me_pool /= 0 .or. my_pool_id /=0) goto 800 
  IF (nat_todo < 0 .OR. nat_todo > nat) CALL errore ('epw_readin', &
       'nat_todo is wrong', 1)
  IF (nat_todo.NE.0) THEN
     IF (meta_ionode)read (5, *, iostat = ios) (atomo (na), na = 1, nat_todo)
     CALL mp_bcast(ios, meta_ionode_id, world_comm  )
700  CALL errore ('epw_readin', 'reading atomo', abs (ios) )
     CALL mp_bcast(atomo, meta_ionode_id, world_comm )
  ENDIF
800 continue
  CALL bcast_ph_input1
  !
  DO it = 1, ntyp
     IF (amass (it) <= 0.d0) CALL errore ('epw_readin', 'Wrong masses', it)
  ENDDO
  ! 
! No k+q anymore, can we cut?
  IF (mod (nks, 2) .ne.0) CALL errore ('epw_readin', &
      'k-points are odd', 0)
  !
  !  broadcast the values of nq1, nq2, nq3
  !
  CALL mp_bcast( nq1, meta_ionode_id, world_comm )
  CALL mp_bcast( nq2, meta_ionode_id, world_comm )
  CALL mp_bcast( nq3, meta_ionode_id, world_comm )
  CALL mp_bcast( nk1, meta_ionode_id, world_comm )
  CALL mp_bcast( nk2, meta_ionode_id, world_comm )
  CALL mp_bcast( nk3, meta_ionode_id, world_comm )
  !
  amass = AMU_RY * amass
  !
  END SUBROUTINE epw_readin

!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE kcw_readin()
  !----------------------------------------------------------------------------
  !
  !!  This routine reads the control variables for the KCW programs.
  !!  from standard input (unit 5).
  !!  Then it calls the readfile routine to reads the variables saved
  !!  on a file by the previous self-consistent calculation.
  !
  USE kinds,              ONLY : DP
  USE io_global,          ONLY : ionode_id, stdout
  USE mp,                 ONLY : mp_bcast
  USE klist,              ONLY : nks, nkstot, lgauss, ltetra
  USE control_flags,      ONLY : gamma_only !, lecrpa
  USE lsda_mod,           ONLY : nspin
  USE run_info,           ONLY : title
  USE control_lr,         ONLY : lgamma, lrpa
  USE qpoint,             ONLY : nksq
  USE io_files,           ONLY : tmp_dir, prefix, check_tempdir
  USE noncollin_module,   ONLY : noncolin
  USE read_cards_module,  ONLY : read_cards
  USE io_global,          ONLY : ionode
  USE mp_global,          ONLY : intra_image_comm 
  USE paw_variables,      ONLY : okpaw
  USE uspp,               ONLY : okvan
  USE control_kcw
  USE control_flags,      ONLY : iverbosity
  USE martyna_tuckerman,  ONLY : do_comp_mt
  USE input_parameters,   ONLY : assume_isolated
  USE martyna_tuckerman,  ONLY : do_comp_mt
  USE exx_base,           ONLY : x_gamma_extrapolation
  USE mp_pools,           ONLY : npool
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER             :: ios 
  CHARACTER (LEN=256) :: outdir
  !
  INTEGER, EXTERNAL   :: atomic_number
  REAL(DP), EXTERNAL  :: atom_weight
  LOGICAL, EXTERNAL   :: imatches
  LOGICAL, EXTERNAL   :: has_xml
  LOGICAL             :: exst, parallelfs
  !
  NAMELIST / CONTROL /  outdir, prefix, read_unitary_matrix, kcw_at_ks, &
                        spread_thr, homo_only, kcw_iverbosity, calculation, &
                        l_vcut, assume_isolated, spin_component, mp1, mp2, mp3, lrpa
  !
  NAMELIST / WANNIER /  num_wann_occ, num_wann_emp, have_empty, has_disentangle, &
                        seedname, check_ks, l_unique_manifold
  !
  NAMELIST / SCREEN /   fix_orb, niter, nmix, tr2, i_orb, eps_inf, check_spread
  !
  NAMELIST / HAM /      qp_symm, kipz_corr, i_orb, do_bands, use_ws_distance, & 
                        write_hr, l_alpha_corr, on_site_only
  !
  !### COTROL
  !! outdir          : directory where input, output, temporary files reside 
  !! prefix          : the prefix of files produced by pwscf
  !! read_unitary_matrix  : if true the unitary matrix relating KS to localized orbital is read
  !! kcw_at_ks       : if true compute screening parameter for KS orbitals
  !! homo_only       : if kcw_at_ks only the screening coefficeint for HOMO is computed 
  !! kcw_iverbosity  : level of verbosity. Form 0 (only relevant infos) to 2 (almost everything) 
  !! calculation     : specify the task: wann2kcw, interface with PW and W90; screen: calculation of the screening coefficients;
  !!                   ham, comoute (interpolate) and diagonalize the KC hamiltonian
  !! l_vcut          : IF true the Gygi-Baldereschi scheme is used to deal with the q->0 divergence in the Coulomb integrals
  !! assume_isolated : scheme to deal with the long rage of Coulomnb for isolated systems (as in PWscf)
  !! spin_component  : set the spin component we are looking at (same as the one from Wannir90) 
  !! mp*             : Monhkost-Pack grid, need to be consistent with PW and W90
  !! lrpa            : If true the response of the system is evaluated at the RPA level (no xc contribution) 
  !! spread_thr      : the tollerance within which two orbital are considered to have the same spread 
  !
  !### WANNIER 
  !! seedname        : seedname for the Wannier calculation
  !! num_wann_occ    : number of occupied wannier
  !! num_wann_emp    : number of empty wannier
  !! have_empty      : TRUE if want to compute empty states 
  !! has_disentangle : TRUE if Wannier functions are generated after a disentangle procedure 
  !! l_unique_manifold : TRUE if a unique wannierization was performed for occupied and empty states (use with caution). 
  !! check_ks        : if TRUE a check of the KS hamiltonian build on the Wannier representation is performed (eigenvalues
  !!                   are compared to the original one from PW.
  !
  !### SCREEN 
  !! fix_orb         : if .true. and kcw_at_ks froze the response of the KS we are looking at (FIXME obsolete, to be removed?) 
  !! niter           : max number of iterations for the SCF problem
  !! nmix            : # of iteration used in the mixing scheme for the SCF potential
  !! tr2             : threshold for the convergence of the SCF problem
  !! i_orb           : IF present in input specify which orbital the LR porblem needs to be solved (used to split the calculation
  !!                   of the screening coefficients) 
  !! eps_inf         : The value of the macroscopic dielectric function (needed to deal with the q->0 limit of the screened KC)
  !! check_spread    : If TRUE the screening calculation is performed only for orbital with different spread (self-hartree)
  !
  !### HAM
  !! do_bands        : if .true. KC electronic bands are computed along the input path
  !! use_ws_distance : as in W90, if .true. the Wannier centers are considered in the interpolation
  !! write_hr        : if .true. KC H(R) is printed into a file
  !! on_site_only    : if .true. only H(R=0) and i=j is computed
  !! qp_symm         : if TRUE make the KI hamitonian hermitian in the spirit of quasiparticle GW scheme 
  !! kipz_corr       : Compute the pKIPZ hamiltonian (only for finite systems: Gamma-only calculation in SC) 
  !! l_alpha_corr    : If true a correction is applied to the screening coefficient to mimick effect beyond the 
  !!                   second order
  ! 
  IF (ionode) THEN
    !
    ! ... Input from file ?
    CALL input_from_file ( )
    !
    ! ... Read the first line of the input file
    READ( 5, '(A)', IOSTAT = ios ) title
    !
  ENDIF
  !
  CALL mp_bcast(ios, ionode_id, intra_image_comm )
  CALL errore( 'KC_SCREEN_reading', 'reading title ', ABS( ios ) )
  ! 
  call mp_bcast ( title, ionode_id, intra_image_comm )
  !
  ! Rewind the input if the title is actually the beginning of inputph namelist
  IF( imatches("&control", title) ) THEN
    IF ( ionode ) THEN
      WRITE(*, '(6x,a)') "Title line not specified: using 'default'."
      title='default'
      REWIND(5, iostat=ios)
    ENDIF
    CALL mp_bcast(ios, ionode_id, intra_image_comm)
    CALL errore('KC_readin', 'Title line missing from input.', abs(ios))
  ENDIF
  !
  ! ... set default values for variables in namelist
  !
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  !
  prefix              = 'pwscf'
  kcw_at_ks           =.TRUE.
  lrpa                =.FALSE.
  fix_orb             =.FALSE.
  spread_thr          = 0.001 !(Rydberg)
  homo_only           =.FALSE.
  read_unitary_matrix =.FALSE.
  qp_symm             =.FALSE.
  kipz_corr           =.FALSE. 
  have_empty          =.FALSE.
  has_disentangle     =.FALSE.
  seedname            = 'wann'
  num_wann_occ        = 0
  num_wann_emp        = 0
  check_ks            = .FALSE.
  kcw_iverbosity       = 1
  spin_component      = 1
  niter               = maxter 
  alpha_mix(:)        = 0.D0
  alpha_mix(1)        = 0.7
  nmix                = 4
  tr2                = 1.0d-14
  i_orb               = -1
  mp1                 = -1
  mp2                 = -1
  mp3                 = -1
  do_bands            = .FALSE.
  use_ws_distance     = .TRUE.
  write_hr            = .TRUE.
  eps_inf             = 1.D0
  l_vcut              = .FALSE.
  x_gamma_extrapolation = .FALSE.
  assume_isolated     = 'none'
  l_alpha_corr        = .FALSE. 
  l_unique_manifold   = .FALSE.
  check_spread        = .FALSE.
  on_site_only        = .FALSE.
  calculation         = " " 
  ! 
  ! ...  reading the namelists (if needed)
  !
  IF (ionode) READ( 5, CONTROL, IOSTAT = ios )
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore( 'kcw_readin', 'reading CONTROL namelist', ABS( ios ) )
  CALL mp_bcast(calculation, ionode_id, intra_image_comm)
  !
  IF (calculation /= 'wann2kcw' .AND. calculation /= 'screen' .AND. calculation /= 'ham' .AND. calculation /= 'cc' ) &
  CALL errore('kcw_readin', 'calculation NOT specified or NOT correct', 1)
  !
  IF (ionode .AND. .NOT. kcw_at_ks .AND. (calculation /= 'cc') )  READ( 5, WANNIER, IOSTAT = ios )
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore( 'kcw_readin', 'reading WANNIER namelist', ABS( ios ) )
  !
  IF (ionode .AND. calculation == 'screen') READ( 5, SCREEN, IOSTAT = ios )
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore( 'kcw_readin', 'reading SCREEN namelist', ABS( ios ) )
  ! 
  IF (ionode .AND. calculation == 'ham' ) READ( 5, HAM, IOSTAT = ios )
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore( 'kcw_readin', 'reading HAM namelist', ABS( ios ) )
  !
  IF (kcw_at_ks) seedname = prefix
  IF (ionode) tmp_dir = trimcheck (outdir)
  !
  ! ... broadcasting all input variables to other nodes
  !
  CALL input_summary ()
  CALL bcast_kcw_input ( ) 
  !
  IF ( do_bands ) THEN
    CALL read_cards( 'PW' )
    CALL convert_kpts_names( )
  ENDIF
  !
  tmp_dir_save=tmp_dir
  tmp_dir_kcw= TRIM (tmp_dir) // 'kcw' //'/'
  CALL check_tempdir ( tmp_dir_kcw, exst, parallelfs )
  tmp_dir_kcwq=tmp_dir_kcw
  !
  ! ... Check all namelist variables
  !
  IF (kcw_iverbosity .gt. 1) iverbosity = 1
  !
  IF (spin_component /= 1 .AND. spin_component /= 2) & 
     CALL errore ('kcw_readin', ' spin_component either 1 (UP) or 2 (DOWN) ', 1)
  !
  IF (kcw_at_ks .AND. read_unitary_matrix) THEN 
     CALL infomsg('kcw_readin','WARNING: "read_unitary_matrix" set to FALSE since kcw_at_ks=.TRUE.')
     read_unitary_matrix = .FALSE.
  ENDIF
  !
  IF (fix_orb .AND. .NOT. kcw_at_ks) &
     CALL errore ('kcw_readin', ' fix_orb only works with kcw_at_ks ', 1)
  !
  IF (homo_only .AND. .NOT. kcw_at_ks) & 
     CALL errore ('kcw_readin', ' homo_only only works with kcw_at_ks ', 1)
  !
  IF (has_disentangle .AND. .NOT. have_empty) &
     CALL errore ('kcw_readin', ' disentangle for empty state only ', 1)
  !
  IF (i_orb .lt. -1 ) & 
     CALL errore('kcw_readin', ' WRONG i_orb, orbital from input must be positive', 1)
  !
  IF ( (mp1 .lt. 1 .OR. mp2 .lt. 1 .OR. mp3 .lt. 1) )&
     CALL errore('kcw_readin', ' WRONG k/q grid: check input for mp1, mp2, mp3', 1)
  !
  IF (calculation == 'ham' .AND. npool .gt. 1) &
     CALL errore('kcw_readin', 'pools not implemented for "ham" calculation', npool)
  !
  IF (trim( assume_isolated ) == 'mt' .OR. trim( assume_isolated ) == 'm-t' .OR. trim(assume_isolated) == 'martyna-tuckerman' ) THEN 
    do_comp_mt =.true. 
  ELSE IF ( trim(assume_isolated) == 'none' ) THEN
    do_comp_mt = .false.
  ELSE 
    CALL errore('kcw_readin', ' "assume isolated" not recognized', 1)
  ENDIF
  ! 
  IF (do_comp_mt .AND. mp1*mp2*mp3 /= 1) THEN 
     CALL infomsg('kcw_readin','WARNING: "do_comp_mt" set to FALSE. "l_vcut" set to TRUE instead')
     WRITE(stdout,'()') 
     do_comp_mt =.false.
     l_vcut = .true.
  ENDIF
  !
  ! read data produced by pwscf
  !
  WRITE( stdout, '(5X,"INFO: Reading pwscf data")')
  CALL read_file ( )
  !
  IF ( lgauss .OR. ltetra ) CALL errore( 'kcw_readin', &
      'KC corrections only for insulators!', 1 )
  !
  IF (calculation /= 'wann2kcw' .AND. (mp1*mp2*mp3 /= nkstot/nspin) ) &
     CALL errore('kcw_readin', ' WRONG number of k points from input, check mp1, mp2, mp3', 1)
  !
  IF (gamma_only) CALL errore('kcw_readin',&
     'cannot start from pw.x data file using Gamma-point tricks',1)
  !
  IF (okpaw.or.okvan) CALL errore('kcw_readin',&
     'The KCW code with US or PAW is not available yet',1)
  !
  IF (noncolin) CALL errore('kcw_readin',&
   'The KCW code with non colliner spin is not available yet',1)
  !
  IF ( nspin == 1 ) THEN
     WRITE(stdout, '(/, 5X, "WARNING: nspin=1. A meaningfull KC requires nspin=2 ALWAYS (even for spin-unpolarized systems)")')
     WRITE(stdout, '(   5X, "WARNING: use nspin=1 only if you know what you are doing")') 
  ENDIF
  !
  IF (l_vcut .AND. do_comp_mt ) THEN
     CALL infomsg('kcw_readin','WARNING: "l_vcut" set to FALSE since do_comp_mt=.TRUE.')
     l_vcut = .false.
  ENDIF
  !
  IF (lgamma) THEN
     nksq = nks
  ELSE
     nksq = nks / 2
  ENDIF
  !
  RETURN
  !
END SUBROUTINE kcw_readin

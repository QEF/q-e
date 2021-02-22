!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE kc_readin()
  !----------------------------------------------------------------------------
  !
  !!  This routine reads the control variables for the KC programs.
  !!  from standard input (unit 5).
  !!  Then it calls the readfile routine to reads the variables saved
  !!  on a file by the previous self-consistent calculation.
  !
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : ionode_id, stdout
  USE mp,            ONLY : mp_bcast
  USE klist,         ONLY : xk, nks, nkstot 
  USE control_flags, ONLY : gamma_only !, lecrpa
  USE lsda_mod,      ONLY : nspin
  USE run_info,      ONLY : title
  USE control_ph,    ONLY : lgamma_gamma, ldisp, niter_ph
  USE control_lr,    ONLY : lgamma, lrpa
  USE qpoint,        ONLY : nksq
  USE disp,          ONLY : nq1, nq2, nq3
  USE io_files,      ONLY : tmp_dir, prefix, check_tempdir
  USE noncollin_module, ONLY : noncolin
  !  USE control_flags, ONLY : twfcollect
  USE read_cards_module, ONLY : read_cards
  USE io_global,     ONLY : ionode
  USE mp_global,     ONLY : intra_image_comm 
  USE paw_variables, ONLY : okpaw
  USE uspp,          ONLY : okvan
  USE control_kc_wann
  USE control_ph,    ONLY : tmp_dir_ph, tmp_dir_phq, alpha_mix, nmix_ph, tr2_ph
  USE mp_images,     ONLY : my_image_id
  USE save_ph,       ONLY : tmp_dir_save
  USE control_flags, ONLY : iverbosity
  USE martyna_tuckerman,    ONLY : do_comp_mt
  USE input_parameters,      ONLY : assume_isolated
  USE martyna_tuckerman, ONLY: do_comp_mt
  !USE mp_pools,      ONLY : npool
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ios 
    ! integer variable for I/O control
    ! counter on polarizations
    ! counter on iterations
    ! counter on atoms
    ! counter on types
  CHARACTER (LEN=256) :: outdir
  !
  CHARACTER(LEN=1), EXTERNAL :: capital
  INTEGER, EXTERNAL  :: atomic_number
  REAL(DP), EXTERNAL :: atom_weight
  LOGICAL, EXTERNAL  :: imatches
  LOGICAL, EXTERNAL  :: has_xml
  CHARACTER(LEN=6)   :: int_to_char
  LOGICAL            :: exst, parallelfs
  !
  ! kc_iverbosity   : verbosity control
  ! modenum      : single mode calculation
  ! fildyn       : output file for the dynamical matrix
  ! fildvscf     : output file containing deltavsc
  ! fildrho      : output file containing deltarho
  ! start_q      : in q list does the q points from start_q to last_q
  ! last_q       : 
  ! start_irr    : does the irred. representation from start_irr to last_irr
  ! last_irr     : 
  ! nogg         : if .true. lgamma_gamma tricks are not used
  !
  NAMELIST / CONTROL /  outdir, prefix, read_unitary_matrix, kc_at_ks, &
                        spread_thr, homo_only, kc_iverbosity, &
                        l_vcut, assume_isolated, spin_component
  !
  NAMELIST / WANNIER /  num_wann_occ, num_wann_emp, have_empty, has_disentangle, &
                        seedname, check_ks, l_unique_manifold
  !
  NAMELIST / SCREEN /   lrpa, fix_orb, niter_ph, nmix_ph, tr2_ph, i_orb, mp1, mp2, mp3, eps_inf, check_spread
  !
  NAMELIST / HAM /      qp_symm, kipz_corr, compute_hf, mp1, mp2, mp3, i_orb, &
                        do_bands, use_ws_distance, write_hr, l_alpha_corr, lrpa, on_site_only
  !! CONTROL NAMELIST:
  !
  !! outdir       : directory where input, output, temporary files reside 
  !
  !! prefix       : the prefix of files produced by pwscf
  !
  !! read_unitary_matrix  : if true the unitary matrix relating KS to localized orbital is read
  !
  !! kc_at_ks     : if true compute screening parameter for KS orbitals
  !
  !! fix_orb      : if .true. and kc_at_ks froze the response of the KS we are looking at
  !
  !! homo_only    : if kc_at_ks only the screening coefficeint for HOMO is computed 
  !
  !! kc_iverbosity: level of verbosity. Form 0 (only relevant infos) to 2 (almost everything) 
  !
  !! l_vcut       : IF true the Gygi-Baldereschi scheme is used to deal with the q->0 divergence in the Coulomb integrals
  !!                USE it only for peridic system, for isolated one, use assume_isolated
  !
  !! assume_isolated : for the moment only Martyna-Tuckerman is used.
  !
  !! spin_component : set the spin component we are looking at (same as the one from Wannir90) 
  !
  !! spread_thr   : the tollerance within which two orbital are considered to have the same spread 
  !
  !! qp_symm      : if TRUE make the KI hamitonian hermitian in the spirit of quasiparticle GW scheme 
  !
  !! compute_hf   : if TRUE compute the expectation value of the Fock operator on the orbitals
  !
  !! seedname     : seedname for the Wannier calculation
  !
  !! num_wann_occ : number of occupied wannier
  !
  !! num_wann_emp : number of empty wannier
  !
  !! do_bands     : if .true. KC electronic bands are computed along the input path
  !
  !! use_ws_distance : as in W90, if .true. the Wannier centers are considered in the interpolation
  !
  !! write_hr     : if .true. KC H(R) is printed into a file
  !
  !! on_site_only : if .true. only H(R=0) and i=j. 
  ! 
  IF (ionode) THEN
  !
  ! ... Input from file ?
  !
     CALL input_from_file ( )
  !
  ! ... Read the first line of the input file
  !
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
  IF( imatches("&control", title)) THEN
    WRITE(*, '(6x,a)') "Title line not specified: using 'default'."
    title='default'
    REWIND(5, iostat=ios)
    CALL errore('KC_readin', 'Title line missing from input.', abs(ios))
  ENDIF
  !
  ! ... set default values for variables in namelist
  !
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  !
  prefix              = 'pwscf'
  kc_at_ks            =.TRUE.
  lrpa                =.FALSE.
  fix_orb             =.FALSE.
  spread_thr          = 0.001 !(Rydberg)
  homo_only           =.FALSE.
  read_unitary_matrix =.FALSE.
  qp_symm             =.FALSE.
  kipz_corr           =.FALSE. 
  have_empty          =.FALSE.
  has_disentangle     =.FALSE.
  compute_hf          =.FALSE.
  seedname            = 'wann'
  num_wann_occ        = 0
  num_wann_emp        = 0
  check_ks            = .FALSE.
  kc_iverbosity       = 1
  spin_component      = 1
  niter_ph            = maxter 
  alpha_mix(:)        = 0.D0
  alpha_mix(1)        = 0.7
  nmix_ph             = 4
  tr2_ph              = 1.0d-14
  i_orb               = -1
  mp1                 = -1
  mp2                 = -1
  mp3                 = -1
  do_bands            = .FALSE.
  use_ws_distance     = .TRUE.
  write_hr            = .TRUE.
  eps_inf             = 1.D0
  l_vcut              = .FALSE.
  assume_isolated     = 'none'
  l_alpha_corr        = .FALSE. 
  l_unique_manifold   = .FALSE.
  check_spread        = .FALSE.
  on_site_only        = .FALSE.
  ! 
  ! ...  reading the namelist inputki
  !
  IF (ionode) READ( 5, CONTROL, IOSTAT = ios )
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore( 'kc_readin', 'reading CONTROL namelist', ABS( ios ) )
  !
  IF (ionode .AND. .NOT. kc_at_ks) READ( 5, WANNIER, IOSTAT = ios )
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore( 'kc_readin', 'reading WANNIER namelist', ABS( ios ) )
  !
  IF (ionode .AND. calculation == 'screen') READ( 5, SCREEN, IOSTAT = ios )
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore( 'kc_readin', 'reading SCREEN namelist', ABS( ios ) )
  ! 
  IF (ionode .AND. calculation == 'ham' ) READ( 5, HAM, IOSTAT = ios )
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore( 'kc_readin', 'reading HAM namelist', ABS( ios ) )
  !
  IF (kc_at_ks) seedname = prefix
  IF (ionode) tmp_dir = trimcheck (outdir)
  !
  ! ... broadcasting all input variables to other nodes
  !
  CALL input_summary ()
  CALL bcast_kc_input ( ) 
  !
  IF ( do_bands ) THEN
    CALL read_cards( 'PW' )
    CALL convert_kpts_names( )
  ENDIF
  !
  tmp_dir_save=tmp_dir
  tmp_dir_ph= TRIM (tmp_dir) // '_ph' // TRIM(int_to_char(my_image_id)) //'/'
  CALL check_tempdir ( tmp_dir_ph, exst, parallelfs )
  tmp_dir_phq=tmp_dir_ph
  !
  ! ... Check all namelist variables
  !
  IF (kc_iverbosity .gt. 1) iverbosity = 1
  
  !IF (npool .gt. 1 .AND. calculation == 'wann2kc') &
  !    CALL errore ('kc_readin','Pool Parallelization not implemented for wann to KC interface. Re-run without pools.',1)
  !
  IF (spin_component /= 1 .AND. spin_component /= 2) & 
     CALL errore ('kc_readin', ' spin_component either 1 (UP) or 2 (DOWN) ', 1)
  !
  IF (kc_at_ks .AND. read_unitary_matrix) THEN 
     CALL infomsg('kc_readin','WARNING: "read_unitary_matrix" set to FALSE since kc_at_ks=.TRUE.')
     read_unitary_matrix = .FALSE.
  ENDIF
  !
  IF (fix_orb .AND. .NOT. kc_at_ks) &
     CALL errore ('kc_readin', ' fix_orb only works with kc_at_ks ', 1)
  !
  IF (homo_only .AND. .NOT. kc_at_ks) & 
     CALL errore ('kc_readin', ' homo_only only works with kc_at_ks ', 1)
  !
  IF (has_disentangle .AND. .NOT. have_empty) &
     CALL errore ('kc_readin', ' disentangle for empty state only ', 1)
  !
  IF (i_orb .lt. -1 ) & 
     CALL errore('kc_readin', ' WRONG i_orb, orbital from input must be positive', 1)
  !
  IF (calculation /= 'wann2kc' .AND. (mp1 .lt. 1 .OR. mp2 .lt. 1 .OR. mp3 .lt. 1) )&
     CALL errore('kc_readin', ' WRONG k/q grid mp1, mp2, mp3', 1)
  !
  !
  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !
  !   amass will also be read from file:
  !   save its content in auxiliary variables
  !
  !
  ! read data produced by pwscf
  !
  IF (trim( assume_isolated ) == 'mt' .OR. trim( assume_isolated ) == 'm-t' .OR. trim(assume_isolated) == 'martyna-tuckerman' ) THEN 
    do_comp_mt =.true. 
  ELSE IF ( trim(assume_isolated) == 'none' ) THEN
    do_comp_mt = .false.
  ELSE 
    CALL errore('kc_readin', ' "assume isolated" not recognized', 1)
  ENDIF
  ! 
  IF (do_comp_mt .AND. mp1*mp2*mp3 /= 1) THEN 
     CALL infomsg('kc_readin','WARNING: "do_comp_mt" set to FALSE. "l_vcut" set to TRUE instead')
     WRITE(stdout,'()') 
     do_comp_mt =.false.
     l_vcut = .true.
  ENDIF
  !
  WRITE( stdout, '(5X,"INFO: Reading pwscf data")')
  CALL read_file ( )
  !
  IF (calculation /= 'wann2kc' .AND. (mp1*mp2*mp3 /= nkstot/nspin) ) &
     CALL errore('kc_readin', ' WRONG number of k points from input, check mp1, mp2, mp3', 1)
  !
  IF (gamma_only) CALL errore('kc_readin',&
     'cannot start from pw.x data file using Gamma-point tricks',1)

  IF (okpaw.or.okvan) CALL errore('kc_readin',&
     'The kc_wann code with US or PAW is not available yet',1)

  IF (noncolin) CALL errore('kc_readin',&
   'The kc_wann code with non colliner spin is not available yet',1)

  ! 
  IF (l_vcut .AND. do_comp_mt ) THEN
     CALL infomsg('kc_readin','WARNING: "l_vcut" set to FALSE since do_comp_mt=.TRUE.')
     l_vcut = .false.
  ENDIF
  !
  lgamma_gamma=.FALSE.
  IF (.NOT.ldisp) THEN
     IF (nkstot==1.OR.(nkstot==2.AND.nspin==2)) THEN
        lgamma_gamma=(lgamma.AND.(ABS(xk(1,1))<1.D-12) &
                            .AND.(ABS(xk(2,1))<1.D-12) &
                            .AND.(ABS(xk(3,1))<1.D-12) )
     ENDIF
!     IF (nogg) lgamma_gamma=.FALSE.
     !
     IF (lgamma) THEN
        nksq = nks
     ELSE
        nksq = nks / 2
     ENDIF
  ENDIF
  !
  lgamma_gamma = .FALSE.
  !
  IF (ldisp .AND. (nq1 .LE. 0 .OR. nq2 .LE. 0 .OR. nq3 .LE. 0)) &
       CALL errore('kc_readin','nq1, nq2, and nq3 must be greater than 0',1)
  !
  RETURN
  !
END SUBROUTINE kc_readin

!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE lr_readin
  !-----------------------------------------------------------------------
  !
  !    This routine reads the control variables from standard input (unit 5).
  !    A second routine read_file reads the variables saved to file
  !    by the self-consistent program.
  !
  USE lr_variables
  USE kinds,               ONLY : DP
  USE io_files,            ONLY : tmp_dir, prefix, wfc_dir
  USE lsda_mod,            ONLY : current_spin, nspin
  USE control_flags,       ONLY : twfcollect,use_para_diag
  USE scf,                 ONLY : vltot, v, vrs, vnew, &
                                  & destroy_scf_type
  USE fft_base,            ONLY : dfftp
  USE gvecs,               ONLY : doublegrid
  USE wvfct,               ONLY : nbnd, et, wg
  USE lsda_mod,            ONLY : isk
  USE ener,                ONLY : ef
  USE io_global,           ONLY : ionode, ionode_id
  USE klist,               ONLY : nks, wk, nelec
  USE fixed_occ,           ONLY : tfixed_occ
  USE input_parameters,    ONLY : degauss, nosym, wfcdir, outdir,&
                                  & max_seconds
  USE ktetra,              ONLY : ltetra
  USE realus,              ONLY : real_space, real_space_debug,&
                                  & init_realspace_vars, qpointlist,&
                                  & betapointlist, read_rs_status, newd_r 
  USE funct,               ONLY : dft_is_meta
  USE io_global,           ONLY : stdout
  USE control_flags,       ONLY : tqr, twfcollect
  USE iotk_module
  USE charg_resp,          ONLY : w_T_prefix, omeg, w_T_npol, epsil
  USE mp,                  ONLY : mp_bcast,mp_barrier
  USE mp_global,           ONLY : my_pool_id, intra_image_comm, &
                                  & intra_pool_comm, nproc_image, &
                                  & nproc_pool, nproc_pool_file, &
                                  & nproc_image_file
  USE io_global,           ONLY : ionode, ionode_id
  USE DFUNCT,              ONLY : newd
  USE vlocal,              ONLY : strf

  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER(LEN=256) :: beta_gamma_z_prefix
          ! fine control of beta_gamma_z file
  CHARACTER(LEN=80) :: disk_io
          ! Specify the amount of I/O activities
  INTEGER :: ios, iunout, ierr, ipol
  LOGICAL :: auto_rs
  !
  NAMELIST / lr_input / restart, restart_step, lr_verbosity, prefix, outdir,&
       & test_case_no, wfcdir, disk_io, max_seconds 
  NAMELIST / lr_control / itermax, ipol, ltammd, real_space, real_space_debug,&
       & charge_response, tqr, auto_rs, no_hxc, n_ipol, project
  NAMELIST / lr_post / omeg, beta_gamma_z_prefix, w_T_npol, plot_type, epsil, itermax_int
  !
  auto_rs = .TRUE.
#ifdef __MPI
  IF (ionode) THEN
#endif
     !
     !   Set default values for variables in namelist
     !
     itermax = 500
     restart = .FALSE.
     restart_step = itermax+1
     lr_verbosity = 1
     prefix = 'pwscf'
     disk_io = 'default'
     ltammd = .FALSE.
     ipol = 1
     n_ipol = 1
     no_hxc = .FALSE.
     real_space = .FALSE.
     real_space_debug = 0
     charge_response = 0
     test_case_no = 0
     tqr = .FALSE.
     auto_rs = .TRUE.
     beta_gamma_z_prefix = 'undefined'
     omeg= 0.0_DP
     epsil = 0.0_DP
     w_T_npol = 1
     plot_type = 1
     project = .FALSE.
     max_seconds = 0.D0
     !
     !   Reading the namelist lr_input
     !
     CALL input_from_file( )
     !
     READ (5, lr_input, err = 200, iostat = ios)
200  CALL errore ('lr_readin', 'reading lr_input namelist', ABS (ios) )
     !
     !
     !   Reading the namelist lr_control
     !
     READ (5, lr_control, err = 201, iostat = ios)
201  CALL errore ('lr_readin', 'reading lr_control namelist', ABS (ios) )
     !
     !
     !   Reading the namelist lr_post
     !
     IF (charge_response == 1) THEN
        READ (5, lr_post, err = 202, iostat = ios)
202     CALL errore ('lr_readin', 'reading lr_post namelist', ABS (ios) )
        bgz_suffix = TRIM ( "-stage2.beta_gamma_z." )
        WRITE(stdout,'(/5x,"Prefix of current run is appended by -stage2")')
        IF ( beta_gamma_z_prefix  == 'undefined' ) THEN
           beta_gamma_z_prefix=TRIM(prefix)
        ENDIF
     ELSE
        bgz_suffix = TRIM ( ".beta_gamma_z." )
     ENDIF
     !
     ! The status of the real space flags should be read manually
     !
     ! Do not mess with already present wfc structure
     twfcollect = .FALSE.
     !
     ! Set-up all the dir and suffix variables.
     !
     outdir = trimcheck(outdir)
     tmp_dir = outdir
     IF ( .NOT. TRIM( wfcdir ) == 'undefined' ) THEN
     !
        wfc_dir = trimcheck ( wfcdir )
        !
     ENDIF
     w_T_prefix = TRIM( tmp_dir ) // TRIM( beta_gamma_z_prefix ) // & 
          & ".beta_gamma_z." 
     !
     ierr = 0
     !
     ! Set-up polarization direction(s). 
     !
     IF ( ipol==4 ) THEN
        !
        n_ipol = 3
        LR_polarization=1
        !
     ELSE
        !
        LR_polarization=ipol
        !
     ENDIF
     IF (itermax_int < itermax) itermax_int=itermax
     !
     ! Limited disk_io support: currently only one setting is supported
     !
     SELECT CASE( TRIM( disk_io ) )
     CASE ( 'reduced' )
        !
        lr_io_level = -1
        restart  = .FALSE.
        !
     CASE DEFAULT
        !
        lr_io_level = 1
     END SELECT
     
#ifdef __MPI
  ENDIF
  !
  CALL bcast_lr_input
  CALL mp_bcast(auto_rs, ionode_id)
#endif
  !
  outdir = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
  IF (auto_rs) CALL read_rs_status( outdir, ierr )
  IF (real_space) real_space_debug=99
  IF (real_space_debug > 0) real_space=.TRUE.
  IF (lr_verbosity > 1) THEN
     WRITE(stdout,'(5x,"Status of real space flags: TQR=", L5 ,& 
          &"  REAL_SPACE=", L5)') tqr, real_space
  ENDIF
  !
  !   Now PWSCF XML file will be read, and various initialisations will be done
  !
  CALL read_file()
  !
  !  Make sure all the features used in the PWscf calculation are actually
  !   supported by TDDFPT.
  !
  CALL input_sanity()
  !
  !  Deallocate some variables created with read_file but not used in TDDFPT
  !
  DEALLOCATE( strf )
  CALL destroy_scf_type(vnew)
  !
  !   Re-initialize all needed quantities from the scf run
  !
  current_spin=1
  !
  CALL init_us_1 ( )
  !
  IF (tqr) THEN
     CALL newd_r()
  ELSE
     CALL newd() !OBM: this is for the ground charge density
  ENDIF
  !
  IF ( real_space_debug > 0 ) THEN
     WRITE(stdout,'(/5x,"Real space implementation V.1 D190908",1x)')
     !  !OBM - correct parellism issues
     CALL init_realspace_vars()
     CALL betapointlist()
     WRITE(stdout,'(5X,"Real space initialisation completed")')
  ENDIF
  !
  ! Now put the potential calculated in read_file into the correct place
  ! and deallocate the now redundant associated variables
  !
  CALL set_vrs ( vrs, vltot, v%of_r, 0, 0, dfftp%nnr, nspin, doublegrid )
  DEALLOCATE( vltot )
  CALL destroy_scf_type(v)
  !
  ! Recalculate the weights of the Kohn-Sham orbitals.
  ! (Should this not be a call to weights() to make this !#!
  ! less insulator specific?)                            !#!
  !
  CALL iweights( nks, wk, nbnd, nelec, et, ef, wg, 0, isk)
  !
  !
  IF ( charge_response == 2 ) CALL lr_set_boxes_density()
  !
  !Scalapack related stuff, 
  !
#ifdef __MPI
  use_para_diag = .TRUE.
  CALL check_para_diag( nbnd )
#else
  use_para_diag = .FALSE.
#endif
  RETURN
CONTAINS
  SUBROUTINE input_sanity()
    !-------------------------------------------------------------------------- 
    ! 
    ! This routine aims to gather together all of the input sanity checks
    !  (features enabled in PWscf which are unsupported in TDDFPT) together in
    !   one place.
    !
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !
    !  Charge response mode 1 is the "do Lanczos chains twice, conserve memory"
    !   scheme
    !
    IF (charge_response == 1 .AND. omeg == 0.D0) & 
         & CALL errore ('lr_readin', & 
         & 'omeg must be defined for charge response mode 1', 1 )  
    IF ( project .AND. charge_response /= 1) &
         & CALL errore ('lr_readin', &
         & 'projection is possible only in charge response mode 1', 1 )

    !
    !  Meta-DFT currently not supported by TDDFPT
    !
    IF (dft_is_meta()) &
         & CALL errore( ' iosys ', ' Meta DFT ' // 'not implemented yet', 1 )
    !
    !  Non-insulating systems currently not supported by TDDFPT
    !
    IF ( (ltetra .OR. tfixed_occ .OR. (degauss /= 0.D0)) ) &
         & CALL errore( ' iosys ', ' Linear response calculation ' // &
         & 'not implemented for non-insulating systems', 1 )
    !
    !  Symmetry not supported
    !
    IF ( .NOT. nosym ) &
         & CALL errore( ' iosys ', ' Linear response calculation ' // &
         & 'not implemented with symmetry', 1 )
    !
    !  K-points implemented but still unsupported (use at your own risk!)
    !
    IF ( .NOT. gamma_only ) &
         & CALL errore(' iosys', 'k-point algorithm is not tested yet',1)
    !
    !  Check that either we have the same numebr of procs as the inital PWscf 
    !  run OR that the wavefunctions were gathered into one file at the end of
    !  the PWscf run.
    !
    IF (nproc_image /= nproc_image_file .AND. .NOT. twfcollect)  &
         CALL errore('lr_readin',&
         & 'pw.x run with a different number of processors. &
         & Use wf_collect=.true.',1)
    IF (nproc_pool /= nproc_pool_file .AND. .NOT. twfcollect)  &
         CALL errore('lr_readin',&
         & 'pw.x run with a different number of pools. &
         & Use wf_collect=.true.',1)

    !
  END SUBROUTINE input_sanity

END SUBROUTINE lr_readin

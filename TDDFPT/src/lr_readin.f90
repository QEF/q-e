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
  USE lr_dav_variables
  USE kinds,               ONLY : DP
  USE io_files,            ONLY : tmp_dir, prefix, wfc_dir
  USE lsda_mod,            ONLY : current_spin, nspin
  USE control_flags,       ONLY : twfcollect,use_para_diag
  USE scf,                 ONLY : vltot, v, vrs, vnew, rho, &
                                  & destroy_scf_type
  USE fft_base,            ONLY : dfftp
  USE gvecs,               ONLY : doublegrid
  USE wvfct,               ONLY : nbnd, et, wg, current_k
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
                                  & betapointlist, newd_r 
  USE funct,               ONLY : dft_is_meta
  USE io_global,           ONLY : stdout
  USE control_flags,       ONLY : tqr, twfcollect, ethr
  USE iotk_module
  USE charg_resp,          ONLY : w_T_prefix, omeg, w_T_npol, epsil
  USE mp,                  ONLY : mp_bcast
  USE mp_world,            ONLY : world_comm
  USE mp_global,           ONLY : my_pool_id, intra_image_comm, &
                                  & intra_bgrp_comm, nproc_image, &
                                  & nproc_pool, nproc_pool_file, &
                                  & nproc_image_file, nproc_bgrp, &
                                  & nproc_bgrp_file
  USE io_global,           ONLY : ionode, ionode_id
  USE DFUNCT,              ONLY : newd
  USE vlocal,              ONLY : strf
  USE exx,                 ONLY : ecutfock
#ifdef __ENVIRON
  USE input_parameters,    ONLY : assume_isolated
  USE environ_base,        ONLY : environ_base_init, ir_end
  USE environ_input,       ONLY : read_environ
  USE environ_base,        ONLY : ifdtype, nfdpoint
  USE ions_base,           ONLY : nsp, ityp, zv, tau, nat
  USE cell_base,           ONLY : at, alat, omega, ibrav
  USE solvent_tddfpt,      ONLY : solvent_initbase_tddfpt
  USE environ_init,        ONLY : environ_initions, environ_initcell,      &
                                  environ_clean, environ_initbase,         &
                                  environ_initions_allocate
  USE environ_main,        ONLY : calc_venviron
  USE mp_bands,            ONLY : me_bgrp
  USE plugin_flags,        ONLY : use_environ
#endif

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
  REAL(kind=dp) :: charge
  !
  NAMELIST / lr_input / restart, restart_step ,lr_verbosity, prefix, outdir, test_case_no, wfcdir, disk_io, max_seconds
  NAMELIST / lr_control / itermax, ipol, ltammd, real_space, real_space_debug, charge_response, tqr, auto_rs, no_hxc, n_ipol, &
       & project, scissor, ecutfock, pseudo_hermitian,d0psi_rs
  NAMELIST / lr_post / omeg, beta_gamma_z_prefix, w_T_npol, plot_type, epsil, itermax_int
  namelist / lr_dav / num_eign, num_init, num_basis_max, residue_conv_thr, precondition,dav_debug, reference,single_pole,&
                          &sort_contr, diag_of_h, close_pre,broadening,print_spectrum,start,finish,step,if_check_orth,&
                          &if_random_init,if_check_her,p_nbnd_occ,p_nbnd_virt,poor_of_ram,poor_of_ram2,max_iter,ecutfock,&
	                  &conv_assistant,if_dft_spectrum,no_hxc,d0psi_rs
  !
  auto_rs = .TRUE.
#ifdef __MPI
  IF (ionode) THEN
#endif
     !
     !   Set default values for variables in namelist
     !
     CALL get_env( 'ESPRESSO_TMPDIR', outdir )
     IF ( trim( outdir ) == ' ' ) outdir = './'
     itermax = 500
     restart = .FALSE.
     restart_step = itermax+1
     lr_verbosity = 1
     prefix = 'pwscf'
     disk_io = 'default'
     ltammd = .FALSE.
     d0psi_rs=.false.
     pseudo_hermitian=.true.
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
     eig_dir='./'
     scissor = 0.d0
     ecutfock = -1d0

     ! For lr_dav
     num_eign=1
     num_init=2
     num_basis_max=20
     broadening=0.005
     residue_conv_thr=1.0E-4
     close_pre=1.0E-5
     turn2planb=1.0E-3
     precondition=.true.
     dav_debug=.false.
     reference=0
     single_pole=.false.
     sort_contr=.true.
     print_spectrum=.true.
     start=0.0d0
     finish=1.0d0
     step=0.001d0
     if_check_orth=.false.
     if_check_her=.false.
     if_random_init=.false.
     p_nbnd_occ=10
     p_nbnd_virt=10
     poor_of_ram=.false.
     poor_of_ram2=.false.
     max_iter=100
     conv_assistant=.false.
     if_dft_spectrum=.false.
     !
     ! ------------------------------------------------------
     ! Reading possible plugin arguments -environ -plumed ...
     ! ------------------------------------------------------
     CALL plugin_arguments()
     ! ------------------------------------------------------
     !
     !   Reading the namelist lr_input
     CALL input_from_file( )
     !
     READ (5, lr_input, err = 200, iostat = ios)
200  CALL errore ('lr_readin', 'reading lr_input namelist', ABS (ios) )
     !
     !
     !   Reading the namelist lr_control
     if(.not. davidson) then
       READ (5, lr_control, err = 201, iostat = ios)
201    CALL errore ('lr_readin', 'reading lr_control namelist', ABS (ios) )
     endif

     if(davidson) then
       READ (5, lr_dav, err = 299, iostat = ios)
299    CALL errore ('lr_readin', 'reading lr_dav namelist', ABS (ios) )
     endif
     !
     !
     !   Reading the namelist lr_post
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
  CALL mp_bcast(auto_rs, ionode_id, world_comm)
#endif
  !
  current_k = 1 ! Required for restart runs as this never gets initalised 
  outdir = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
  IF (auto_rs) CALL read_rs_status( outdir, tqr, real_space, ierr )
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
  ! Copy data read from input file (in subroutine "read_input_file") and
  ! stored in modules input_parameters into internal modules of Environ module
  !
#ifdef __ENVIRON
  !
  IF ( use_environ ) THEN
     !
     !!!!!!!!!!!!!!!!!!!!!!!!!!! Initialisation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     ! Copied from PW/src/input.f90
     ! Note: in the routine "environ_base_init" the variable use_environ (from 
     ! environ_base) is defined according to use_environ (from input_parameters).
     ! In the Environ code the variable use_environ (from environ_base) is used.
     !
     ! Warning: There is something strange with the variable 'assume_isolated'!
     ! It is not used currently.
     !
     CALL read_environ( nat, nsp, assume_isolated, ibrav )
     !
     ! Taken from PW/src/init_run.f90
     !
     ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1)) 
     CALL environ_initbase( dfftp%nnr )
     !
     ! Taken from PW/src/electrons.f90
     !
     CALL environ_initions( dfftp%nnr, nat, nsp, ityp, zv, tau, alat )
     CALL environ_initcell( dfftp%nnr, dfftp%nr1, dfftp%nr2, dfftp%nr3, ibrav, omega, alat, at )
     !
     ! Compute additional unperturbed potentials due to the presence of the environment
     ! and add them to the SCF (HXC) potential.
     !
     WRITE( stdout, '(/5x,"Computing and adding the polarization potentials to HXC potential")' )
     !
     CALL calc_venviron( .TRUE., dfftp%nnr, nspin, 0.0d0, rho%of_r, v%of_r(:,1))
     !
     ! Now, once the Environ potentials were computed we can deallocate numerous 
     ! Environ arrays because they are not needed any longer.
     !
     CALL environ_clean( .TRUE. )
     !
     ! Allocations for TDDFPT
     !
     CALL solvent_initbase_tddfpt(ifdtype, nfdpoint, dfftp%nnr)
     !
  ENDIF
  !
#endif
  !
  !   Set wfc_dir - this is done here because read_file sets wfc_dir = tmp_dir
  !   FIXME:,if wfcdir is not present in input, wfc_dir is set to "undefined"
  !   instead of tmp_dir, because of the logic used in the rest of TDDFPT   
  !
  wfc_dir = trimcheck ( wfcdir )
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
  IF ( charge_response == 2 ) CALL lr_set_boxes_density()
  !
  !   Checking
  !
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
    USE fft_base,         ONLY : dffts
    USE paw_variables,    ONLY : okpaw
    USE uspp,             ONLY : okvan
    USE funct,            ONLY : dft_is_hybrid

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
    IF (nproc_bgrp /= nproc_bgrp_file .AND. .NOT. twfcollect)  &
         CALL errore('lr_readin',&
         & 'pw.x run with a different number of band groups. &
         & Use wf_collect=.true.',1)
    !
    ! No taskgroups and EXX.
    !
    IF (dffts%have_task_groups .AND. dft_is_hybrid()) &
         & CALL errore( ' iosys ', ' Linear response calculation ' // &
         & 'not implemented for EXX+Task groups', 1 )
    !
    ! Experimental task groups warning.
    !
    IF (dffts%have_task_groups) CALL infomsg( 'lr_readin','Usage of task &
         &groups with TDDFPT is still experimental. Use at your own risk.' )
    !
    ! No PAW support.
    !
    IF (okpaw) &
         & CALL errore( ' iosys ', ' Linear response calculation ' // &
         & 'not implemented for PAW', 1 )
    !
    ! No USPP+EXX support.
    !
    IF (okvan .AND. dft_is_hybrid()) &
         & CALL errore( ' iosys ', ' Linear response calculation ' // &
         & 'not implemented for EXX+Ultrasoft', 1 )
    !
    RETURN
    !
  END SUBROUTINE input_sanity
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_rs_status( dirname, tqr, real_space, ierr )
    !------------------------------------------------------------------------
    !
    ! This subroutine reads the real space control flags from a pwscf punch card
    ! OBM 2009 - FIXME: should be moved to qexml.f90
    !
      USE iotk_module
      USE io_global,     ONLY : ionode,ionode_id
      USE io_files,      ONLY : iunpun, xmlpun
      USE mp,            ONLY : mp_bcast
      USE mp_images,     ONLY : intra_image_comm
      !
      IMPLICIT NONE
      !
      CHARACTER(len=*), INTENT(in)  :: dirname
      LOGICAL,          INTENT(out) :: tqr, real_space
      INTEGER,          INTENT(out) :: ierr
      !
      !
      IF ( ionode ) THEN
          !
          ! ... look for an empty unit
          !
          CALL iotk_free_unit( iunpun, ierr )
          !
          CALL errore( 'realus->read_rs_status', 'no free units to read real space flags', ierr )
          !
          CALL iotk_open_read( iunpun, FILE = trim( dirname ) // '/' // &
                            & trim( xmlpun ), IERR = ierr )
          !
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         CALL iotk_scan_begin( iunpun, "CONTROL" )
         !
         CALL iotk_scan_dat( iunpun, "Q_REAL_SPACE", tqr )
         CALL iotk_scan_dat( iunpun, "BETA_REAL_SPACE", real_space )
         !
         CALL iotk_scan_end( iunpun, "CONTROL" )
         !
         CALL iotk_close_read( iunpun )
      ENDIF
      CALL mp_bcast( tqr,  ionode_id, intra_image_comm )
      CALL mp_bcast( real_space,    ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE read_rs_status

END SUBROUTINE lr_readin

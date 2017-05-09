!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
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
  !    by PW scf (and PW nscf for EELS).
  !
  USE lr_variables
  USE lr_dav_variables
  USE kinds,               ONLY : DP
  USE io_files,            ONLY : tmp_dir, prefix, wfc_dir
  USE lsda_mod,            ONLY : current_spin, nspin, isk, lsda
  USE control_flags,       ONLY : twfcollect,use_para_diag, &
                                  & tqr, lkpoint_dir, gamma_only, &
                                  & do_makov_payne
  USE scf,                 ONLY : vltot, v, vrs, vnew, &
                                  & destroy_scf_type, rho
  USE fft_base,            ONLY : dfftp, dffts, dtgs
  USE gvecs,               ONLY : doublegrid
  USE wvfct,               ONLY : nbnd, et, wg, current_k
  USE lsda_mod,            ONLY : isk
  USE ener,                ONLY : ef
  USE io_global,           ONLY : ionode, ionode_id, stdout
  USE klist,               ONLY : nks, wk, nelec, lgauss, ltetra
  USE fixed_occ,           ONLY : tfixed_occ
  USE input_parameters,    ONLY : degauss, nosym, wfcdir, outdir,&
                                  & max_seconds, assume_isolated
  USE realus,              ONLY : real_space, real_space_debug,&
                                  & init_realspace_vars, qpointlist,&
                                  & betapointlist
  USE funct,               ONLY : dft_is_meta
  USE iotk_module
  USE charg_resp,          ONLY : w_T_prefix, omeg, w_T_npol, epsil
  USE mp,                  ONLY : mp_bcast
  USE mp_world,            ONLY : world_comm
  USE mp_global,           ONLY : my_pool_id, intra_image_comm, &
                                  & intra_bgrp_comm, nproc_image, &
                                  & nproc_pool, nproc_pool_file, &
                                  & nproc_image_file, nproc_bgrp, &
                                  & nproc_bgrp_file, my_image_id
  USE DFUNCT,              ONLY : newd
  USE vlocal,              ONLY : strf
  USE exx,                 ONLY : ecutfock
  USE martyna_tuckerman,   ONLY : do_comp_mt
  USE esm,                 ONLY : do_comp_esm
  USE qpoint,              ONLY : xq
  USE xml_io_base,         ONLY : create_directory
  USE io_rho_xml,          ONLY : write_scf
  USE noncollin_module,    ONLY : noncolin
  USE mp_bands,            ONLY : ntask_groups
  USE constants,           ONLY : eps4
  USE control_lr,          ONLY : lrpa
#if defined(__ENVIRON)
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
  ! Fine control of beta_gamma_z file
  CHARACTER(LEN=80) :: disk_io
  ! Specify the amount of I/O activities
  INTEGER :: ios, iunout, ierr, ipol
  LOGICAL :: auto_rs
  CHARACTER(LEN=6) :: int_to_char
  !
  NAMELIST / lr_input /   restart, restart_step ,lr_verbosity, prefix, outdir, &
                        & test_case_no, wfcdir, disk_io, max_seconds
  NAMELIST / lr_control / itermax, ipol, ltammd, real_space, real_space_debug, lrpa,   &
                        & charge_response, tqr, auto_rs, no_hxc, n_ipol, project,      &
                        & scissor, ecutfock, pseudo_hermitian, d0psi_rs, lshift_d0psi, &
                        & q1, q2, q3, approximation 
  NAMELIST / lr_post /    omeg, beta_gamma_z_prefix, w_T_npol, plot_type, epsil, itermax_int,sum_rule
  namelist / lr_dav /     num_eign, num_init, num_basis_max, residue_conv_thr, precondition,         &
                        & dav_debug, reference,single_pole, sort_contr, diag_of_h, close_pre,        &
                        & broadening,print_spectrum,start,finish,step,if_check_orth, if_random_init, &
                        & if_check_her,p_nbnd_occ,p_nbnd_virt,poor_of_ram,poor_of_ram2,max_iter,     &
                        & ecutfock, conv_assistant,if_dft_spectrum,no_hxc,d0psi_rs,lshift_d0psi,     &
                        & lplot_drho, vccouple_shift, ltammd
  !
  auto_rs = .TRUE.
  !
#if defined(__MPI)
  IF (ionode) THEN
#endif
     !
     ! Checking for the path to the output directory.
     !
     CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
     IF ( trim( outdir ) == ' ' ) outdir = './'
     !
     ! Set default values for variables in namelist.
     !
     itermax = 500
     restart = .FALSE.
     restart_step = itermax+1
     lr_verbosity = 1
     prefix = 'pwscf'
     disk_io = 'default'
     ltammd = .FALSE.
     d0psi_rs=.false.
     lshift_d0psi = .true.
     pseudo_hermitian=.true.
     ipol = 1
     n_ipol = 1
     no_hxc = .FALSE.      
     lrpa = .false.         
     real_space = .FALSE.
     real_space_debug = 0
     charge_response = 0
     sum_rule = -99
     test_case_no = 0
     tqr = .FALSE.
     auto_rs = .TRUE.
     beta_gamma_z_prefix = 'undefined'
     omeg= 0.0_DP
     epsil = 0.0_DP
     w_T_npol = 1
     plot_type = 1
     project = .FALSE.
     max_seconds = 1.0E+7_DP
     scissor = 0.d0
     ecutfock = -1d0
     !
     ! For EELS
     !
     q1 = 1.0d0         
     q2 = 1.0d0         
     q3 = 1.0d0
     approximation = 'TDDFT'
     !
     ! For lr_dav (Davidson program)
     !
     num_eign=1
     num_init=2
     num_basis_max=20
     broadening=0.005d0
     residue_conv_thr=1.0E-4
     close_pre=1.0E-5
     turn2planb=1.0E-3
     precondition=.true.
     dav_debug=.false.
     reference=0.0d0
     vccouple_shift=0.0d0
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
     lplot_drho=.false.
     !
     ! 
     !   Reading possible plugin arguments (-environ).
     !
     CALL plugin_arguments()
     !
     !   Reading the namelist lr_input
     !
     CALL input_from_file( )
     !
     READ (5, lr_input, err = 200, iostat = ios)
200  CALL errore ('lr_readin', 'reading lr_input namelist', ABS (ios) )
     !
     !   Reading the namelist lr_dav or lr_control
     !
     IF (davidson) THEN
        READ (5, lr_dav, err = 201, iostat = ios)
201     CALL errore ('lr_readin', 'reading lr_dav namelist', ABS (ios) )
     ELSE
        READ (5, lr_control, err = 202, iostat = ios)
202     CALL errore ('lr_readin', 'reading lr_control namelist', ABS (ios) )
     ENDIF
     !
     !   Reading the namelist lr_post (only for optical case)
     !
     IF (charge_response == 1 .AND. .NOT.eels) THEN
        !
        READ (5, lr_post, err = 203, iostat = ios)
203     CALL errore ('lr_readin', 'reading lr_post namelist', ABS (ios) )
        !
        bgz_suffix = TRIM ( "-stage2.beta_gamma_z." )
        WRITE(stdout,'(/5x,"Prefix of current run is appended by -stage2")')
        !
        IF ( beta_gamma_z_prefix  == 'undefined' ) THEN
             beta_gamma_z_prefix = TRIM(prefix)
        ENDIF
        !
     ELSE
        bgz_suffix = TRIM ( ".beta_gamma_z." )
     ENDIF
     !
     ! The status of the real space flags should be read manually
     !
     ! Do not mess with already present wfc structure
     !
     twfcollect = .FALSE.
     !
     ! Set-up all the dir and suffix variables.
     !
     outdir = trimcheck(outdir)
     tmp_dir = outdir
     !
     !IF ( .NOT. TRIM( wfcdir ) == 'undefined' ) THEN
     !   wfc_dir = trimcheck ( wfcdir )
     !ENDIF
     !
     IF (.NOT.eels) THEN
        w_T_prefix = TRIM( tmp_dir ) // &
                   & TRIM( beta_gamma_z_prefix ) // ".beta_gamma_z." 
     ENDIF
     !
     ierr = 0
     !
     IF (eels) THEN
        !
        IF (n_ipol /= 1) THEN
           WRITE(stdout,'(5X,"n_ipol /= 1 is not allowed for EELS.")') 
           WRITE(stdout,'(5X,"Setting n_ipol = 1 ...")') 
           n_ipol = 1
        ENDIF
        IF (ipol /= 1) THEN
           WRITE(stdout,'(5X,"ipol /= 1 is not allowed for EELS.")') 
           WRITE(stdout,'(5X,"Setting ipol = 1 ...")')                
           ipol = 1
        ENDIF
        LR_polarization = 1
        !
     ELSE
        !
        ! Optics: set up polarization direction(s) 
        !
        IF ( ipol==4 ) THEN
           !
           n_ipol = 3
           LR_polarization = 1
           !
        ELSE
           !
           LR_polarization = ipol
           !
       ENDIF
       !
     ENDIF
     !
     IF (itermax_int < itermax) itermax_int = itermax
     !
     ! Limited disk_io support: currently only one setting is supported
     !
     SELECT CASE( TRIM( disk_io ) )
     !
     CASE ( 'reduced' )
        !
        lr_io_level = -1
        restart  = .FALSE.
        !
     CASE DEFAULT
        !
        lr_io_level = 1
        !
     END SELECT
     !
     IF (eels) THEN
        !
        ! Level of approximation in turboEELS.
        !
        SELECT CASE( trim(approximation) )
         !
         CASE ( 'TDDFT' )
           !
           no_hxc = .FALSE.
           lrpa   = .FALSE.
           !
         CASE ( 'IPA' )
           !
           no_hxc = .TRUE.
           lrpa   = .TRUE.
           !
         CASE ( 'RPA_with_CLFE' )
           !
           no_hxc = .FALSE.
           lrpa   = .TRUE.
           !
         CASE DEFAULT
           !
           CALL errore( 'lr_readin', 'Approximation ' // &
                & trim( approximation ) // ' not implemented', 1 )
           !
        END SELECT
        !
        ! We do this trick because xq is used in LR_Modules/dv_of_drho.f90
        ! in the Hartree term ~1/|xq+k|^2
        !
        xq(1) = q1
        xq(2) = q2
        xq(3) = q3
        !
        IF ( (q1.lt.eps4) .AND. (q2.lt.eps4) .AND. (q3.lt.eps4) ) &
           CALL errore( 'lr_readin', 'The transferred momentum |q| is too small, the limit is not implemented.', 1 )
        !
     ENDIF
     !
#if defined(__MPI)
  ENDIF
  !
  CALL bcast_lr_input
  CALL mp_bcast(auto_rs, ionode_id, world_comm)
#endif
  !
  ! Required for restart runs as this never gets initialized.
  !
  current_k = 1     
  !
  outdir = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
  !
  IF (.NOT.eels) THEN
     !
     IF (auto_rs) CALL read_rs_status( outdir, tqr, real_space, ierr )
     IF (real_space) real_space_debug=99
     IF (real_space_debug > 0) real_space=.TRUE.
     IF (lr_verbosity > 1) THEN
        WRITE(stdout,'(5x,"Status of real space flags: TQR=", L5 ,& 
                      &"  REAL_SPACE=", L5)') tqr, real_space
     ENDIF
     !
  ENDIF
  !
  ! EELS: Create a temporary directory for nscf files, and for
  ! writing of the turboEELS restart files.
  !
  IF (eels) THEN
     tmp_dir_lr = TRIM (tmp_dir) // 'tmp_eels/'
     CALL create_directory(tmp_dir_lr)
  ENDIF
  !
  ! EELS: If restart=.true. read the initial information from the file
  ! where the turboEELS code saved its own data (including the data about
  ! the nscf calculation)
  !
  IF (eels .AND. restart) tmp_dir = tmp_dir_lr
  !
  ! Now PWSCF XML file will be read, and various initialisations will be done.
  ! I. Timrov: Allocate space for PW scf variables (EELS: for PW nscf files,
  ! if restart=.true.), read and check them.
  !
  ! Optical case: the variables igk_k and ngk are set up through this path:
  ! read_file -> init_igk.
  ! EELS: the variables igk_k and ngk will be re-set up later (because there 
  ! will be not only poins k but also points k+q) through the path:
  ! lr_run_nscf -> init_run -> hinit0 -> init_igk 
  !
  CALL read_file()
  !
  !   Set wfc_dir - this is done here because read_file sets wfc_dir = tmp_dir
  !   FIXME:,if wfcdir is not present in input, wfc_dir is set to "undefined"
  !   instead of tmp_dir, because of the logic used in the rest of TDDFPT   
  !
  wfc_dir = trimcheck ( wfcdir )
  !
  IF (eels) THEN
     !
     ! Specify the temporary directory.
     !
     tmp_dir = tmp_dir_lr
     !
     ! Copy the scf-charge-density to the tmp_dir (PH/check_initial_status.f90).
     ! Needed for the nscf calculation.
     !
     IF (.NOT.restart) CALL write_scf( rho, nspin )
     !
     ! If a band structure calculation needs to be done, do not open a file
     ! for k point (PH/phq_readin.f90)
     !
     lkpoint_dir = .FALSE.
     !
  ENDIF
  !
  ! Make sure all the features used in the PWscf calculation 
  ! are actually supported by TDDFPT.
  !
  CALL input_sanity()
  !
#if defined(__ENVIRON)
  !
  ! Self-consistent continuum solvation model
  !
  IF ( use_environ ) THEN
     !
     ! Periodic boundary corrections, which were possibly activated in PW
     !
#if defined(__MPI)
  IF (ionode) THEN
#endif
     if (do_makov_payne) then
        assume_isolated = 'makov-payne'
     elseif (do_comp_mt) then
        assume_isolated = 'martyna-tuckerman'
     elseif (do_comp_esm) then
        assume_isolated = 'esm'
     else
        assume_isolated = 'none'
     endif
#if defined(__MPI)
  ENDIF
  CALL mp_bcast(assume_isolated, ionode_id, world_comm)
#endif
     !
     ! Copied from PW/src/input.f90
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
     ! Compute additional unperturbed potentials due to the presence of the
     ! environment and add them to the SCF (HXC) potential.
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
  !  Deallocate some variables created with read_file but not used in TDDFPT
  !
  DEALLOCATE( strf )
  CALL destroy_scf_type(vnew)
  !
  ! Re-initialize all needed quantities from the scf run
  ! I. Timrov: this was already done in read_file.
  current_spin = 1
  !
  ! I. Timrov: The routine init_us_1 was already called in read_file above.
  CALL init_us_1 ( )
  !
  ! I. Timrov: The routine newd was already called in read_file above.
  !
  CALL newd() !OBM: this is for the ground charge density
  !
  IF ( real_space_debug > 0 .AND. .NOT.eels) THEN
     !
     WRITE(stdout,'(/5x,"Real space implementation V.1 D190908",1x)')
     ! OBM - correct parellism issues
     CALL init_realspace_vars()
     CALL betapointlist()
     WRITE(stdout,'(5X,"Real space initialisation completed")')
  ENDIF
  !
  ! Now put the potential calculated in read_file into the correct place
  ! and deallocate the redundant associated variables.
  ! Set the total local potential vrs on the smooth mesh 
  ! adding the scf (Hartree + XC) part and the sum of 
  ! all the local pseudopotential contributions.
  ! vrs = vltot + v%of_r
  !
  CALL set_vrs ( vrs, vltot, v%of_r, 0, 0, dfftp%nnr, nspin, doublegrid )
  !
  DEALLOCATE( vltot )
  CALL destroy_scf_type(v)
  !
  ! Recalculate the weights of the Kohn-Sham orbitals.
  ! (Should this not be a call to weights() to make this 
  ! less insulator specific?)                            
  !
  IF (.NOT.eels) CALL iweights( nks, wk, nbnd, nelec, et, ef, wg, 0, isk)
  !
  IF ( charge_response == 2 .AND. .NOT.eels) CALL lr_set_boxes_density()
  !
  ! Scalapack related stuff.
  !
#if defined(__MPI)
  use_para_diag = .TRUE.
  CALL check_para_diag( nbnd )
#else
  use_para_diag = .FALSE.
#endif
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE input_sanity()
    !-------------------------------------------------------------------------- 
    ! 
    ! This subroutine aims to gather all of the input sanity checks
    ! (features enabled in PWscf which are unsupported in TDDFPT).
    ! Written by Simone Binnie 2011
    ! Modified by Iurii Timrov 2015
    !
    USE paw_variables,    ONLY : okpaw
    USE uspp,             ONLY : okvan
    USE funct,            ONLY : dft_is_hybrid
    USE ldaU,             ONLY : lda_plus_u

    IMPLICIT NONE
    !
    !  Charge response mode 1 is the "do Lanczos chains twice, conserve memory" scheme.
    !
    IF (.NOT.eels) THEN
       !
       IF (charge_response == 1 .AND. ( omeg == 0.D0 .AND. sum_rule == -99 ) ) &
           & CALL errore ('lr_readin', &
           & 'omeg must be defined for charge response mode 1', 1 )
       !
       IF ( project .AND. charge_response /= 1) &
           & CALL errore ('lr_readin', &
           & 'projection is possible only in charge response mode 1', 1 )
       !
    ENDIF
    !
    !  Meta-DFT currently not supported by TDDFPT
    !
    IF (dft_is_meta()) CALL errore( 'lr_readin', 'Meta DFT is not implemented yet', 1 )
    !
    ! Hubbard U is not supported
    !
    IF (lda_plus_u) CALL errore('lr_readin', 'TDDFPT with Hubbard U is not implemented',1)
    !
    !  Tetrahedron method and fixed occupations are not implemented.
    !
    IF (ltetra)                 CALL errore( 'lr_readin', 'ltetra is not implemented', 1 )
    IF (tfixed_occ)             CALL errore( 'lr_readin', 'tfixed_occ is not implemented', 1 )
    !
    ! Some limitations of turboTDDFT (and not of turboEELS).
    !
    IF (.NOT. eels) THEN
       !
       !  Non-insulating systems currently not supported by turboTDDFPT, but
       !  supported by turboEELS.
       !
       IF (lgauss) CALL errore( 'lr_readin', 'turboTDDFT is not extended to metals', 1 )
       !
       ! Symmetry is not supported.
       !
       IF (.NOT.nosym ) CALL errore( 'lr_readin', 'Linear response calculation' // &
                                    & 'is not implemented with symmetry', 1 )
       !
       ! K-points are implemented but still unsupported (use at your own risk!)
       !
       IF (.NOT. gamma_only ) CALL errore('lr_readin', 'k-point algorithm is not tested yet',1)
       !
    ENDIF
    !
    !  Check that either we have the same number of procs as the initial PWscf run 
    !  OR that the wavefunctions were gathered into one file at the end of
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
    IF (dtgs%have_task_groups .AND. dft_is_hybrid()) &
         & CALL errore( 'lr_readin', ' Linear response calculation ' // &
         & 'not implemented for EXX+Task groups', 1 )
    !
    ! Experimental task groups warning.
    !
    IF (dtgs%have_task_groups) &
         & CALL infomsg( 'lr_readin','Usage of task &
         &groups with TDDFPT is still experimental. Use at your own risk.' )
    !
    ! No PAW support.
    !
    IF (okpaw) &
         & CALL errore( 'lr_readin', ' Linear response calculation ' // &
         & 'not implemented for PAW', 1 )
    !
    ! No USPP+EXX support.
    !
    IF (okvan .AND. dft_is_hybrid()) &
         & CALL errore( 'lr_readin', ' Linear response calculation ' // &
         & 'not implemented for EXX+Ultrasoft', 1 )
    !
    ! Spin-polarised case is not implemented, but partially accounted in
    ! some routines.
    !
    IF (lsda) CALL errore( 'lr_readin', 'LSDA is not implemented', 1 )
    !
    ! EELS-related restrictions
    !
    IF (eels) THEN
       !
       IF (okvan .AND. noncolin) CALL errore( 'lr_readin', 'Ultrasoft PP + noncolin is not fully implemented', 1 )
       IF (gamma_only)  CALL errore( 'lr_readin', 'gamma_only is not supported', 1 )
       !
       ! Tamm-Dancoff approximation is not recommended to be used with EELS, and
       ! thus it was not implemented.
       !
       IF (ltammd)      CALL errore( 'lr_readin', 'EELS + Tamm-Dancoff approximation is not supported', 1 )
       !
       IF (project)     CALL errore( 'lr_readin', 'project is not allowed', 1 )
       IF (tqr)         CALL errore( 'lr_readin', 'tqr is not supported', 1 )
       IF (real_space)  CALL errore( 'lr_readin', 'real_space is not supported', 1 )
       IF (charge_response /= 0) CALL errore( 'lr_readin', 'charge_response /= 0 is not allowed', 1 )
       IF (dft_is_hybrid())    CALL errore( 'lr_readin', 'EXX is not supported', 1 )
       IF (do_comp_mt)  CALL errore( 'lr_readin', 'Martyna-Tuckerman PBC is not supported.', 1 )
       IF (d0psi_rs)    CALL errore( 'lr_readin', 'd0psi_rs is not allowed', 1 )
       !
       ! Note, all variables of the turboDavidson code cannot be used by turboEELS.
       !
#if defined(__ENVIRON)
       !
       ! EELS + implicit solvent model is not supported.
       !
       IF ( use_environ ) CALL errore( 'lr_readin', 'Implicit solvent model cannot be used.', 1)
#endif
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE input_sanity
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_rs_status( dirname, tqr, real_space, ierr )
    !------------------------------------------------------------------------
    !
    ! This subroutine reads the real space control flags from a PWscf punch card
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

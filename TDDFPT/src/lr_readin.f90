!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
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
  USE io_files,            ONLY : tmp_dir, prefix, wfc_dir, create_directory
  USE lsda_mod,            ONLY : current_spin, nspin, isk, lsda
  USE control_flags,       ONLY : use_para_diag, tqr, gamma_only,&
                                  & do_makov_payne, noinv
  USE scf,                 ONLY : vltot, v, vrs, vnew, &
                                  & destroy_scf_type, rho
  USE fft_base,            ONLY : dfftp, dffts
  USE gvect,               ONLY : gcutm
  USE gvecs,               ONLY : doublegrid
  USE wvfct,               ONLY : nbnd, et, wg, current_k
  USE lsda_mod,            ONLY : isk
  USE ener,                ONLY : ef
  USE io_global,           ONLY : ionode, ionode_id, stdout, meta_ionode, meta_ionode_id
  USE klist,               ONLY : nks, wk, nelec, lgauss, ltetra
  USE fixed_occ,           ONLY : tfixed_occ
  USE symm_base,           ONLY : nosym
  USE cell_base,           ONLY : at, alat
  USE ions_base,           ONLY : nat, tau
  USE check_stop,          ONLY : max_seconds
  USE realus,              ONLY : real_space, init_realspace_vars, generate_qpointlist, &
                                  betapointlist
  USE xc_lib,              ONLY : xclib_dft_is
  USE charg_resp,          ONLY : w_T_prefix, omeg, w_T_npol, epsil
  USE mp,                  ONLY : mp_bcast
  USE mp_global,           ONLY : my_pool_id, intra_image_comm, &
                                  & intra_bgrp_comm, nproc_image, &
                                  & nproc_pool, nproc_pool_file, &
                                  & nproc_image_file, nproc_bgrp, &
                                  & nproc_bgrp_file, my_image_id
  USE DFUNCT,              ONLY : newd
  USE vlocal,              ONLY : strf
  USE martyna_tuckerman,   ONLY : do_comp_mt
  USE esm,                 ONLY : do_comp_esm
  USE qpoint,              ONLY : xq
  USE io_rho_xml,          ONLY : write_scf
  USE mp_bands,            ONLY : ntask_groups
  USE constants,           ONLY : eps4, rytoev
  USE control_lr,          ONLY : lrpa, alpha_mix, ethr_nscf
  USE mp_world,            ONLY : world_comm
#if defined (__ENVIRON)
  USE plugin_flags,          ONLY : use_environ
  USE environ_base_module,   ONLY : read_environ_input, init_environ_setup, &
                                    init_environ_base, update_environ_ions, &
                                    update_environ_cell, clean_environ, &
                                    check_environ_compatibility
  USE environ_pw_module,     ONLY : update_environ_potential, calc_environ_potential
#endif

  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: wfcdir = 'undefined', outdir
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER(LEN=256) :: beta_gamma_z_prefix
  ! Fine control of beta_gamma_z file
  CHARACTER(LEN=80) :: disk_io
  ! Specify the amount of I/O activities
  CHARACTER(LEN=6) :: int_to_char
  INTEGER :: ios, iunout, ierr
  !
  CHARACTER(LEN=80)          :: card
  INTEGER :: i
  !
  NAMELIST / lr_input /   restart, restart_step ,lr_verbosity, prefix, outdir, &
                        & test_case_no, wfcdir, disk_io, max_seconds
  NAMELIST / lr_control / itermax, ipol, ltammd, lrpa,   &
                        & charge_response, no_hxc, n_ipol, project,      &
                        & scissor, pseudo_hermitian, d0psi_rs, lshift_d0psi, &
                        & q1, q2, q3, approximation, calculator, alpha_mix, start, &
                        & end, increment, epsil, units, ethr_nscf, force_real_gamma, &
                        & force_real_alpha, force_zero_alpha, lan_precondition 
  NAMELIST / lr_post /    omeg, beta_gamma_z_prefix, w_T_npol, plot_type, epsil, itermax_int,sum_rule
  namelist / lr_dav /     num_eign, num_init, num_basis_max, residue_conv_thr, precondition,         &
                        & reference,single_pole, sort_contr, diag_of_h, close_pre,        &
                        & broadening,print_spectrum,start,finish,step, if_random_init, &
                        & p_nbnd_occ,p_nbnd_virt,poor_of_ram,poor_of_ram2,max_iter,     &
                        & conv_assistant,if_dft_spectrum,no_hxc,d0psi_rs,lshift_d0psi,               &
                        & lplot_drho, vccouple_shift, ltammd
  !
#if defined (__ENVIRON)
  REAL(DP) :: at_scaled(3, 3)
  REAL(DP) :: gcutm_scaled
  REAL(DP), ALLOCATABLE :: tau_scaled(:, :)
#endif
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
     charge_response = 0
     sum_rule = -99
     test_case_no = 0
     beta_gamma_z_prefix = 'undefined'
     omeg= 0.0_DP
     w_T_npol = 1
     plot_type = 1
     project = .FALSE.
     max_seconds = 1.0E+7_DP
     scissor = 0.d0
     ethr_nscf = 1.D-11
     !
     ! For EELS
     !
     q1 = 1.0d0
     q2 = 1.0d0
     q3 = 1.0d0
     approximation = 'TDDFT'
     calculator = 'lanczos'
     !
     ! Sternheimer
     !
     start_freq=1
     last_freq=0
     alpha_mix(:) = 0.0D0
     alpha_mix(1) = 0.7D0
     units = 0
     end = 2.5D0
     epsil = 0.02D0
     increment = 0.001D0
     !
     ! For Magnons
     !
     force_real_gamma = .FALSE.
     force_real_alpha = .FALSE.
     force_zero_alpha = .FALSE.
     lan_precondition = .FALSE.
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
     reference=0.0d0
     vccouple_shift=0.0d0
     single_pole=.false.
     sort_contr=.true.
     print_spectrum=.true.
     start=0.0d0
     finish=1.0d0
     step=0.001d0
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
     !   Reading possible plugin arguments
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
     IF (charge_response == 1 .AND. .NOT.eels .AND. .NOT. magnons) THEN
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
     ! Set-up all the dir and suffix variables.
     !
     outdir = trimcheck(outdir)
     tmp_dir = outdir
     !
     !IF ( .NOT. TRIM( wfcdir ) == 'undefined' ) THEN
     !   wfc_dir = trimcheck ( wfcdir )
     !ENDIF
     !
     IF (.NOT.eels .AND. .NOT. magnons) THEN
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
        ! Set up polarization direction(s) for the electric field (optics)
        ! or for the magnetic field (magnons)
        !  
        IF (ipol==4) THEN
           n_ipol = 3
           LR_polarization = 1
        ELSEIF (ipol==1 .OR. ipol==2 .OR. ipol==3) THEN
           n_ipol = 1
           LR_polarization = ipol
        ELSE
           CALL errore( 'lr_readin', 'ipol must be 1, 2, 3, or 4',1)
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
     IF (magnons) THEN
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
     IF (davidson) THEN
        !
        ! check and set num_init and num_basis_max
        !
        IF (num_init < num_eign ) THEN
           WRITE(stdout,'(5X,"num_init is too small, set to num_init = 2*num_eign")')
           num_init = 2 * num_eign
        ENDIF
        IF (num_basis_max < 2*num_init ) THEN
           WRITE(stdout,'(5X,"num_basis_max is too small, set to num_basis_max = 4*num_init")')     
           num_basis_max = 4 * num_init
        ENDIF   
        !
     ENDIF
     !
#if defined(__MPI)
  ENDIF
  !
  CALL bcast_lr_input
  !
#endif
  !
  IF ( trim(calculator)=='sternheimer' ) THEN
     nfs=0
     nfs = ((end-start) / increment) + 1
     if (nfs < 1) call errore('lr_readin','Too few frequencies',1)
     ALLOCATE(fiu(nfs))
     ALLOCATE(fru(nfs))
     ALLOCATE(comp_f(nfs))
     comp_f=.TRUE.
     IF (units == 0) THEN
        fru(1) =start
        fru(nfs) = end
        deltaf = increment
     ELSEIF (units == 1) THEN 
        fru(1) =start/rytoev
        fru(nfs) = end/rytoev
        deltaf = increment/rytoev
     ENDIF
     fiu(:) = epsil
     DO i=2,nfs-1
        fru(i)=fru(1) + (i-1) * deltaf
     ENDDO
     IF (start_freq<=0) start_freq=1
     IF (last_freq<=0.OR.last_freq>nfs) last_freq=nfs
  ELSE
     nfs=1
     ALLOCATE(fru(1))
     ALLOCATE(fiu(1))
     ALLOCATE(comp_f(1))
     fru=0.0_DP
     fiu=0.0_DP
     comp_f=.TRUE.
  END IF
  !
  ! Required for restart runs as this never gets initialized.
  !
  current_k = 1
  !
  outdir = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
  !
  ! EELS: Create a temporary directory for nscf files, and for
  ! writing of the turboEELS restart files.
  !
  IF (eels) THEN
     tmp_dir_lr = TRIM (tmp_dir) // 'tmp_eels/'
     CALL create_directory(tmp_dir_lr)
  ENDIF
  ! same for magnons
  IF (magnons) THEN
     tmp_dir_lr = TRIM (tmp_dir) // 'tmp_magnons/'
     CALL create_directory(tmp_dir_lr)
  ENDIF
  !
  ! Now PWSCF XML file will be read, and various initialisations will be done.
  ! Allocate space for PW scf variables, read and check them.
  ! Optical case: the variables igk_k and ngk are set up through this path:
  ! read_file -> init_igk.
  ! EELS: the variables igk_k and ngk will be re-set up later (because there
  ! will be not only poins k but also points k+q) through the path:
  ! lr_run_nscf -> init_run -> allocate_wfc_k -> init_igk
  !
  CALL read_file()
  !
  IF (.NOT.eels .AND. .NOT. magnons .AND. (tqr .OR. real_space)) &
     WRITE(stdout,'(/5x,"Status of real space flags: TQR=", L5 , &
                      & "  REAL_SPACE=", L5)') tqr, real_space
  !
  !   Set wfc_dir - this is done here because read_file sets wfc_dir = tmp_dir
  !   FIXME:,if wfcdir is not present in input, wfc_dir is set to "undefined"
  !   instead of tmp_dir, because of the logic used in the rest of TDDFPT
  !
  wfc_dir = trimcheck ( wfcdir )
  !
  IF (eels .OR. magnons) THEN
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
  ENDIF
  !
  ! Make sure all the features used in the PWscf calculation
  ! are actually supported by TDDFPT.
  !
  CALL input_sanity()
  !
  ! Compute and add plugin contributions to SCF potential
  !
#if defined (__ENVIRON)
  IF (use_environ) THEN
     at_scaled = at * alat
     gcutm_scaled = gcutm / alat**2
     ALLOCATE (tau_scaled(3, nat))
     tau_scaled = tau * alat
     CALL read_environ_input()
     CALL init_environ_setup('TD')
     CALL init_environ_base(at_scaled, gcutm_scaled, do_comp_mt)
     CALL update_environ_ions(tau_scaled)
     DEALLOCATE (tau_scaled)
     CALL update_environ_cell(at_scaled)
     CALL update_environ_potential(v%of_r(:, 1))
     CALL calc_environ_potential(rho, .FALSE., -1.D0, v%of_r(:, 1))
     CALL clean_environ('TD', .FALSE.)
  END IF
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
  IF (.NOT.eels .AND. .NOT. magnons) CALL iweights( nks, wk, nbnd, nelec, et, ef, wg, 0, isk)
  !
  IF ( charge_response == 2 .AND. .NOT.eels .AND. .NOT. magnons) CALL lr_set_boxes_density()
  !
  ! Scalapack related stuff.
  !
  CALL set_para_diag( nbnd, use_para_diag )
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
    USE xc_lib,           ONLY : xclib_dft_is
    USE ldaU,             ONLY : lda_plus_u
    USE noncollin_module, ONLY : domag

    IMPLICIT NONE
    !
    !  Charge response mode 1 is the "do Lanczos chains twice, conserve memory" scheme.
    !
    IF (.NOT.eels .AND. .NOT. magnons) THEN
       !
       IF (charge_response == 1 .AND. ( omeg == 0.D0 .AND. sum_rule == -99 ) ) &
           & CALL errore ('lr_readin', &
           & 'omeg must be defined for charge response mode 1', 1 )
       !
       IF ( project .AND. charge_response /= 1) &
           & CALL errore ('lr_readin', &
           & 'projection is possible only in charge response mode 1', 1 )
       !
       IF (gamma_only) THEN
          nosym=.true.
          WRITE(stdout,*) "Symmetries are disabled for the gamma_only case"
       ENDIF
       !
    ENDIF
    !
    !  Meta-DFT currently not supported by TDDFPT
    !
    IF (xclib_dft_is('meta')) CALL errore( 'lr_readin', 'Meta DFT is not implemented yet', 1 )
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
    IF (.NOT. eels .and. .NOT. magnons) THEN
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
    IF (eels .AND. domag) CALL errore('lr_readin', 'EELS for magnetic systems is not implemented',1)
    !
    ! No taskgroups and EXX.
    !
    IF (dffts%has_task_groups .AND. xclib_dft_is('hybrid')) &
         & CALL errore( 'lr_readin', ' Linear response calculation ' // &
         & 'not implemented for EXX+Task groups', 1 )
    !
    ! Experimental task groups warning.
    !
    IF (dffts%has_task_groups) &
         & CALL infomsg( 'lr_readin','Usage of task &
         & groups with TDDFPT is still experimental. Use at your own risk.' )
    !      & CALL errore( 'lr_readin', ' Linear response calculation ' // &
    !      & 'not implemented for task groups', 1 )
    !
    ! No PAW support.
    !
    IF (okpaw) &
         & CALL errore( 'lr_readin', ' Linear response calculation ' // &
         & 'not implemented for PAW', 1 )
    !
    ! No USPP+EXX support.
    !
    IF (okvan .AND. xclib_dft_is('hybrid')) &
         & CALL errore( 'lr_readin', ' Linear response calculation ' // &
         & 'not implemented for EXX+Ultrasoft', 1 )
    !
    ! Spin-polarised case is not implemented, but partially accounted in
    ! some routines.
    !
    IF (lsda) CALL errore( 'lr_readin', 'LSDA is not implemented', 1 )
    !
    IF (real_space)  THEN
       IF (eels .OR. magnons) THEN
          CALL errore( 'lr_readin', 'Option real_space=.true. is not implemented', 1 )
       ELSE
          CALL errore( 'lr_readin', 'Option real_space=.true. is not tested', 1 )
       ENDIF
    ENDIF  
    !
    ! EELS-related restrictions
    !
    IF (eels .OR. magnons) THEN
       !
       IF (gamma_only)  CALL errore( 'lr_readin', 'gamma_only is not supported', 1 )
       !
       ! Tamm-Dancoff approximation is not recommended to be used with EELS, and
       ! thus it was not implemented.
       !
       IF (ltammd)      CALL errore( 'lr_readin', 'EELS + Tamm-Dancoff approximation is not supported', 1 )
       !
       IF (project)     CALL errore( 'lr_readin', 'project is not allowed', 1 )
       IF (tqr)         CALL errore( 'lr_readin', 'tqr is not supported', 1 )
       IF (charge_response /= 0) CALL errore( 'lr_readin', 'charge_response /= 0 is not allowed', 1 )
       IF (xclib_dft_is('hybrid'))    CALL errore( 'lr_readin', 'EXX is not supported', 1 )
       IF (do_comp_mt)  CALL errore( 'lr_readin', 'Martyna-Tuckerman PBC is not supported.', 1 )
       IF (d0psi_rs)    CALL errore( 'lr_readin', 'd0psi_rs is not allowed', 1 )
       !
       ! Note, all variables of the turboDavidson code cannot be used by turboEELS.
       !
       ! EELS + plugins is not supported.
       !
#if defined (__ENVIRON)
       IF (use_environ) CALL check_environ_compatibility('lr_readin')
#endif
       !
    ENDIF
    !
    ! MAgnons restrictions
    !
    IF (magnons) THEN
       IF (okvan.OR.okpaw) &     
          CALL errore ('lr_readin', ' Magnons linear response calculation ' // &
                      & 'not implemented for USPP and PAW', 1 )            
       IF (xclib_dft_is('gradient')) &
          call errore('lr_readin', 'Magnons linear response calculation ' // &
                     & 'does not support GGA', 1 )
       IF ( (.not. noinv) .or. (.not. nosym)) THEN
          call errore('lr_readin', 'Magnons linear response calculation ' // &
                     & 'is not implemented with symmetry', 1 )
       ENDIF
       IF (.not. domag) &
          CALL errore ('lr_readin', ' Magnons linear response calculation ' // &
                      & 'non-magnetic system', 1 )
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE input_sanity
  !
END SUBROUTINE lr_readin

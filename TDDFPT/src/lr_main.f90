!
! Copyright (C) 2004-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------

PROGRAM lr_main
  !---------------------------------------------------------------------
  ! Brent Walker, ICTP, 2004
  ! Dario Rocca, SISSA, 2006
  ! O. Baris Malcioglu, SISSA, 2008
  !---------------------------------------------------------------------
  ! ... overall driver routine for applying lanczos algorithm
  ! ... to the matrix of equations coming from tddft
  ! ... spectrum version
  !---------------------------------------------------------------------
  !
  !
  USE lr_lanczos,            ONLY : one_lanczos_step
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : dp
  USE lr_variables,          ONLY : restart, restart_step,&
       itermax, lr_verbosity,&
       evc1, evc1_old,norm0, charge_response,&
       n_ipol, d0psi, &
       LR_iteration, LR_polarization, &
       plot_type, no_hxc, nbnd_total, project, F,R, &
       itermax_int, revc0, lr_io_level, code
  USE io_files,              ONLY : nd_nmbr
  USE global_version,        ONLY : version_number
  USE charg_resp,            ONLY : lr_calc_w_T, read_wT_beta_gamma_z, &
       lr_dump_rho_tot_compat1, lr_dump_rho_tot_cube,&
       lr_dump_rho_tot_xyzd,lr_dump_rho_tot_xcrys,&
       lr_dump_rho_tot_pxyd,chi,lr_calc_R,w_t_norm0_store,&
       resonance_condition, lr_dump_rho, lr_calc_project, lr_project_init
  USE ions_base,             ONLY : tau,nat,atm,ityp
  USE environment,           ONLY: environment_start
  USE mp_global,             ONLY : nimage, mp_startup, init_index_over_band, inter_bgrp_comm

  USE wvfct,                 ONLY : nbnd
  USE wavefunctions_module,  ONLY : psic
  USE control_flags,         ONLY : tddfpt
  USE check_stop,            ONLY : check_stop_now, check_stop_init
  USE funct,                 ONLY : dft_is_hybrid
  USE fft_base,              ONLY : dffts
#ifdef __ENVIRON
  USE environ_base,          ONLY : do_environ
  USE environ_info,          ONLY : environ_summary
#endif

  !Debugging
  USE lr_variables, ONLY: check_all_bands_gamma, check_density_gamma,check_vector_gamma
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER            :: ip,pol_index,ibnd_occ,ibnd_virt,ibnd
  INTEGER            :: iter_restart,iteration
  LOGICAL            :: rflag, nomsg, tg_tmp
  COMPLEX(kind=dp)   :: sum_F,sum_c,temp
  !
  !
  pol_index=1
  !
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( code )
  !
  CALL start_clock('lr_main')
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_main>")')
  ENDIF
  !
  ! Let the PHonon and Environ routines know that 
  ! they are doing tddfpt.
  !
  tddfpt=.TRUE.
  !
  !   Reading input file and PWSCF xml, some initialisation
  !
  CALL lr_readin ( )
  !
  ! Writing a summary to the standard output 
  ! about Environ variables
  !
#ifdef __ENVIRON
  IF ( do_environ ) CALL environ_summary()
#endif
  !
  CALL check_stop_init()
  !
  !   Initialisation of degauss/openshell related stuff
  CALL lr_init_nfo()
  !
  !   Allocate and zero lr variables
  CALL lr_alloc_init()
  !
  !   Now print some preamble info about the run to stdout
  CALL lr_print_preamble()
  !
  !   Read in ground state wavefunctions
  CALL lr_read_wf()
  !
  CALL init_index_over_band(inter_bgrp_comm,nbnd)
  !
  tg_tmp = dffts%have_task_groups
  !   Set up initial response orbitals
  IF ( test_restart(1) ) THEN
     CALL lr_read_d0psi()
  ELSE
     CALL lr_solve_e()
  ENDIF

  !do ip = 1, n_ipol
  !  temp=wfc_dot(ibnd)
  !enddo

  
  dffts%have_task_groups = tg_tmp
  !
  DEALLOCATE( psic )
  !
  IF(project) THEN
     call lr_project_init
  ENDIF

  !
  !   Set up initial stuff for derivatives
  CALL lr_dv_setup()
  !
  !   Coordinates of the read atom, just in case
  IF (lr_verbosity > 1) THEN
     WRITE(stdout,'(/,5X,"Positions of atoms in internal coordinates")')
     DO ip=1,nat ! I am using ip here as a counter over atoms
        WRITE(stdout,'(5X,A3,2X,3F15.8)') atm(ityp(ip)),tau(1:3,ip)
     ENDDO
  ENDIF
  !
  !   Lanczos loop where the real work happens
  DO ip=1, n_ipol
     IF (n_ipol/=1) THEN
        LR_polarization=ip
        pol_index=LR_polarization
     ENDIF
     IF (charge_response == 1 ) THEN
        ! Read precalculated beta gamma z
        CALL read_wT_beta_gamma_z()
        CALL lr_calc_w_T()
     ENDIF
     IF (test_restart(2)) THEN
        CALL lr_restart(iter_restart,rflag)
        CALL lr_read_d0psi()
        WRITE(stdout,'(/5x,"Restarting Lanczos loop",1x,i8)') LR_polarization
     ELSE
        CALL lr_read_d0psi()
        evc1(:,:,:,1) = d0psi(:,:,:,pol_index)
        evc1(:,:,:,2) = d0psi(:,:,:,pol_index)
        evc1_old(:,:,:,1) = cmplx(0.0d0,0.0d0)
        evc1_old(:,:,:,2) = cmplx(0.0d0,0.0d0)
        ! The new structure doesn't need it here
        ! CALL lr_normalise( evc1(:,:,:,1), norm0(pol_index) ) 
        ! evc1(:,:,:,2) = evc1(:,:,:,1)
        !
        iter_restart=1
        !
        WRITE(stdout,'(/5x,"Starting Lanczos loop",1x,i8)')   LR_polarization
     ENDIF
     !
     CALL sd0psi() !after this d0psi is Sd0psi 
     !OBM:Check if this is really necessary
     lancz_loop1 : DO iteration = iter_restart, itermax
        LR_iteration=iteration
        WRITE(stdout,'(/5x,"Lanczos iteration:",1x,i8,3x"Pol:",i1,i8)') LR_iteration, ip
        CALL one_lanczos_step()
        IF ( lr_io_level > 0 .and. (mod(LR_iteration,restart_step)==0 .or. &
             LR_iteration==itermax .or. LR_iteration==1) )&
             CALL lr_write_restart()
        !
        ! Check to see if the wall time limit has been exceeded.
        ! if it has exit gracefully saving the last set of Lanczos
        ! iterations.
        IF ( check_stop_now() ) THEN
           CALL lr_write_restart()
           !   Deallocate pw variables
           CALL clean_pw( .FALSE. )
           CALL stop_clock('lr_main')
           CALL print_clock_lr()
           CALL stop_lr( .FALSE. )
        ENDIF
     ENDDO lancz_loop1
     !
     IF (charge_response == 1 ) THEN
        !
        !  Write the apropriate charge density plot.
        CALL lr_dump_rho (plot_type)
     ENDIF
     IF (project) THEN
        !
        !  Calculate projections onto virtual states if required.
        CALL lr_calc_project(ip)
     ENDIF
  ENDDO
  !
  !
  WRITE(stdout,'(5x,"End of Lanczos iterations")')
  !
  IF (project .and. n_ipol == 3) THEN
     !
     !  Final projection and report (if required).
     CALL lr_calc_project(4)
  ENDIF

  !
  !   Deallocate pw variables
  CALL clean_pw( .false. )
  WRITE(stdout,'(5x,"Finished linear response calculation...")')
  CALL stop_clock('lr_main')
  CALL print_clock_lr()
  CALL stop_lr( .TRUE. )
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<end of lr_main>")')
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Additional small-time subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
  SUBROUTINE lr_print_preamble()

    USE lr_variables, ONLY : no_hxc
    USE uspp,         ONLY : okvan
    USE funct,        only : dft_is_hybrid

    IMPLICIT NONE

!    WRITE( stdout, '(/5x,"----------------------------------------")' )
!    WRITE( stdout, '(/5x,"")' )
!    WRITE( stdout, '(/5x,"Please cite this project as:  ")' )
!    WRITE( stdout, '(/5x,"O.B. Malcioglu, R. Gebauer, D. Rocca, S. Baroni,")' )
!    WRITE( stdout, '(/5x,"""turboTDDFT – a code for the simulation of molecular")' )
!    WRITE( stdout, '(/5x,"spectra using the Liouville-Lanczos approach to")' )
!    WRITE( stdout, '(/5x,"time-dependent density-functional perturbation theory""")' )
!    WRITE( stdout, '(/5x,"CPC, 182, 1744 (2011)")' )
!    WRITE( stdout, '(/5x,"----------------------------------------")' )
    !
    IF(okvan) WRITE( stdout, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")' )
    !
    WRITE(stdout,'(/,5X,"Lanczos linear response spectrum calculation")')
    WRITE(stdout,'(5x,"Number of Lanczos iterations = ",i6)') itermax
    !
    IF (no_hxc)  THEN
       WRITE(stdout,'(5x,"No Hartree/Exchange/Correlation")')
    ELSEIF (dft_is_hybrid()) THEN
       WRITE(stdout, '(/5x,"Use of exact-exchange enabled. Note the EXX correction to the [H,X]", &
            &/5x,"commutator is NOT included hence the f-sum rule will be violated.")')
    ENDIF
    !
  END SUBROUTINE lr_print_preamble


  LOGICAL FUNCTION test_restart(test_this)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This tests whether the restart flag is applicable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    USE lr_variables,     ONLY : n_ipol,LR_polarization,restart,bgz_suffix
    USE io_files,         ONLY: prefix, tmp_dir, nd_nmbr, wfc_dir
    USE mp,               ONLY : mp_bcast, mp_barrier,mp_sum
    USE mp_world,         ONLY : world_comm
    USE io_global,        ONLY : ionode, ionode_id

    IMPLICIT NONE
    INTEGER, INTENT(in) :: test_this
    CHARACTER(len=256) :: tempfile, filename, tmp_dir_saved
    LOGICAL :: exst
    CHARACTER(len=6), EXTERNAL :: int_to_char
    INTEGER :: i, temp_restart

    !
    !test_this= 1 : d0psi files
    !test_this= 2 : lanczos restart files

    temp_restart=0
    !print *, "test_restart with restart=",restart
    IF (.not.restart) THEN
       test_restart = .false.
       RETURN
    ENDIF
    test_restart=.true.
    IF (test_this == 1) THEN
       !
       !Check for parallel i/o files that are in wfc_dir
       tmp_dir_saved = tmp_dir
       IF ( wfc_dir /= 'undefined' ) tmp_dir = wfc_dir
       !
       IF ( n_ipol == 1 ) THEN
          filename = trim(prefix)//'.d0psi.'//trim(int_to_char(LR_polarization))
          tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
          INQUIRE (file = tempfile, exist = exst)
          !print *, tempfile," exst=",exst
          IF (.not. exst) THEN
             temp_restart=1
          ENDIF
       ELSE
          DO i=1, n_ipol
             filename = trim(prefix)//'.d0psi.'//trim(int_to_char(i))
             tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
             INQUIRE (file = tempfile, exist = exst)
             !print *, tempfile," exst=",exst
             IF (.not. exst) THEN
                temp_restart=1
             ENDIF
          ENDDO
       ENDIF

       tmp_dir = tmp_dir_saved

       IF ( wfc_dir /= 'undefined' ) THEN
          ! check if these files can be read from outdir instead of wfcdir
          !
          IF ( n_ipol == 1 ) THEN
             filename = trim(prefix)//'.d0psi.'//trim(int_to_char(LR_polarization))
             tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
             INQUIRE (file = tempfile, exist = exst)
             IF (exst) THEN
                temp_restart=0
             ENDIF
          ELSE
             DO i=1, n_ipol
                filename = trim(prefix)//'.d0psi.'//trim(int_to_char(i))
                tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
                INQUIRE (file = tempfile, exist = exst)
                IF (exst) THEN
                   temp_restart=0
                ENDIF
             ENDDO
          ENDIF
       ENDIF
    ENDIF !for test_this = 1
    IF (test_this == 2) THEN

       !Restart files are always written in outdir
       IF ( n_ipol == 1 ) THEN
          filename = trim(prefix)//'.restart_lanczos.'//trim(int_to_char(LR_polarization))
          tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
       ELSE
          filename = trim(prefix)//'.restart_lanczos.'//trim(int_to_char(LR_polarization))
          tempfile = trim(tmp_dir) // trim(filename)//nd_nmbr
       ENDIF
       INQUIRE (file = tempfile, exist = exst)
       !print *, tempfile," exst=",exst
       IF (.not. exst) THEN
          temp_restart=1
       ENDIF
       !
       !End of parallel file i/o
       !
       IF ( n_ipol == 1 ) THEN
          filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(LR_polarization))
          tempfile = trim(tmp_dir) // trim(filename)
       ELSE
          filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(LR_polarization))
          tempfile = trim(tmp_dir) // trim(filename)
       ENDIF
       INQUIRE (file = tempfile, exist = exst)
       !print *, tempfile," exst=",exst
       IF (.not. exst) THEN
          temp_restart=1
       ENDIF
    ENDIF !for test_this = 2

    !print *,"temp_restart",temp_restart
#ifdef __MPI
    CALL mp_sum(temp_restart,world_comm)
#endif
    !print *, "current temp_restart", temp_restart
    IF (temp_restart > 0 ) THEN
       !print *,"restart falsified",nd_nmbr
       !WRITE(stdout,'(5X,A,3X,"is missing, unable to restart.")') offender
       WRITE(stdout,'(5X,"There are missing files!")')
       IF (test_this==1) WRITE(stdout,'(5X,"d0psi files can not be found,trying to recompansate")')
       IF (test_this==2) WRITE(stdout,'(5X,"lanczos restart files can not be found, starting run from scratch")')
       test_restart=.false.
    ENDIF

    RETURN
  END FUNCTION test_restart
END PROGRAM lr_main
!-----------------------------------------------------------------------

!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM lr_eels_main
  !---------------------------------------------------------------------
  !
  ! This is the main driver of the turboEELS code
  ! for Electron Energy Loss Spectroscopy.
  ! It applys the Lanczos algorithm to the matrix 
  ! of equations coming from TDDFPT. It can calculate:
  !
  ! Iurii Timrov, Ecole Polytechnique and SISSA, 2010-2015
  !
  USE lr_lanczos,            ONLY : one_lanczos_step
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : dp
  USE lr_variables,          ONLY : restart, restart_step, itermax, lr_verbosity,  &
                                  & evc1, evc1_old, norm0, n_ipol, &
                                  & d0psi, d0psi2, LR_iteration, LR_polarization, &
                                  & plot_type, nbnd_total, pseudo_hermitian, &
                                  & itermax_int, revc0, lr_io_level, code2, &
                                  & eels, approximation 
  USE io_files,              ONLY : nd_nmbr
  USE global_version,        ONLY : version_number
  USE ions_base,             ONLY : tau,nat,atm,ityp
  USE environment,           ONLY : environment_start
  USE mp_global,             ONLY : nimage, mp_startup, set_bgrp_indices, &
                                    ibnd_start, ibnd_end
  USE wvfct,                 ONLY : nbnd
  USE wavefunctions_module,  ONLY : psic
  USE check_stop,            ONLY : check_stop_now, check_stop_init
  USE fft_base,              ONLY : dffts
  USE uspp,                  ONLY : okvan
  USE iso_c_binding,         ONLY : c_int
  USE mp_bands,              ONLY : ntask_groups
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER             :: ip, na, pol_index, ibnd
  INTEGER             :: iter_restart, iteration
  LOGICAL             :: rflag
  INTEGER(kind=c_int) :: kilobytes
  LOGICAL, EXTERNAL   :: test_restart
  !
  pol_index = 1
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  !
  CALL environment_start ( code2 )
  !
  CALL start_clock('lr_eels_main')
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_eels_main>")')
  ENDIF
  !
  ! Let the TDDFPT routines know that they are doing EELS.
  !
  eels   = .TRUE.
  !
  ! Reading input file and PWSCF xml, some initialisation
  ! Read the input variables for TDDFPT;
  ! allocate space for all quantities already computed
  ! in the PWscf(nscf) program, and read them from the data file;
  ! define the tmp_dir directory.
  !
  CALL lr_readin ( )
  !
  CALL check_stop_init()
  !
  ! Print a preamble info about the run
  !
  CALL lr_print_preamble_eels()
  !
  ! Non-scf calculation at k and k+q
  !
  IF (.NOT.restart) THEN
     !
     WRITE( stdout, '(/,5X,"------------ Nscf calculation ---------------")')
     !
     CALL lr_run_nscf( )
     !
  ENDIF
  !
  ! Initialisation, and read the wfct's at k and k+q
  !
  CALL lr_init_nfo()
  !
  ! Output summary of the main variables of the TDDFPT code
  !
  CALL lr_summary()
  !
  ! Allocate the arrays
  !
  CALL lr_alloc_init()
  !
  ! Memory usage
  !
  CALL memstat( kilobytes )
  IF ( kilobytes > 0 ) WRITE(stdout,'(5X,"lr_eels_main, & 
      & per-process dynamical memory:",f7.1,"Mb")' ) kilobytes/1000.0
  !
  IF ( ntask_groups > 1 ) WRITE(stdout,'(5X,"Task groups is activated...")' )
  !
  ! Band groups parallelization (if activated)
  !
  CALL set_bgrp_indices(nbnd, ibnd_start, ibnd_end)
  !
  ! Set up initial response orbitals (starting Lanczos vectors)
  !
  IF ( test_restart(1) ) THEN
     CALL lr_read_d0psi()
  ELSE
     CALL lr_solve_e()
  ENDIF
  !
  DEALLOCATE( psic )
  !
  ! Calculate a derivative of the XC potential
  !
  CALL lr_dv_setup()
  !
  WRITE(stdout,'(/,5X,"LANCZOS LINEAR-RESPONSE SPECTRUM CALCULATION")')
  WRITE(stdout,'(5X," ")')
  WRITE(stdout,'(5x,"Number of Lanczos iterations = ",i6)') itermax
  !
  ! Lanczos loop where the real work happens
  !
  DO ip = 1, n_ipol
     !
     IF (n_ipol/=1) THEN
        LR_polarization = ip
        pol_index = LR_polarization
     ENDIF
     !
     ! Read the starting Lanczos vectors d0psi (EELS: and d0psi2) from the file,
     ! which was written above by lr_solve_e.
     !
     CALL lr_read_d0psi()
     !
     ! Normalization of the starting Lanczos vectors,
     ! or reading of the data from the restart file.
     !
     IF (test_restart(2)) THEN 
        !
        CALL lr_restart(iter_restart,rflag)
        !
        WRITE(stdout,'(/5x,"Restarting Lanczos loop",1x,i8)') LR_polarization
        !
     ELSE
        !
        ! The two starting Lanczos vectors are equal.
        !
        evc1(:,:,:,1) = d0psi(:,:,:,pol_index)
        !
        ! The new structure of the Lanczos algorithm
        ! does not need the normalisation of the starting Lanczos 
        ! vectors here.
        !
        evc1(:,:,:,2) = evc1(:,:,:,1)
        !
        evc1_old(:,:,:,1) = cmplx(0.0d0,0.0d0)
        evc1_old(:,:,:,2) = cmplx(0.0d0,0.0d0)
        !
        iter_restart = 1
        !
        IF (.NOT. eels) WRITE(stdout,'(/5x,"Starting Lanczos loop",1x,i8)') LR_polarization
        !
     ENDIF
     !
     ! Loop on the Lanczos iterations
     ! 
     lancz_loop1 : DO iteration = iter_restart, itermax
        !
        LR_iteration = iteration
        !
        WRITE(stdout,'(/5x,"Lanczos iteration:",1x,i6)') LR_iteration
        !
        CALL one_lanczos_step()
        !
        IF ( lr_io_level > 0 .and. (mod(LR_iteration,restart_step)==0 .or. &
                           & LR_iteration==itermax .or. LR_iteration==1) ) &
                           CALL lr_write_restart()
        !
        ! Check to see if the wall time limit has been exceeded.
        ! if it has exit gracefully saving the last set of Lanczos
        ! iterations.
        !
        IF ( check_stop_now() ) THEN
           !
           CALL lr_write_restart()
           !
           ! Deallocate PW variables.
           !
           CALL clean_pw( .FALSE. )
           CALL stop_clock('lr_main')
           CALL print_clock_lr()
           CALL stop_lr( .FALSE. )
           !
        ENDIF
        !
     ENDDO lancz_loop1
     !
  ENDDO
  ! 
  WRITE(stdout,'(5x,"End of Lanczos iterations")')
  !
  ! Deallocate PW variables
  !
  CALL clean_pw( .FALSE. )
  !
  WRITE(stdout,'(5x,"Finished linear response calculation...")')
  !
  CALL stop_clock('lr_eels_main')
  !
  CALL print_clock_lr()
  !
  CALL stop_lr( .TRUE. )
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<end of lr_eels_main>")')
  ENDIF

CONTAINS
 
SUBROUTINE lr_print_preamble_eels()
    
    USE uspp,           ONLY : okvan

    IMPLICIT NONE

    WRITE( stdout, '(/5x,"----------------------------------------")' )
    WRITE( stdout, '(/5x,"Please cite this project as:")' )
    WRITE( stdout, '(/5x,"I. Timrov, N. Vast, R. Gebauer, and S. Baroni,",                       &
                   & /5x,"Electron energy loss and inelastic x-ray scattering cross sections",   &
                   & /5x,"from time-dependent density-functional perturbation theory",           &
                   & /5x,"Phys. Rev. B 88, 064301 (2013); ibid. 91, 139901 (2015). ")' )
    WRITE( stdout, '(/5x,"and")' )
    WRITE( stdout, '(/5x,"I. Timrov, N. Vast, R. Gebauer, and S. Baroni,",                            &
                   & /5x,"turboEELS - A code for the simulation of the electron energy loss and",     &
                   & /5x,"inelastic X-ray scattering spectra using the Liouville - Lanczos approach", &
                   & /5x,"to time-dependent density-functional perturbation theory", &
                   & /5x,"Comp. Phys. Commun. 196, 460 (2015). ")' )
    WRITE( stdout, '(/5x,"----------------------------------------")' )
    !
    WRITE( stdout, '(/5x,"Using the ' // trim(approximation) // ' approximation.")' )
    !
    If (pseudo_hermitian) THEN
       WRITE( stdout, '(/5x,"Using the pseudo-Hermitian Lanczos algorithm.")' )
    ELSE
       WRITE( stdout, '(/5x,"Using the non-Hermitian Lanczos algorithm.")' )
    ENDIF
    !
    IF (okvan) WRITE( stdout, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")')
    !
    RETURN
    !
END SUBROUTINE lr_print_preamble_eels

END PROGRAM lr_eels_main
!-----------------------------------------------------------------------

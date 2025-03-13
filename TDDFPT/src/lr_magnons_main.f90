!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM lr_magnons_main
  !---------------------------------------------------------------------
  !
  ! This is the main driver of the turboMAGNONS code for spin-fluctation spectra in magnetic system.
  ! It applys the Lanczos algorithm to the matrix of equations coming from TDDFPT. 
  !
  ! Created by Tommaso Gorni (2018)
  ! Modified by Oscar Baseggio (2019)
  !
  USE lr_lanczos,            ONLY : one_lanczos_step
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : dp
  USE lr_variables,          ONLY : restart_step, itermax, lr_verbosity,  &
                                  & evc1, evc1_old, norm0, n_ipol, &
                                  & d0psi, d0psi2, LR_iteration, LR_polarization, &
                                  & plot_type, nbnd_total, pseudo_hermitian, &
                                  & itermax_int, revc0, lr_io_level, code3, &
                                  & magnons, approximation, V0psi, ipol, n_op, &
                                  & evc1_rgt, evc1_lft, evc1_rgt_old, evc1_lft_old
  USE io_files,              ONLY : nd_nmbr
  USE global_version,        ONLY : version_number
  USE ions_base,             ONLY : tau,nat,atm,ityp
  USE environment,           ONLY : environment_start
  USE mp_global,             ONLY : mp_startup
  USE mp_world,              ONLY : world_comm
  USE mp_pools,              ONLY : intra_pool_comm
  USE mp_bands,              ONLY : intra_bgrp_comm, inter_bgrp_comm, &
                                    ntask_groups
  USE mp_bands_TDDFPT,       ONLY : ibnd_start, ibnd_end
  USE command_line_options,  ONLY : ndiag_
  USE wvfct,                 ONLY : nbnd
  USE wavefunctions,         ONLY : psic
  USE check_stop,            ONLY : check_stop_now, check_stop_init
  USE fft_base,              ONLY : dffts
  USE uspp,                  ONLY : okvan
  USE clib_wrappers,         ONLY : memstat
  USE klist,                 ONLY : igk_k
  USE control_flags,         ONLY : use_gpu
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  ! Local variables
  !
  INTEGER             :: ip, na, pol_index, ibnd
  INTEGER             :: iter_restart, iteration
  LOGICAL             :: rflag
  INTEGER             :: kilobytes
  LOGICAL, EXTERNAL   :: test_restart
  LOGICAL, EXTERNAL   :: check_gpu_support
  !
  pol_index = 1
  !
  CALL mp_startup ( )
  !
  CALL environment_start ( code3 )
  !
  CALL start_clock('lr_magnons_main')
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_magnons_main>")')
  ENDIF
  !
  ! Let the TDDFPT routines know that they are doing MAGNONS.
  !
  magnons  = .TRUE.
  !
  use_gpu = check_gpu_support()
  ! Reading input file and PWSCF xml, some initialisation
  ! Read the input variables for TDDFPT;
  ! allocate space for all quantities already computed
  ! in the PWscf(nscf) program, and read them from the data file;
  ! define the tmp_dir directory.
  !
  CALL lr_readin ()
  !
  CALL check_stop_init()
  !
  ! Print a preamble info about the run
  !
  CALL lr_print_preamble_magnons()
  !
  ! NSCF calculation at k, k+q, k-q, -k, -k+q, -k-q
  !
  CALL lr_run_nscf()
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
  IF ( kilobytes > 0 ) WRITE(stdout,'(5X,"lr_magnons_main, & 
      & per-process dynamical memory:",f7.1,"Mb")' ) kilobytes/1000.0
  !
  IF ( ntask_groups > 1 ) WRITE(stdout,'(5X,"Task groups is activated...")' )
  !
  ! Band groups parallelization (if activated)
  !
  CALL divide(inter_bgrp_comm, nbnd, ibnd_start, ibnd_end)
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
     ELSE
        LR_polarization = ipol
     ENDIF
     !
     ! Read the starting Lanczos vectors V0psi and O_psi for magnons from the file,
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
        evc1_rgt(:,:,:,:) = V0psi(:,:,:,:,ip)
        evc1_lft(:,:,:,:) = evc1_rgt(:,:,:,:)
        !
        ! The new structure of the Lanczos algorithm
        ! does not need the normalisation of the starting Lanczos 
        ! vectors here.
        !
        evc1_rgt_old(:,:,:,:) = cmplx(0.0d0,0.0d0)
        evc1_lft_old(:,:,:,:) = cmplx(0.0d0,0.0d0)
        !
        iter_restart = 1
        !
        WRITE(stdout,'(/5x,"Starting Lanczos loop",1x,i8)') LR_polarization
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
           CALL stop_clock('lr_magnons_main')
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
  CALL stop_clock('lr_magnons_main')
  !
  CALL print_clock_lr()
  !
  CALL stop_lr( .TRUE. )
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<end of lr_magnons_main>")')
  ENDIF

CONTAINS
 
SUBROUTINE lr_print_preamble_magnons()
    
    USE uspp,           ONLY : okvan

    IMPLICIT NONE

    WRITE( stdout, '(/5x,"--------------------------------------------------------------------------------------------------")' )
    WRITE( stdout, '(/5x,"Please cite this project as:")' )
    WRITE( stdout, '(/5x,"T. Gorni, I. Timrov, S. Baroni, Eur. Phys. J. B 91, 249 (2018).")')
    WRITE( stdout, '(/5x,"T. Gorni, O. Baseggio, P. Delugas, S. Baroni, I. Timrov, Comput. Phys. Commun. 280, 108500 (2022).")')
    WRITE( stdout, '(/5x,"--------------------------------------------------------------------------------------------------")' )
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
END SUBROUTINE lr_print_preamble_magnons

END PROGRAM lr_magnons_main
!-----------------------------------------------------------------------

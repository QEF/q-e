!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM lr_main
  !---------------------------------------------------------------------
  !
  ! This is the main driver of the TDDFPT code
  ! for Absorption Spectroscopy. 
  ! It applys the Lanczos algorithm to the matrix 
  ! of equations coming from TDDFPT.
  !
  ! Brent Walker, ICTP, 2004
  ! Dario Rocca, SISSA, 2006
  ! Osman Baris Malcioglu, SISSA, 2008
  ! Simone Binnie, SISSA, 2011
  ! Xiaochuan Ge, SISSA, 2013
  ! Iurii Timrov, SISSA, 2015
  !
  USE lr_lanczos,            ONLY : one_lanczos_step
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : dp
  USE lr_variables,          ONLY : restart, restart_step, itermax, lr_verbosity,  &
                                  & evc1, evc1_old,norm0, charge_response, n_ipol, &
                                  & d0psi, LR_iteration, LR_polarization, &
                                  & plot_type, no_hxc, nbnd_total, project, F,R, &
                                  & itermax_int, revc0, lr_io_level, code1
  USE io_files,              ONLY : nd_nmbr
  USE global_version,        ONLY : version_number
  USE charg_resp,            ONLY : lr_calc_w_T, read_wT_beta_gamma_z, lr_project_init,&
                                  & lr_dump_rho_tot_compat1, lr_dump_rho_tot_cube,&
                                  & lr_dump_rho_tot_xyzd,lr_dump_rho_tot_xcrys,&
                                  & lr_dump_rho_tot_pxyd,chi,lr_calc_R,w_t_norm0_store,&
                                  & resonance_condition, lr_dump_rho, lr_calc_project
  USE ions_base,             ONLY : tau,nat,atm,ityp
  USE environment,           ONLY : environment_start
  USE mp_global,             ONLY : nimage, mp_startup, set_bgrp_indices, &
                                    ibnd_start, ibnd_end
  USE wvfct,                 ONLY : nbnd
  USE wavefunctions_module,  ONLY : psic
  USE control_flags,         ONLY : tddfpt
  USE check_stop,            ONLY : check_stop_now, check_stop_init
  USE funct,                 ONLY : dft_is_hybrid
  USE fft_base,              ONLY : dffts
  USE uspp,                  ONLY : okvan
  USE mp_bands,              ONLY : ntask_groups
#if defined(__ENVIRON)
  USE plugin_flags,          ONLY : use_environ
  USE environ_info,          ONLY : environ_summary
#endif
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER            :: ip,na,pol_index,ibnd
  INTEGER            :: iter_restart,iteration
  LOGICAL            :: rflag
  COMPLEX(kind=dp)   :: temp
  LOGICAL, EXTERNAL  :: test_restart
  !
  pol_index = 1
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  !
  CALL environment_start ( code1 )
  !
  CALL start_clock('lr_main')
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_main>")')
  ENDIF
  !
  ! Let the routines of the Environ plugin know that 
  ! they are doing TDDFPT.
  !
  tddfpt = .TRUE.
  !
  ! Reading input file and PWSCF xml, some initialisation;
  ! Read the input variables for TDDFPT;
  ! Allocate space for all quantities already computed
  ! by PWscf, and read them from the data file;
  ! Define the tmp_dir directory.
  !
  CALL lr_readin ( )
  !
  ! Writing a summary to the standard output 
  ! about Environ variables
  !
#if defined(__ENVIRON)
  IF ( use_environ ) CALL environ_summary()
#endif
  !
  CALL check_stop_init()
  !
  ! Initialisation
  !
  CALL lr_init_nfo()
  !
  ! Allocate the arrays
  !
  CALL lr_alloc_init()
  !
  ! Print a preamble info about the run
  !
  CALL lr_print_preamble()
  !
  IF ( ntask_groups > 1 ) WRITE(stdout,'(5X,"Task groups is activated...")' )
  !
  ! Read ground state wavefunctions from PWscf.
  !
  CALL lr_read_wf()
  !
  ! Band groups parallelization (if activated)
  !
  CALL set_bgrp_indices(nbnd,ibnd_start,ibnd_end)
  !
  ! Set up initial response orbitals (starting Lanczos vectors)
  !
  IF ( test_restart(1) ) THEN
     CALL lr_read_d0psi()
  ELSE
     CALL lr_solve_e()
  ENDIF
  !
  !do ip = 1, n_ipol
  !  temp = wfc_dot(ibnd)
  !enddo
  !
  DEALLOCATE( psic )
  !
  ! Projection analysis
  !
  IF (project) CALL lr_project_init()
  !
  ! Calculate a derivative of the XC potential
  !
  CALL lr_dv_setup()
  !
  ! S. Binnie: Write coordinates of the read atom (just in case).
  !
  IF (lr_verbosity > 1) THEN
     !
     WRITE(stdout,'(/,5X,"Positions of atoms in internal coordinates")')
     !
     DO na = 1, nat 
        WRITE(stdout,'(5X,A3,2X,3F15.8)') atm(ityp(na)), tau(1:3,na)
     ENDDO
     !
  ENDIF
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
     ! This is for the response of the charge density
     !
     IF (charge_response == 1) THEN
        !
        ! Read precalculated beta, gamma, and z.
        !
        CALL read_wT_beta_gamma_z()
        CALL lr_calc_w_T()
        !
     ENDIF
     !
     ! Read the starting Lanczos vectors d0psi from the file which
     ! was written above by lr_solve_e.
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
        ! Xiaochuan Ge: The new structure of the Lanczos algorithm
        ! does not need the normalisation of the starting Lanczos 
        ! vectors here.
        !
        !CALL lr_normalise( evc1(:,:,:,1), norm0(pol_index) ) 
        !
        evc1(:,:,:,2) = evc1(:,:,:,1)
        !
        evc1_old(:,:,:,1) = cmplx(0.0d0,0.0d0)
        evc1_old(:,:,:,2) = cmplx(0.0d0,0.0d0)
        !
        iter_restart = 1
        !
        WRITE(stdout,'(/5x,"Starting Lanczos loop",1x,i8)') LR_polarization
        !
     ENDIF
     !
     ! d0psi = S * d0psi 
     ! This is needed in lr_lanczos for the dot product
     ! in the calculation of the zeta-coefficients.
     !
     IF (okvan) CALL sd0psi() 
     !
     ! Loop on the Lanczos iterations
     ! 
     lancz_loop1 : DO iteration = iter_restart, itermax
        !
        LR_iteration = iteration
        !
        WRITE(stdout,'(/5x,"Lanczos iteration:",1x,i6,3x,"Pol:",i1,i8)') LR_iteration, ip
        !
        CALL one_lanczos_step()
        !
        IF ( lr_io_level > 0 .and. (mod(LR_iteration,restart_step)==0 .or. &
                           & LR_iteration==itermax .or. LR_iteration==1) ) &
                           CALL lr_write_restart()
        !
        ! Check to see if the wall time limit has been exceeded.
        ! If it has exit gracefully saving the last set of Lanczos
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
     IF (charge_response == 1) THEN
        !
        ! Write the apropriate charge density plot.
        !
        CALL lr_dump_rho (plot_type)
        !
     ENDIF
     !
     IF (project) THEN
        !
        ! Calculate projections onto virtual states if required.
        !
        CALL lr_calc_project(ip)
        !
     ENDIF
     !
  ENDDO
  ! 
  WRITE(stdout,'(5x,"End of Lanczos iterations")')
  !
  IF (project .AND. n_ipol == 3) THEN
     !
     !  Final projection and report (if required).
     !
     CALL lr_calc_project(4) 
     !
  ENDIF
  !
  ! Deallocate PW variables
  !
  CALL clean_pw( .FALSE. )
  !
  WRITE(stdout,'(5x,"Finished linear response calculation...")')
  !
  CALL stop_clock('lr_main')
  !
  CALL print_clock_lr()
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<end of lr_main>")')
  ENDIF
  !
  CALL stop_lr( .TRUE. )

CONTAINS
 
SUBROUTINE lr_print_preamble()
    
    USE lr_variables,        ONLY : no_hxc, d0psi_rs
    USE uspp,                ONLY : okvan
    USE funct,               ONLY : dft_is_hybrid
    USE martyna_tuckerman,   ONLY : do_comp_mt
    USE control_flags,       ONLY : do_makov_payne
#if defined(__ENVIRON)
    USE plugin_flags,        ONLY : use_environ
#endif    

    IMPLICIT NONE
    !
    WRITE( stdout, '(/5x,"=-----------------------------------------------------------------=")')
    WRITE( stdout, '(/5x,"Please cite the TDDFPT project as:")')
    WRITE( stdout, '(7x,"O. B. Malcioglu, R. Gebauer, D. Rocca, and S. Baroni,")')
    WRITE( stdout, '(7x,"Comput. Phys. Commun. 182, 1744 (2011)")')
    WRITE( stdout, '(5x,"and")' )
    WRITE( stdout, '(7x,"X. Ge, S. J. Binnie, D. Rocca, R. Gebauer, and S. Baroni,")')
    WRITE( stdout, '(7x,"Comput. Phys. Commun. 185, 2080 (2014)")')
#if defined(__ENVIRON)
    IF ( use_environ ) THEN
      WRITE( stdout, '(5x,"and the TDDFPT+Environ project as:")' )
      WRITE( stdout, '(7x,"I. Timrov, O. Andreussi, A. Biancardi, N. Marzari, and S. Baroni,")' )
      WRITE( stdout, '(7x,"J. Chem. Phys. 142, 034111 (2015)")' )
    ENDIF
#endif
    WRITE( stdout, '(5x,"in publications and presentations arising from this work.")' )
    WRITE( stdout, '(/5x,"=-----------------------------------------------------------------=")')
    !
    IF (okvan) WRITE( stdout, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")' )
    !
    IF (do_comp_mt) THEN
       WRITE( stdout, '(/5x,"Martyna-Tuckerman periodic-boundary correction is used")' )
    ELSEIF (do_makov_payne) THEN
       WRITE( stdout, '(/5x,"WARNING! Makov-Payne periodic-boundary correction was activated in PWscf,",  &
                      & /5x,"but it is of no use for TDDFPT. It just corrects the total energy in PWscf", &
                      & /5x,"(post-processing correction) and nothing more, thus no effect on TDDFPT.", &
                      & /5x,"You can try to use the Martyna-Tuckerman correction scheme instead.")' )
    ENDIF
    !
    IF (no_hxc)  THEN
       WRITE(stdout,'(5x,"No Hartree/Exchange/Correlation")')
    ELSEIF (dft_is_hybrid() .AND. .NOT.d0psi_rs) THEN
       WRITE(stdout, '(/5x,"Use of exact-exchange enabled. Note the EXX correction to the [H,X]", &
                     & /5x,"commutator is NOT included hence the f-sum rule will be violated.",   &
                     & /5x,"You can try to use the variable d0psi_rs=.true. (see the documentation).")' )
    ENDIF
    !
    RETURN
    !
END SUBROUTINE lr_print_preamble

END PROGRAM lr_main
!-----------------------------------------------------------------------

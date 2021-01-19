!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM lr_eels_main
  !---------------------------------------------------------------------
  !
  ! This is the main driver of the turboEELS code for Electron Energy Loss Spectroscopy.
  ! It applys the Lanczos or Sternheimer algorithm to the matrix of 
  ! equations coming from TDDFPT. 
  !
  ! Iurii Timrov (Ecole Polytechnique, SISSA, and EPFL) 2010-2018
  ! Oscar Baseggio (SISSA) 2020
  !
  USE lr_lanczos,            ONLY : one_lanczos_step
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : dp
  USE lr_variables,          ONLY : restart_step, itermax, lr_verbosity,  &
                                  & evc1, evc1_old, norm0, n_ipol, &
                                  & d0psi, d0psi2, LR_iteration, LR_polarization, &
                                  & plot_type, nbnd_total, pseudo_hermitian, &
                                  & itermax_int, revc0, lr_io_level, code2, &
                                  & eels, approximation, calculator, fru, fiu, &
                                  & current_w, nfs, start_freq, last_freq, chirr, &
                                  & chirz, chizz, chizr, epsm1 
  USE ions_base,             ONLY : tau,nat,atm,ityp
  USE environment,           ONLY : environment_start
  USE mp_global,             ONLY : mp_startup
  USE mp_bands,              ONLY : inter_bgrp_comm, ntask_groups
  USE mp_bands_TDDFPT,       ONLY : ibnd_start, ibnd_end
  USE wvfct,                 ONLY : nbnd
  USE wavefunctions,         ONLY : psic
  USE check_stop,            ONLY : check_stop_now, check_stop_init
  USE fft_base,              ONLY : dffts
  USE uspp,                  ONLY : okvan
  USE wrappers,              ONLY : memstat
  USE lr_sternheimer,        ONLY : one_sternheimer_step
  USE control_lr,            ONLY : flmixdpot
  USE control_flags,         ONLY : use_para_diag
  USE qpoint,                ONLY : xq
  USE uspp_param,            ONLY : nhm
  USE noncollin_module,      ONLY : noncolin
  USE lsda_mod,              ONLY : nspin
  USE lrus,                  ONLY : intq, intq_nc
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER             :: ip, na, pol_index, ibnd, iu
  INTEGER             :: iter_restart, iteration
  LOGICAL             :: rflag
  INTEGER             :: kilobytes
  LOGICAL, EXTERNAL   :: test_restart
  !
  pol_index = 1
  !
  CALL mp_startup ( )
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
  ! NSCF calculation at k and k+q
  !
  CALL lr_run_nscf( )
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
  CALL divide(inter_bgrp_comm,nbnd, ibnd_start, ibnd_end)
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
  IF (trim(calculator)=='sternheimer') THEN
     !
     ! Sternheimer algorithm
     !
     WRITE(stdout,'(/,5X,"STERNHEIMER LINEAR-RESPONSE SPECTRUM CALCULATION")')
     !
     WRITE(stdout,'(/5x,"Number of frequencies = ",i6)') nfs
     !
     IF (okvan) THEN
        ALLOCATE (intq (nhm, nhm, nat) )
        IF (noncolin) ALLOCATE(intq_nc( nhm, nhm, nat, nspin))
        call lr_compute_intq()
     ENDIF
     !
     ! Set flmixdpot
     !
     flmixdpot = 'mixd'
     !
     ALLOCATE(chirr(nfs))
     ALLOCATE(chirz(nfs))
     ALLOCATE(chizr(nfs))
     ALLOCATE(chizz(nfs))
     ALLOCATE(epsm1(nfs))
     chirr=(0.0_DP,0.0_DP)
     chizr=(0.0_DP,0.0_DP)
     chirz=(0.0_DP,0.0_DP)
     chizz=(0.0_DP,0.0_DP)
     epsm1=(0.0_DP,0.0_DP)
     !
     ! frequencies loop
     !
     DO iu = start_freq, last_freq
        !
        current_w=CMPLX(fru(iu), fiu(iu))
        !
        WRITE( stdout, '(/,5x,70("-"))')
        WRITE( stdout, '(/,20x," Step ",i5,"  /",i5)') iu, nfs
        WRITE( stdout, '(/,10x," &
              &  Frequency = (",f8.4,", ",f8.4,") Ry")') current_w
        !
        CALL one_sternheimer_step( iu, 1 )
        ! 
        CALL write_chi_on_disk(iu)
        !  
     ENDDO
     !
     DEALLOCATE(chirr)
     DEALLOCATE(chirz)
     DEALLOCATE(chizr)
     DEALLOCATE(chizz)
     DEALLOCATE(epsm1)
     !
     IF (okvan) THEN
        DEALLOCATE (intq )
        IF (noncolin) DEALLOCATE(intq_nc)
     ENDIF
     !
  ELSEIF (trim(calculator)=='lanczos') THEN
     !
     ! Lanczos algorithm
     !
     WRITE(stdout,'(/,5X,"LANCZOS LINEAR-RESPONSE SPECTRUM CALCULATION")')

     IF (pseudo_hermitian) THEN
        WRITE( stdout, '(/5x,"Using the pseudo-Hermitian Lanczos algorithm")' )
     ELSE
        WRITE( stdout, '(/5x,"Using the non-Hermitian Lanczos algorithm")' )
     ENDIF

     WRITE(stdout,'(/5x,"Number of Lanczos iterations = ",i6)') itermax
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
           evc1_old(:,:,:,1) = (0.0d0,0.0d0)
           evc1_old(:,:,:,2) = (0.0d0,0.0d0)
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
  ELSE
     !
     CALL errore('lr_eels_main', 'Not known type of the calculator',1)
     !
  ENDIF
  !
  ! Deallocate PW variables
  !
  CALL clean_pw( .FALSE. )
  !
  WRITE(stdout,'(/5x,"Finished linear response calculation...")')
  !
  CALL stop_clock('lr_eels_main')
  !
  CALL print_clock_lr()
  !
  IF ( use_para_diag ) CALL laxlib_end()
  CALL stop_lr( .TRUE. )
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<end of lr_eels_main>")')
  ENDIF

CONTAINS
 
SUBROUTINE lr_print_preamble_eels()
    
    USE uspp,           ONLY : okvan

    IMPLICIT NONE
    
    WRITE( stdout, '(/5x,74("-"))')
    WRITE( stdout, '(/5x,"Please cite this project as:")' )
    WRITE( stdout, '(/5x,"I. Timrov, N. Vast, R. Gebauer, and S. Baroni,",                       &
                   & /5x,"Electron energy loss and inelastic x-ray scattering cross sections",   &
                   & /5x,"from time-dependent density-functional perturbation theory",           &
                   & /5x,"Phys. Rev. B 88, 064301 (2013); ibid. 91, 139901 (2015). ")' )
    WRITE( stdout, '(/5x,"I. Timrov, N. Vast, R. Gebauer, and S. Baroni,",                            &
                   & /5x,"turboEELS - A code for the simulation of the electron energy loss and",     &
                   & /5x,"inelastic X-ray scattering spectra using the Liouville - Lanczos approach", &
                   & /5x,"to time-dependent density-functional perturbation theory", &
                   & /5x,"Comp. Phys. Commun. 196, 460 (2015). ")' )
    WRITE( stdout, '(/5x,"O. Motornyi, N. Vast, I. Timrov, O. Baseggio, S. Baroni, and A. Dal Corso", &
                   & /5x,"Electron energy loss spectroscopy of bulk gold with ultrasoft",             &
                   & /5x,"pseudopotentials and the Liouville-Lanczos method",                         &
                   & /5x,"Phys. Rev. B  102, 035156 (2020). ")' )
    WRITE( stdout, '(/5x,74("-"))')
    !
    WRITE( stdout, '(/5x,"Using the ' // trim(approximation) // ' approximation")' )
    !
    IF (okvan) WRITE( stdout, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")')
    !
    WRITE( stdout, '(/5x,"Transferred momentum: q = (",3f10.5,"  ) in units 2*pi/a")') xq(1:3)
    !
    RETURN
    !
END SUBROUTINE lr_print_preamble_eels

SUBROUTINE write_chi_on_disk(iu)

    USE kinds,            ONLY : DP
    USE lr_variables,     ONLY : current_w, fru, fiu, &
                                 chirr, chirz, chizz, &
                                 chizr, epsm1, units
    USE lsda_mod,         ONLY : nspin, lsda
    USE mp_images,        ONLY : my_image_id
    USE io_global,        ONLY : stdout, ionode
    USE io_files,         ONLY : prefix
    USE noncollin_module, ONLY : noncolin, nspin_mag
    USE constants,        ONLY : rytoev

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iu
    INTEGER :: iu_epsil
    COMPLEX(DP) :: epsi
    CHARACTER(LEN=6) :: int_to_char
    CHARACTER(LEN=256) :: filename

    LOGICAL :: exst

    IF (ionode) THEN
       IF (nspin_mag==1) THEN
          !
          !   nonmagnetic case or noncollinear with time reversal
          !
          iu_epsil=2
          filename = trim(prefix) // '.plot_chi.dat'
          IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
          INQUIRE(FILE=TRIM(filename), exist=exst)
          IF (exst.AND.iu>1) THEN
             OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', &
                                  POSITION='append', FORM='formatted')
          ELSE
             OPEN (UNIT=iu_epsil, FILE=TRIM(filename), &
                                    STATUS='unknown', FORM='formatted')
             IF (units == 0) THEN
                WRITE(iu_epsil,'("# \hbar \omega(Ry)       Re(chi)         Im(chi) ")')
             ELSEIF (units == 1) THEN
                WRITE(iu_epsil,'("# \hbar \omega(eV)       Re(chi)         Im(chi) ")')
             ENDIF
          END IF
          !
          IF (units == 0) THEN
             WRITE(iu_epsil,'(4x,f12.8,4x,e15.7,2x,e15.7)') fru(iu), &
                  DREAL(chirr(iu)), DIMAG(chirr(iu))
          ELSEIF (units == 1) THEN
             WRITE(iu_epsil,'(4x,f12.8,4x,e15.7,2x,e15.7)') fru(iu)*rytoev, &
                  DREAL(chirr(iu)), DIMAG(chirr(iu))
          ENDIF
          CLOSE(iu_epsil)
       ELSEIF (nspin_mag==2) THEN
          !
          !  lsda case
          !
          iu_epsil=2
          filename = trim(prefix) // '.plot_re_chi_mag.dat'
          IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
          INQUIRE(FILE=TRIM(filename), exist=exst)
          IF (exst.AND.iu>1) THEN
             OPEN (UNIT=iu_epsil, FILE=TRIM(filename), &
                             STATUS='old', POSITION='append', FORM='formatted')
          ELSE
             OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='unknown', &
                                                FORM='formatted')
             WRITE(iu_epsil,'("#  \hbar \omega(Ry)        rr       rz        zz&
                              &         +-        -+       xx        xy")')
          ENDIF
          !
          WRITE(iu_epsil,'(4x,f12.8,4x,3e15.7)') fru(iu), &
                     DREAL(chirr(iu)), DREAL(chirz(iu)),  &
                     DREAL(chizz(iu))
          !
          CLOSE(iu_epsil)
          !
          iu_epsil=2
          filename = trim(prefix) // '.plot_im_chi_mag.dat'
          IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
          INQUIRE(FILE=TRIM(filename), exist=exst)
          IF (exst.AND.iu>1) THEN
             OPEN (UNIT=iu_epsil, FILE=TRIM(filename), &
                           STATUS='old', POSITION='append', FORM='formatted')
          ELSE
             OPEN (UNIT=iu_epsil, FILE=TRIM(filename), &
                               STATUS='unknown', FORM='formatted')
             WRITE(iu_epsil,'("#  \hbar \omega(Ry)     rr       rz        zz &
                              &         +-        -+       xx        xy")')
          ENDIF
          !
          WRITE(iu_epsil,'(4x,f12.8,4x,3e15.7)') fru(iu),   &
                        DIMAG(chirr(iu)), DIMAG(chirz(iu)), &
                        DIMAG(chizz(iu))
          !
          CLOSE(iu_epsil)
       ELSE
          !
          ! Noncollinear case: not yet implemented
          CALL errore('write_chi_on_disk', 'Noncollinear case not implemented yet',1)
          !
       ENDIF
       !
       IF (.NOT.lsda) THEN
          iu_epsil=2
          filename = trim(prefix) // '.plot_eps.dat'
          IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
          INQUIRE(FILE=TRIM(filename), exist=exst)
          IF (exst.AND.iu>1) THEN
             OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', &
                                  POSITION='append', FORM='formatted')
          ELSE
             OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='unknown', &
                                                            FORM='formatted')
             IF (units == 0) THEN
                WRITE(iu_epsil,'("#  \hbar \omega(Ry)     Re(1/eps)     -Im(1/eps)&
                              &      Re(eps)      Im(eps) ")')
             ELSEIF (units == 1) THEN
                WRITE(iu_epsil,'("#  \hbar \omega(eV)     Re(1/eps)     -Im(1/eps)&
                              &      Re(eps)      Im(eps) ")')
             ENDIF
          ENDIF
          !
          epsi = (0.0_DP,0.0_DP)
          IF (ABS(epsm1(iu))>1.D-10) epsi = CMPLX(1.0_DP,0.0_DP) / epsm1(iu)
          IF (units == 0 ) THEN
             WRITE(iu_epsil,'(4x,f12.8,4x,4e14.6)') fru(iu),  &
                    DREAL(epsm1(iu)), -DIMAG(epsm1(iu)), &
                    DREAL(epsi), DIMAG(epsi)
          ELSEIF (units == 1) THEN
             WRITE(iu_epsil,'(4x,f12.8,4x,4e14.6)') fru(iu)*rytoev,  &
                    DREAL(epsm1(iu)), -DIMAG(epsm1(iu)), &
                    DREAL(epsi), DIMAG(epsi)
          ENDIF
          !
          CLOSE(iu_epsil)
       ENDIF
       !
    ENDIF
    !
    RETURN
END SUBROUTINE write_chi_on_disk

END PROGRAM lr_eels_main
!-----------------------------------------------------------------------

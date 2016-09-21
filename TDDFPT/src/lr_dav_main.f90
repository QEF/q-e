!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------

PROGRAM lr_dav_main
  !---------------------------------------------------------------------
  ! Xiaochuan Ge, SISSA, 2013
  !---------------------------------------------------------------------
  ! ... overall driver routine for applying davidson algorithm
  ! ... to the matrix of equations coming from tddft
  !---------------------------------------------------------------------

  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : dp
  USE lr_variables,          ONLY : restart, restart_step,&
       evc1,n_ipol, d0psi, &
       no_hxc, nbnd_total, &
       revc0, lr_io_level, code1,davidson
  USE io_files,              ONLY : nd_nmbr
  USE global_version,        ONLY : version_number
  USE ions_base,             ONLY : tau,nat,atm,ityp
  USE environment,           ONLY : environment_start
  USE mp_global,             ONLY : nimage, mp_startup, set_bgrp_indices, &
                                    ibnd_start, ibnd_end
  USE wvfct,                 ONLY : nbnd
  USE wavefunctions_module,  ONLY : psic
  USE control_flags,         ONLY : tddfpt, do_makov_payne
  USE check_stop,            ONLY : check_stop_now, check_stop_init
  USE funct,                 ONLY : dft_is_hybrid

  use lr_dav_routines
  use lr_dav_variables
  use lr_dav_debug
#if defined(__ENVIRON)
  USE plugin_flags,          ONLY : use_environ
  USE environ_info,          ONLY : environ_summary
#endif
  !
  IMPLICIT NONE
  INTEGER            :: ibnd_occ,ibnd_virt,ibnd,ip
  LOGICAL            :: rflag, nomsg
  complex(dp)            :: temp

#if defined(__MPI)
  CALL mp_startup ( )
#endif
  !
  ! Let the routines of the Environ plugin know that 
  ! they are doing TDDFPT. 
  !
  tddfpt = .true.
  !
  ! Tell to the code that we are using the Davidson method
  !
  davidson = .true.
  !
  CALL environment_start ( code1 )
  CALL start_clock('lr_dav_main')

  !   Reading input file and PWSCF xml, some initialisation
  CALL lr_readin ( )

  ! Writing a summary to the standard output about Environ variables
#if defined(__ENVIRON)
  IF ( use_environ ) CALL environ_summary()
#endif
 
  CALL check_stop_init()

  CALL lr_init_nfo() !Initialisation of degauss/openshell related stuff

  n_ipol = 3 ! Davidson automaticly calculates all three polarizations
  CALL lr_alloc_init()   ! Allocate and zero lr variables

  !   Now print some preamble info about the run to stdout
  CALL lr_print_preamble()

  !   Read in ground state wavefunctions
  CALL lr_read_wf()
  !
  CALL set_bgrp_indices(nbnd,ibnd_start,ibnd_end)

  !   Set up initial response orbitals
  CALL lr_solve_e()
  DEALLOCATE( psic )

  if( if_dft_spectrum) call dft_spectrum() ! If we just want to calculate the dft_spectrum
                                           ! it is actually not necessary to do Davidson
                                           ! iteration
  call lr_dav_alloc_init() ! allocate for davidson algorithm
  CALL lr_dav_set_init()
  
  !   Set up initial stuff for derivatives
  CALL lr_dv_setup()

  !   Davidson loop
  if (precondition) write(stdout,'(/5x,"Precondition is used in the algorithm,")')
  do while (.not. dav_conv .and. dav_iter .lt. max_iter)
    dav_iter=dav_iter+1
      if(if_check_orth) call check_orth()
      ! In one david step, M_C,M_D and M_CD are first constructed;then will be
        ! solved rigorously; then the solution in the subspace left_sub() will
        ! be transformed into full space left_full()
      call one_dav_step()
      call dav_calc_residue()
      call dav_expan_basis()
      ! 
      ! Check to see if the wall time limit has been exceeded.
      if ( check_stop_now() ) then
         call lr_write_restart_dav() 
         goto 100
      endif
      !
  enddo
  ! call check_hermitian()
  ! Extract physical meaning from the solution
  
  call interpret_eign('END')
  ! The check_orth at the end may take quite a lot of time in the case of 
  ! USPP because we didn't store the S* vector basis. Turn this step on only
  ! in cases of debugging
  ! call check_orth() 
  if(lplot_drho) call plot_drho()

100 continue
  !   Deallocate pw variables
  CALL clean_pw( .false. )
  WRITE(stdout,'(5x,"Finished linear response calculation...")')
  CALL stop_clock('lr_dav_main')
  CALL print_clock_lr()
  CALL stop_lr( .false. )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Additional small-time subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
  SUBROUTINE lr_print_preamble()

    USE lr_variables,        ONLY : no_hxc, d0psi_rs
    USE uspp,                ONLY : okvan
    USE funct,               ONLY : dft_is_hybrid
    USE martyna_tuckerman,   ONLY : do_comp_mt
#if defined(__ENVIRON)
    USE plugin_flags,        ONLY : use_environ
#endif

    IMPLICIT NONE
    !
    WRITE( stdout, '(/5x,"=-----------------------------------------------------------------=")')
    WRITE( stdout, '(/5x,"Please cite the TDDFPT project as:")')
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
    IF(okvan) WRITE( stdout, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")' )
    
    IF (do_comp_mt) THEN
       WRITE( stdout, '(/5x,"Martyna-Tuckerman periodic-boundary correction is used")' )
    ELSEIF (do_makov_payne) THEN
       WRITE( stdout, '(/5x,"WARNING! Makov-Payne periodic-boundary correction was activated in PWscf,",  &
                      & /5x,"but it is of no use for TDDFPT. It just corrects the total energy in PWscf", &
                      & /5x,"(post-processing correction) and nothing more, thus no effect on TDDFPT.", &
                      & /5x,"You can try to use the Martyna-Tuckerman correction scheme instead.")' )
    ENDIF
 
    IF (no_hxc)  THEN
       WRITE(stdout,'(5x,"No Hartree/Exchange/Correlation")')
    ELSEIF (dft_is_hybrid() .AND. .NOT.d0psi_rs) THEN
       WRITE(stdout, '(/5x,"Use of exact-exchange enabled. Note the EXX correction to the [H,X]", &
                     & /5x,"commutator is NOT included hence the f-sum rule will be violated.",   &
                     & /5x,"You can try to use the variable d0psi_rs=.true. (see the documentation).")' )
    ENDIF

  END SUBROUTINE lr_print_preamble

END PROGRAM lr_dav_main
!-----------------------------------------------------------------------

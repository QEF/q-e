!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM pwscf
  !----------------------------------------------------------------------------
  !
  ! ... Plane Wave Self-Consistent Field code 
  !
  USE io_global,        ONLY : stdout, ionode
  USE parameters,       ONLY : ntypx, npk, lmaxx, nchix, ndmx, nqfx, nbrx
  USE global_version,   ONLY : version_number
  USE wvfct,            ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin
  USE control_flags,    ONLY : nstep, istep, conv_elec, conv_ions, &
                               lpath, lmetadyn
  USE io_files,         ONLY : nd_nmbr
  USE ions_base,        ONLY : tau
  USE path_variables,   ONLY : conv_path
  USE check_stop,       ONLY : check_stop_init
  USE path_base,        ONLY : initialize_path, search_mep
  USE metadyn_base,     ONLY : metadyn_init
  USE path_io_routines, ONLY : io_path_start, io_path_stop
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  CHARACTER (LEN=9) :: code = 'PWSCF'
  !
  !
  ! ... use ".FALSE." to disable all clocks except the total cpu time clock
  ! ... use ".TRUE."  to enable clocks
  !
  CALL init_clocks( .TRUE. )
  CALL start_clock( code )
  !
  CALL startup( nd_nmbr, code, version_number )
  !
#if defined (EXX)
  WRITE( UNIT = stdout, &
       & FMT = '(/,5X,"!!! EXPERIMENTAL VERSION WITH EXX STUFF  !!!", &
       &         /,5X,"!!! DO NOT USE IT FOR ANY PRODUCTION RUN !!!")' )

#endif
  !
  IF ( ionode ) THEN
     !
     WRITE( UNIT = stdout, &
            FMT = '(/5X,"Ultrasoft (Vanderbilt) Pseudopotentials")')
     !
     WRITE( unit = stdout, FMT = 9010 ) &
         ntypx, npk, lmaxx, nchix, ndmx, nbrx, nqfx
     !
  END IF   
  !
  CALL iosys()
  !
  CALL check_stop_init()
  !
  IF ( ionode .AND. noncolin ) &
    WRITE( UNIT = stdout, &
         & FMT = '(/,5X,"non-colinear magnetization allowed",/)' )
  IF ( ionode .AND. gamma_only ) &
    WRITE( UNIT = stdout, &
         & FMT = '(/,5X,"gamma-point specific algorithms are used",/)' )
  !
  IF ( lpath ) THEN
     !
     CALL io_path_start()
     !
     CALL initialize_path( 'PW' )
     !
     ! ... this routine does all the "string" job
     !   
     CALL search_mep()
     !
     CALL io_path_stop()
     !  
     CALL stop_run( conv_path )
     !
  ELSE
     !
     istep = 0
     !
     IF ( lmetadyn ) THEN
        !
        ! ... meta-dynamics
        !
        CALL metadyn_init( 'PW', tau )
        !
        CALL init_run()
        !
        CALL metadyn()
        !
     ELSE
        !
        CALL init_run()
        !
        main_loop: DO
           !
           istep = istep + 1
           !
           ! ... electronic self-consistentcy
           !
           CALL electrons()
           !
           IF ( .NOT. conv_elec ) CALL stop_run( conv_elec )
           !
           ! ... if requested ions are moved
           !
           CALL ions()
           !
           ! ... exit condition (ionic convergence) is checked here
           !
           IF ( conv_ions .OR. ( istep >= nstep ) ) EXIT main_loop
           !
           ! ... the ionic part of the hamiltonian is reinitialized
           !
           CALL hinit1()
           !
        END DO main_loop
        !
     END IF
     !
     CALL stop_run( conv_ions )
     !
  END IF      
  !
  STOP
  !
9010 FORMAT( /,5X,'Current dimensions of program pwscf are:', /, &
           & /,5X,'ntypx = ',I2,'   npk = ',I5,'  lmax = ',I2   &
           & /,5X,'nchix = ',I2,'  ndmx = ',I5,'  nbrx = ',I2,'  nqfx = ',I2 )
  !
END PROGRAM pwscf

!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
PROGRAM pwscf
  !----------------------------------------------------------------------------
  !
  ! ... Plane Wave Self-Consistent Field code 
  !
  USE io_global,        ONLY : stdout
  USE parameters,       ONLY : ntypx, npk, lmaxx, nchix, ndm, nqfm, nbrx
  USE global_version,   ONLY : version_number
  USE wvfct,            ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin
  USE varie,            ONLY : nstep, istep, conv_elec, conv_ions, lneb
  USE io_files,         ONLY : nd_nmbr, iunneb
  USE neb_variables,    ONLY : conv_neb
  USE neb_variables,    ONLY : neb_deallocation
  USE input_parameters, ONLY : deallocate_input_parameters
  USE neb_routines,     ONLY : initialize_neb, search_mep
  USE io_routines,      ONLY : write_output
#ifdef __PARA
  USE para,             ONLY : me, mypool
#endif
  !
  IMPLICIT NONE
  !  
  CHARACTER (LEN=9)  :: code = 'PWSCF'
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
  WRITE( UNIT = stdout, &
         FMT = '(/5X,"Ultrasoft (Vanderbilt) Pseudopotentials")')
  !
  WRITE( unit = stdout, FMT = 9010 ) ntypx, npk, lmaxx, nchix, ndm, nbrx, nqfm
  !
  CALL iosys()
  !
  IF ( noncolin ) &
    WRITE( UNIT = stdout, &
         & FMT = '(/,5X,"non-colinear magnetization allowed",/)' )
  IF ( gamma_only ) &
    WRITE( UNIT = stdout, &
         & FMT = '(/,5X,"gamma-point specific algorithms are used",/)' )
  !
  CALL show_memory()
  !
  IF ( lneb ) THEN
     !
     ! ... stdout is connected to a file ( specific for each image ) 
     ! ... via unit 17
     !
#ifdef __PARA
     IF ( me == 1 .AND. mypool == 1 ) THEN
#endif
        !
        stdout = 17
        !
#ifdef __PARA
     END IF
#endif     
     !
     CALL initialize_neb()
     !
     ! ... this routine does all the NEB job
     !   
     CALL search_mep()
     !
     ! ... output is written 
     !
     CALL write_output()
     !
     ! ... stdout is reconnected to standard output
     !
     stdout = 6
     !  
     CALL stop_pw( conv_neb )
     !
  ELSE
     !
     CALL init_run()
     !
     istep = 0
     !
     main_loop : DO WHILE ( istep < nstep )
        !
        istep = istep + 1
        !
        ! ... electronic self-consistentcy
        !
        CALL electrons()
        !
        IF ( .NOT. conv_elec ) CALL stop_pw( conv_elec )
        !
        ! ... if requested ions are moved
        !
        CALL ions()
        !
        IF ( conv_ions ) EXIT main_loop
        !
        ! ... the ionic part of the hamiltonian is reinitialized
        !
        CALL hinit1()
        !
     END DO main_loop
     !
     CALL punch()
     !
     CALL stop_pw( conv_ions )
     !
  END IF      
  !
  STOP
  !
9010 FORMAT( /5X,'Current dimensions of program pwscf are:' &
             /5X,'ntypx =',I2,'   npk =',I5,'  lmax =',I2   &
             /5X,'nchix =',I2,'  ndim =',I5,'  nbrx =',I2,' nqfm =',I2 )
  !
END PROGRAM pwscf

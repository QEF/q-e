!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM pwscf
  !! Author: Paolo Giannozzi
  !
  !! Version: v6.1
  !
  !! License: GNU
  !
  !! Summary: Main program calling one or more instances of Plane Wave Self-Consistent Field code
  !
  !! This is the main program for executable "pw.x".  |
  !! * If called as "pw.x -ipi server-address" or "pw.x --ipi server-address",
  !! works in "server" mode, calls [[run_driver]].
  !! * If called as "manypw.x" via a link, works in "manypw" mode, runs many
  !! instances (images) of pw.x (see [[run_manypw]])
  !! * If called as "dist.x" via a link, works in "dry run" mode, computes
  !! distances, angles, neighbors, writes to file "dist.out" and stops. 
  !! Otherwise: see [[run_pwscf]]
  !!
  !! @Note
  !! 10/01/17 Samuel Ponce: Add Ford documentation
  !! @endnote
  !!
  !! @warning
  !! Example of Warning
  !!
  !! @todo
  !! Have automatic parallelisation. 
  !!
  !! @bug
  !! No bug.
  !!
  USE environment,          ONLY : environment_start
  USE mp_global,            ONLY : mp_startup
  USE read_input,           ONLY : read_input_file
  USE command_line_options, ONLY : input_file_, command_line, nimage_
  !
  IMPLICIT NONE
  !
  CHARACTER(len=256) :: srvaddress
  !! Get the address of the server 
  CHARACTER(len=256) :: get_server_address
  !! Get the address of the server 
  INTEGER :: exit_status
  !! Status at exit
  LOGICAL :: use_images
  !! true if running "manypw.x"
  LOGICAL, EXTERNAL :: matches
  !! checks if first string is contained in the second
  !
  CALL mp_startup( start_images=.TRUE., images_only=.TRUE. )
  !
  !
  CALL environment_start( 'PWSCF' )
  !
  ! ... Check if running standalone or in "driver" mode
  !
  srvaddress = get_server_address( command_line ) 
  !
  ! ... Check if running standalone or in "manypw" mode
  !
  use_images = matches( 'manypw.x', command_line )
  !
  ! ... Perform actual calculation
  !
  IF ( TRIM(srvaddress) == ' ' ) THEN
    ! When running standalone:
    IF ( use_images ) THEN
       ! as manypw.x
       CALL run_manypw( )
       CALL run_pwscf( exit_status )
       !
     ELSE
       ! as pw.x
       IF ( nimage_ > 1 ) CALL errore('run_pwscf', &
                          'image parallelization not allowed',1)
       CALL read_input_file( 'PW', input_file_ )
       CALL run_pwscf( exit_status )
       !
    ENDIF
    !
  ELSE
     ! When running as library
     !
     CALL read_input_file('PW+iPi', input_file_ )
     CALL run_driver( srvaddress, exit_status )
     !
  ENDIF
  !
  CALL laxlib_end()
  CALL stop_run( exit_status )
  CALL do_stop( exit_status )
  !
  STOP
  !
END PROGRAM pwscf

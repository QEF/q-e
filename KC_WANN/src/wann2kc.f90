!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------
PROGRAM wann2kc
  !-----------------------------------------------------------------
  !
  !!  This is one of the main drivers of the KC_WANN code 
  !!  It reads the PWSCF and Wannier90 outputs and prepare the 
  !!  subsequent KC calculations. Call kc_Setup.f90 to compute 
  !!  the periodic part of the wannier functions and save these
  !!  on file. 
  !!   
  !!  Code written by Nicola Colonna (EPFL April 2019) 
  !
  USE mp_global,             ONLY : mp_startup,mp_global_end 
  USE environment,           ONLY : environment_start, environment_end
  USE check_stop,            ONLY : check_stop_init
  USE control_kc_wann,       ONLY : calculation
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=9) :: code='KC_WANN'
  !
  ! 1) Initialize MPI, clocks, print initial messages
  CALL mp_startup ( )
  CALL environment_start ( code )
  !
  CALL check_stop_init ( )
  !
  calculation="wann2kc"
  !
  ! 2) Read the input file and the PW output
  CALL kc_readin( ) 
  !
  ! 3) Set up for the KC calculation. 
  CALL kc_setup( )
  !
  CALL clean_pw( .TRUE. )
  CALL close_kc ( ) 
  !
  CALL print_clock_kc ( )
  !
  ! 5) Clean and Close 
  CALL mp_global_end()
  CALL environment_end( code )
  !
END PROGRAM wann2kc

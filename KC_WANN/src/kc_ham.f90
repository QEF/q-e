!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------
PROGRAM kc_ham
  !-----------------------------------------------------------------
  !
  !!  This is one the main drivers of the KC_WANN code to build up the 
  !!  KC hamiltonian in Real Space. It reads the output of 
  !!  a PWscf calculation and the U matrices from W90
  !!   
  !!  Code written by Nicola Colonna and Riccardo de Gennaro (EPFL April 2019) 
  !
  USE mp_global,             ONLY : mp_startup,mp_global_end
  USE environment,           ONLY : environment_start, environment_end
  USE check_stop,            ONLY : check_stop_init
  USE control_kc_wann,       ONLY : calculation, do_bands, write_hr
  USE interpolation,         ONLY : interpolate_ham, dealloc_interpolation
  !
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=9) :: code='KC_WANN'
  !
  ! 1) Initialize MPI, clocks, print initial messages+
  CALL mp_startup ( )
  CALL environment_start ( code )
  !
  CALL check_stop_init ( )
  !
  calculation="ham"
  !
  ! 2) Read the input file and the PW output
  CALL kc_readin( ) 
  !
  ! 3) Set up for the KC calculation. 
  CALL kc_setup_ham( )
  !
  ! 4) Build up the Hamiltonian
  ! 4a) Diagonal term only to 2nd order
  CALL ham_R0_2nd ( )
  ! 4b) Full Hamiltonian to 2nd order 
  CALL koopmans_ham ( )
  !
  ! 5) If do_bands=TRUE interpolate H(k) and prints bands
  IF ( do_bands ) CALL interpolate_ham( )
  !
  ! 6) If write_hr=TRUE write H(R) to file
  IF ( write_hr ) CALL write_hr_to_file( )
  !
  IF (do_bands) CALL dealloc_interpolation( )
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
END PROGRAM kc_ham

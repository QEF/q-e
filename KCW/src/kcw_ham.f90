!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------
SUBROUTINE kcw_ham
  !-----------------------------------------------------------------
  !
  !!  This is one the main subroutines of the KCW code to build up the 
  !!  KC hamiltonian in Real Space. It reads the output of 
  !!  a PWscf calculation and the U matrices from W90
  !!   
  !!  Code written by Nicola Colonna and Riccardo de Gennaro (EPFL April 2019) 
  !
  USE control_kcw,               ONLY : do_bands, write_hr
  USE interpolation,         ONLY : interpolate_ham, dealloc_interpolation
  !
  USE io_rho_xml,            ONLY : write_scf
  USE io_files,              ONLY : prefix, iunwfc
  USE scf,                   ONLY : rho
  USE lsda_mod,              ONLY : nspin
  USE units_lr,              ONLY : iuwfc
  !
  !
  IMPLICIT NONE
  !
  ! 1) Set up for the KC calculation. 
  CALL kcw_setup_ham( )
  !
  ! 2) Build up the Hamiltonian
  ! 2a) Diagonal term only to 2nd order
  !     OBSOLETE: inside koopmans_ham this is triggered by "on_site_only": FIXME
  !CALL ham_R0_2nd ( )
  ! 2b) Full Hamiltonian to 2nd order 
  CALL koopmans_ham ( )
  !
  ! 3) If do_bands=TRUE interpolate H(k) and prints bands
  IF ( do_bands ) CALL interpolate_ham( )
  !
  ! 4) If write_hr=TRUE write H(R) to file
  IF ( write_hr ) CALL write_hr_to_file( )
  !
  IF (do_bands) CALL dealloc_interpolation( )
  ! 
  ! WRITE data file
  iunwfc = iuwfc
  prefix = TRIM(prefix)//"_kcw"
  CALL write_scf(rho, nspin)
  !CALL punch('config-only')
  CALL punch ('all')
  !
  CALL clean_pw( .TRUE. )
  CALL close_kcw ( ) 
  !
END SUBROUTINE kcw_ham

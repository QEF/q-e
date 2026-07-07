!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------
SUBROUTINE wann2kcw
  !-----------------------------------------------------------------
  !
  !!  This is one of the main drivers of the KCW code 
  !!  It reads the PWSCF and Wannier90 outputs and prepare the 
  !!  subsequent KC calculations. Call kcw_setup.f90 to compute 
  !!  the periodic part of the wannier functions and save these
  !!  on file. 
  !!   
  !!  Code written by Nicola Colonna (EPFL April 2019) 
  !
  USE io_files,           ONLY : iunwfc, nwordwfc
  USE wvfct,              ONLY : nbnd, npwx
  USE control_kcw,        ONLY : spin_component, iuwfc_wann, num_wann
  USE noncollin_module,   ONLY : npol
  USE pw_restart_new,     ONLY : write_collected_wfc
  !
  IMPLICIT NONE
  !
  INTEGER  :: iunwfc_, nwordwfc_, nbnd_
  !
  !
  ! 1) Set up for the KC calculation. 
  CALL kcw_setup( )
  ! 
  ! 2) Save MLWF on file 
  ! To use write_collected_wfc, we need to redefine some global variables
  ! Here we store ...
  nbnd_     = nbnd
  nwordwfc_ = nwordwfc
  iunwfc_   = iunwfc
  ! ... overwrite ...
  nbnd      = num_wann
  nwordwfc  = num_wann * npwx * npol
  iunwfc    = iuwfc_wann
  CALL write_collected_wfc (spin_component)
  ! ... restore original values
  nbnd      = nbnd_
  nwordwfc  = nwordwfc_
  iunwfc    = iunwfc_
  !
  CALL clean_pw( .TRUE. )
  CALL close_kcw ( ) 
  !
END SUBROUTINE wann2kcw

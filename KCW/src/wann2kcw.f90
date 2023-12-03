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
  USE io_kcw,             ONLY : write_mlwf
  !
  IMPLICIT NONE
  !
  !
  ! 1) Set up for the KC calculation. 
  CALL kcw_setup( )
  ! 
  ! 2) Save MLWF on file 
  CALL write_mlwf( ) 
  !
  CALL clean_pw( .TRUE. )
  CALL close_kcw ( ) 
  !
END SUBROUTINE wann2kcw

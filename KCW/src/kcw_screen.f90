!
! Copyright (C) 2003-2024 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------
SUBROUTINE kcw_screen
  !-----------------------------------------------------------------
  !
  !!  This is the one main driver of the kcw.x code. It computes the 
  !!  screaning parameters alpha for KOOPMANS. It reads the ouput of 
  !!  a PWscf calculation and the U matrices from W90 and computes
  !!  the orbital dependent screening coefficients as described in 
  !!  N. Colonna et al. JCTC 14, 2549 (2018) 
  !!  https://pubs.acs.org/doi/10.1021/acs.jctc.7b01116
  !!  and
  !! 
  !!  Code written by Nicola Colonna (EPFL April 2019) 
  !!
  !!  Non-collinear code written in 2022-24 by 
  !!  Antimo Marrazzo (SISSA, UniTS) and Nicola Colonna (PSI)
  !
  USE klist,                 ONLY : nkstot
  USE lsda_mod,              ONLY : nspin
  USE control_kcw,           ONLY : nkstot_eff
  !
  IMPLICIT NONE
  !
  ! 3) Set up for the KC calculation. 
  CALL kcw_setup_screen( )
  !
  ! 4) Compute the screening coefficient via Linear Respoonse
  CALL screen_coeff ( )
  ! 
  CALL clean_pw( .TRUE. )
  CALL close_kcw ( ) 
  !
  IF (nkstot_eff .gt. 1) CALL print_clock_pw ( )
  !
END SUBROUTINE kcw_screen

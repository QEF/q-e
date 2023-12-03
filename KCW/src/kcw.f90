!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------
PROGRAM kcw
  !-----------------------------------------------------------------
  !
  !! This is the main driver of the kcw.x code, an implementation of koopmans
  !! functionals based on DFPT and Wannier functions. The code reads the output
  !! of PWSCF and Wannier90 and performe a Koopmans calculation in a perturbative
  !! way. It performe several task depending on the vaule of the variable "calculation". 
  !! 1) calculation=wann2kcw: interface between PWSCF and W90, and KCW. 
  !! 2) calculation=screen: calculates the screening coefficients as described in 
  !!    N. Colonna et al. J. Chem. Theory Comput. 14, 2549 (2018) 
  !!    N. Colonna et al. J. Chem. Theory Comput. 18, 5435 (2022)
  !! 3) calculation=ham: compute, interpolate and diagonalize the KC hamiltonian 
  !!
  !!  Code written by Nicola Colonna (EPFL April 2019) 
  !
  USE mp_global,         ONLY : mp_startup,mp_global_end 
  USE environment,       ONLY : environment_start, environment_end
  USE control_kcw,       ONLY : calculation
  USE mp_global,         ONLY : mp_startup
  USE check_stop,        ONLY : check_stop_init
  USE coulomb,           ONLY : setup_coulomb
  USE control_flags,     ONLY : use_gpu
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=9) :: code='KCW'
  LOGICAL,EXTERNAL :: check_gpu_support 
  !
  use_gpu = check_gpu_support()
  IF(use_gpu) Call errore('KCW', 'KCW with GPU NYI', 1)
  !
  ! 1) Initialize MPI, clocks, print initial messages
  CALL mp_startup ( )
  !
  CALL header ()
  !
  CALL environment_start ( code )
  !
  CALL check_stop_init ()
  !
  ! 2) Read the input file and the PW outputs
  CALL kcw_readin( ) 
  !
  IF (calculation == 'cc') call setup_coulomb()
  IF (calculation == 'wann2kcw') CALL wann2kcw ( )
  IF (calculation == 'screen')   CALL kcw_screen ( )
  IF (calculation == 'ham' )     CALL kcw_ham ()
  !
  CALL print_clock_kcw ( )
  !
  ! 5) Clean and Close 
  CALL mp_global_end()
  CALL environment_end( code )
  !
END PROGRAM kcw


!-------------------------------------------------------------------
SUBROUTINE header 
  !-----------------------------------------------------------------
  !
  USE io_global,         ONLY : stdout, ionode
  IMPLICIT NONE
 
  IF (ionode) THEN 
  WRITE( stdout, '(/5x,"=--------------------------------------------------------------------------------=")')
  WRITE( stdout,*) "                     :::    :::           ::::::::         :::       ::: "
  WRITE( stdout,*) "                    :+:   :+:           :+:    :+:        :+:       :+:  "
  WRITE( stdout,*) "                   +:+  +:+            +:+               +:+       +:+   "
  WRITE( stdout,*) "                  +#++:++             +#+               +#+  +:+  +#+    "
  WRITE( stdout,*) "                 +#+  +#+            +#+               +#+ +#+#+ +#+     "
  WRITE( stdout,*) "                #+#   #+#           #+#    #+#         #+#+# #+#+#       " 
  WRITE( stdout,*) "               ###    ###           ########           ###   ###         "
  WRITE( stdout, '(/5x,"  Koopmans functional implementation based on DFPT; please cite this program as")')
  WRITE( stdout, '(/5x,"   N.Colonna, R. De Gannaro, E. Linscott, and N. Marzari, JCTC 18, 5435 (2022) ")')
  WRITE( stdout, '( 5x,"=--------------------------------------------------------------------------------=")')  
  ENDIF
  !
END SUBROUTINE header

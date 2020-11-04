!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE engine_to_path_nat()
  !-----------------------------------------------------------------------------
  !
  !
  USE kinds,         ONLY : DP
  !
  USE input_parameters, ONLY : nat 
  USE path_input_parameters_module, ONLY : nat_ => nat 
  !
  !
  IMPLICIT NONE
  !
  nat_ = nat
  !
  !
  RETURN
  !
END SUBROUTINE engine_to_path_nat
!

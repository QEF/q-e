!
! Copyright (C) 2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE engine_to_path_tot_charge(idx)
  !-----------------------------------------------------------------------------
  !
  !
  USE kinds,         ONLY : DP
  !
  USE klist, ONLY : tot_charge
  USE path_input_parameters_module, ONLY : tot_charge_ => tot_charge
  !
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: idx
  !
  tot_charge_(idx) = tot_charge
  !
  !
  RETURN
  !
END SUBROUTINE engine_to_path_tot_charge
!

!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE engine_to_path_alat()
  !-----------------------------------------------------------------------------
  !
  !
  USE kinds,         ONLY : DP
  !
  USE cell_base, ONLY : alat 
  USE path_input_parameters_module, ONLY : alat_ => alat
  !
  !
  IMPLICIT NONE
  !
  alat_ = alat
  !
  !
  RETURN
  !
END SUBROUTINE engine_to_path_alat
!

!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE engine_to_path_fix_atom_pos()
  !-----------------------------------------------------------------------------
  !
  !
  USE kinds,         ONLY : DP
  !
  USE ions_base, ONLY : if_pos
  USE path_variables, ONLY : fix_atom_pos
  USE path_input_parameters_module, ONLY : nat 
  !
  ! ... "path" specific
  !
  !
  IMPLICIT NONE
  !
  ! set_my_if_pos
  !
  allocate(fix_atom_pos(3,nat))
  fix_atom_pos(:,:) = 1
  fix_atom_pos(:,:) = if_pos(:,:)
  !
  RETURN
  !
END SUBROUTINE engine_to_path_fix_atom_pos
!

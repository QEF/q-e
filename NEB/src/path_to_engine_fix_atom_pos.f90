!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE path_to_engine_fix_atom_pos()
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
  IMPLICIT NONE
  !
  ! set_my_if_pos
  !
  if(.not.allocated(if_pos)) allocate(if_pos(3,nat))
  if_pos(:,:) = 1
  if_pos(:,:) = fix_atom_pos(:,:)
  !
  !
  RETURN
  !
END SUBROUTINE path_to_engine_fix_atom_pos
!

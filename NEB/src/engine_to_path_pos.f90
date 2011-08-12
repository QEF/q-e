!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE engine_to_path_pos(index)
  !-----------------------------------------------------------------------------
  !
  !
  USE kinds,         ONLY : DP
  !
  USE path_input_parameters_module, ONLY : input_images
  !
  USE path_input_parameters_module, ONLY : nat, alat 
  !
  !
  USE path_input_parameters_module, ONLY : pos
  USE ions_base, ONLY : tau
  !
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: index
  !
  ! set_my_tau
  ! il tau e gia in units di pw
  !
  if(.not.allocated(pos)) allocate(pos(3*nat,input_images))
  pos(:,index) = 0.0_dp
  !
  ! ... note that this positions array is in Bohr
  !
  pos(1:3*nat,index) = reshape( tau, (/ 3 * nat /) ) * alat
  !
  RETURN
  !
END SUBROUTINE engine_to_path_pos
!

!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE engine_to_path_pos(idx)
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
  USE path_input_parameters_module, ONLY : pos, typ
  USE ions_base, ONLY : tau, ityp
  !
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: idx
  !
  ! set_my_tau
  ! il tau e gia in units di pw
  !
  ! is this really necessary? (GS)
  if(.not.allocated(pos)) allocate(pos(3*nat,input_images))
  pos(:,idx) = 0.0_dp
  !
  ! ... note that this positions array is in Bohr
  !
  pos(1:3*nat,idx) = reshape( tau, (/ 3 * nat /) ) * alat
  !
  ! consistency check on atomic type, just to be sure... (GS)
  !
  IF ( idx == 1 ) THEN
     ! is this really necessary? (GS)
     if(.not.allocated(typ)) allocate(typ(nat))
     typ = ityp(:)
  ELSE
     IF ( ANY( typ .NE. ityp ) ) CALL errore("engine_to_path_pos", &
        "inconsistency of atomic species", idx )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE engine_to_path_pos
!

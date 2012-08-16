!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine set_irr_nosym_new (u, npert, nirr)
  !---------------------------------------------------------------------
  !
  !     This routine substitutes set_irr when there are no symmetries.
  !     The irreducible representations are all one dimensional and
  !     we set them to the displacement of a single atom in one direction
  !
  USE kinds, only : DP
  USE ions_base, ONLY : nat
  USE modes, ONLY : num_rap_mode, name_rap_mode
  USE control_ph, ONLY : search_sym
  IMPLICIT NONE
  !
  INTEGER, INTENT(OUT) ::  npert (3 * nat), nirr
  ! output: the dimension of each representation
  ! output: the number of representation
  COMPLEX(DP), INTENT(OUT) :: u( 3 * nat, 3 * nat )
  !
  integer :: imode, irr
  ! counter on modes
  ! counter on representations
  !
  !
  nirr = 3 * nat
  npert = 1

  u = (0.d0, 0.d0)
  do imode = 1, 3 * nat
     u (imode, imode) = (1.d0, 0.d0)
  enddo
  IF (search_sym) THEN
     DO imode = 1, 3 * nat
        num_rap_mode(imode)=1
        name_rap_mode(imode)='A'
     END DO
  ENDIF

  return
end subroutine set_irr_nosym_new

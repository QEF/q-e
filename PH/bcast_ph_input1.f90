!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine bcast_ph_input1
  !-----------------------------------------------------------------------
  !
#ifdef __PARA
#include "machine.h"

  use pwcom
  use phcom
  use mp, only: mp_bcast
  implicit none
  integer :: root = 0

  !
  ! integers
  !
  call mp_bcast (nat_todo, root)
  if (nat_todo.gt.0) then
     call mp_bcast (atomo, root)
  endif
  call mp_bcast (nrapp, root)
  if (nrapp.gt.0) then
     call mp_bcast (list, root)
  endif
#endif
  return
end subroutine bcast_ph_input1

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
#include "f_defs.h"

  use pwcom
  use phcom
  use mp, only: mp_bcast
  USE io_global, ONLY : ionode_id
  implicit none

  !
  ! integers
  !
  call mp_bcast (nat_todo, ionode_id)
  if (nat_todo.gt.0) then
     call mp_bcast (atomo, ionode_id)
  endif
  call mp_bcast (nrapp, ionode_id)
  if (nrapp.gt.0) then
     call mp_bcast (list, ionode_id)
  endif
#endif
  return
end subroutine bcast_ph_input1

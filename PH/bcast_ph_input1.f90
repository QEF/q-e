!
! Copyright (C) 2001 PWSCF group
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
#ifdef PARA
#include "machine.h"

use pwcom
use parameters, only : DP
use phcom
implicit none

include 'mpif.h'

integer :: root, errcode
root = 0
call MPI_barrier (MPI_COMM_WORLD, errcode)
call error ('bcast_ph_input1', 'at barrier ', errcode)
!
! integers
!
call MPI_bcast (nat_todo, 1, MPI_INTEGER, root, MPI_COMM_WORLD, &
 errcode)
call error ('bcast_ph_input1', 'bcasting nat_todo ', errcode)

call MPI_barrier (MPI_COMM_WORLD, errcode)
if (nat_todo.gt.0) then
   call MPI_bcast (atomo, nat_todo, MPI_INTEGER, root, &
    MPI_COMM_WORLD, errcode)
   call error ('bcast_ph_input1', 'bcasting atomo ', errcode)

endif
call MPI_bcast (nrapp, 1, MPI_INTEGER, root, MPI_COMM_WORLD, &
 errcode)
call error ('bcast_ph_input1', 'bcasting nrapp ', errcode)

call MPI_barrier (MPI_COMM_WORLD, errcode)
if (nrapp.gt.0) then
   call MPI_bcast (list, nrapp, MPI_INTEGER, root, MPI_COMM_WORLD, &
    errcode)
   call error ('bcast_ph_input1', 'bcasting list ', errcode)
endif
#endif
return
end subroutine bcast_ph_input1

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine d3_recover (ilab, isw)
  !-----------------------------------------------------------------------
  !
  !  isw = +1 Writes d3dyn in a file for possible recover
  !  isw = -1 Starts a recover run
#include "machine.h"
  use pwcom
  use phcom
  use d3com
#ifdef PARA
  use para
#endif
  implicit none

#ifdef PARA
  include 'mpif.h'
#endif
  integer :: ilab, isw, root, iaux, errcode
  logical :: exst

  iunrec = 98
  if (isw.eq.1) then
#ifdef PARA
     call MPI_barrier (MPI_COMM_WORLD, errcode)
     call error ('d3_recover', 'at barrier', errcode)
     if (me.ne.1.or.mypool.ne.1) return
#endif

     call seqopn (iunrec, 'recv_d3', 'unformatted', exst)
     if (ilab.le.4) then
        write (iunrec) ilab
     else
        write (iunrec) ilab, d3dyn

     endif

     close (unit = iunrec, status = 'keep')
  elseif (isw.eq. - 1) then
#ifdef PARA
     call MPI_barrier (MPI_COMM_WORLD, errcode)
     call error ('d3_recover', 'at barrier', errcode)
     if (me.ne.1.or.mypool.ne.1) goto 100
#endif

     call seqopn (iunrec, 'recv_d3', 'unformatted', exst)
     read (iunrec) ilab
     if (ilab.ge.5) then
        rewind (iunrec)
        read (iunrec) ilab, d3dyn

     endif

     close (unit = iunrec, status = 'keep')
#ifdef PARA

100  continue
     root = 0
     call MPI_barrier (MPI_COMM_WORLD, errcode)

     call error ('d3_recover', 'at barrier2', errcode)
     iaux = 2 * 27 * nat * nat * nat
     call MPI_bcast (d3dyn, iaux, MPI_REAL8, root, MPI_COMM_WORLD, &
          errcode)
     call error ('d3_recover', 'at bcast1', errcode)
     call MPI_bcast (ilab, 1, MPI_INTEGER, root, MPI_COMM_WORLD, &
          errcode)

     call error ('d3_recover', 'at bcast2', errcode)
#endif

  endif
  return
end subroutine d3_recover

!
! Copyright (C) 2001-2003 PWSCF group
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
#ifdef __PARA
  use para
  use mp, only: mp_bcast
#endif
  implicit none
  integer :: ilab, isw
  integer :: root = 0
  logical :: exst

  iunrec = 98
  if (isw.eq.1) then
#ifdef __PARA
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
#ifdef __PARA
     if (me.ne.1.or.mypool.ne.1) goto 100
#endif

     call seqopn (iunrec, 'recv_d3', 'unformatted', exst)
     read (iunrec) ilab
     if (ilab.ge.5) then
        rewind (iunrec)
        read (iunrec) ilab, d3dyn

     endif

     close (unit = iunrec, status = 'keep')
#ifdef __PARA

100  continue

     call mp_bcast (d3dyn, root)
     call mp_bcast (ilab, root)
#endif

  endif
  return
end subroutine d3_recover

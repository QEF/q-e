!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine stop_ph (flag)
  !-----------------------------------------------------------------------
  !
  ! Close all files and synchronize processes before stopping.
  ! Called at the end of the run with flag=.true. (removes 'recover')
  ! or during execution with flag=.false. (does not remove 'recover')
  !

  use pwcom
  use parameters, only : DP
  use phcom
#ifdef __PARA
  use para
#endif
  implicit none
#ifdef __PARA
  include 'mpif.h'
  integer :: info
#endif
  logical :: flag, exst

  close (unit = iuwfc, status = 'keep')
  close (unit = iudwf, status = 'keep')
  close (unit = iubar, status = 'keep')
#ifdef __PARA
  if (me.ne.1) goto 100
#endif
  if (fildrho.ne.' ') close (unit = iudrho, status = 'keep')
#ifdef __PARA
100 continue
#endif

  if (flag) then
     call seqopn (iunrec, 'recover','unformatted',exst)
     close (unit=iunrec,status='delete')
  end if
  close (unit = iunigk, status = 'delete')
  call print_clock_ph
  call show_memory ()
#ifdef __PARA
  call mpi_barrier (MPI_COMM_WORLD, info)

  call mpi_finalize (info)
#endif
#ifdef __T3E
  !
  ! set streambuffers off
  !

  call set_d_stream (0)
#endif
  if (flag) then
     stop
  else
     stop 1
  endif

end subroutine stop_ph

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
subroutine stop_pw (flag)
  !--------------------------------------------------------------------
  !
  ! Close all files and synchronize processes before stopping.
  ! Called at the end of the run with flag=.true. (removes 'restart')
  ! or during execution with flag=.false. (does not remove 'restart')
  !
  use pwcom
  use io, only : prefix
  use mp, only : mp_end
  logical :: flag
#ifdef __PARA
  include 'mpif.h'
  integer :: info
#endif
  logical exst
  !
  !  iunwfc contains wavefunctions and is kept open during
  !  the execution - close and save the file
  !
  close (unit = iunwfc, status = 'keep')
  if (flag) then
     !
     !  all other files must be reopened and removed
     !
     call seqopn (iunres, 'restart','unformatted',exst)
     close (unit=iunres,status='delete')
     call seqopn (4, trim(prefix)//'.bfgs','unformatted',exst)
     close (unit=4,status='delete')
     call seqopn (4, trim(prefix)//'.md','formatted',exst)
     close (unit=4,status='delete')
  endif
  !
  !  iunigk is kept open during the execution - close and remove
  !
  close (unit = iunigk, status = 'delete')
  call print_clock_pw

  call show_memory ()
#ifdef __PARA
  call mpi_barrier (MPI_COMM_WORLD, info)

  ! call mpi_finalize (info)
#endif
 
  call mp_end()

#ifdef __T3E
  !
  ! set streambuffers off
  !

  call set_d_stream (0)
#endif

  call clean_pw

  if (flag) then
     stop
  else
     stop 1
  endif
end subroutine stop_pw
!
!--------------------------------------------------------------------

subroutine closefile
  !--------------------------------------------------------------------
  !
  ! Close all files and synchronize processes before stopping
  ! Called by "sigcatch" when it receives a signal
  !
  write (6, '(5x,"Signal Received, stopping ... ")')

  call stop_pw (.false.)
  return
end subroutine closefile

!--------------------------------------------------------------------

subroutine cpflush
  !--------------------------------------------------------------------
  !
  ! TEMP: compatibility with Car-Parrinello code
  !
  print *, "what am i doing in cpflush ?"
  call stop_pw (.false.)
  return
end subroutine cpflush





!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine stop_d3 (flag)
!-----------------------------------------------------------------------
!
!    This routine closes all files before stopping
!    flag is no longer used
!
  use pwcom
  use phcom
  use d3com
#ifdef __PARA
  use para
#endif
  implicit none
#ifdef __PARA
  include 'mpif.h'
  integer :: info
#endif
  logical :: flag

  close (unit = iuwfc, status = 'keep')
  close (unit = iubar, status = 'keep')
  close (unit = iudwf, status = 'keep')
#ifdef __PARA
  if (me.ne.1) goto 100
#endif
  close (unit = iudrho, status = 'keep')
  if (.not.lgamma) close (unit = iud0rho, status = 'keep')
#ifdef __PARA
100 continue
#endif
  close (unit = iunigk, status = 'delete')
  if (.not.lgamma) then
     close (unit = iud0qwf, status = 'keep')
     close (unit = iudqwf, status = 'keep')
  endif
  close (unit = iupdqvp, status = 'keep')
  if (.not.lgamma) close (unit = iupd0vp, status = 'keep')
  if (degauss.ne.0.d0) then
     close (unit = iudpdvp_1, status = 'keep')
     if (.not.lgamma) then
        close (unit = iudpdvp_2, status = 'keep')
        close (unit = iudpdvp_3, status = 'keep')
     endif
  endif
  call print_clock_d3
  call show_memory ()

#ifdef __PARA
  call mpi_barrier (MPI_COMM_WORLD, info)
  call mpi_finalize (info)
#endif

  stop
  return
end subroutine stop_d3

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine bcast_d3_input
  !-----------------------------------------------------------------------
  !
  !     In this routine the first processor send the input to all
  !     the other processors
  !
  !
#ifdef PARA
#include "machine.h"
  use pwcom
  use phcom
  use d3com
  use para
  use io

  implicit none
  include 'mpif.h'
  integer :: root, errcode
  character (len=512) :: buffer
  integer :: buffer_t3d (64)

  equivalence (buffer, buffer_t3d)
  root = 0
  call MPI_barrier (MPI_COMM_WORLD, errcode)
  call error ('bcast_ph_input', 'at barrier', errcode)
  ! logicals
  call MPI_bcast (lgamma, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, &
       errcode)
  call error ('bcast_ph_input', 'at bcast1', errcode)
  call MPI_bcast (wraux, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, &
       errcode)
  call error ('bcast_d3_input', 'at bcast51', errcode)
  call MPI_bcast (recv, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, &
       errcode)
  call error ('bcast_d3_input', 'at bcast52', errcode)
  call MPI_bcast (testflag, 50, MPI_LOGICAL, root, MPI_COMM_WORLD, &
       errcode)
  call error ('bcast_d3_input', 'at bcast53', errcode)
  !
  ! integers
  !
  call MPI_bcast (iverbosity, 1, MPI_INTEGER, root, MPI_COMM_WORLD, &
       errcode)
  call error ('bcast_ph_input', 'at bcast9', errcode)
  call MPI_bcast (testint, 50, MPI_INTEGER, root, MPI_COMM_WORLD, &
       errcode)
  call error ('bcast_d3_input', 'at bcast91', errcode)
  call MPI_bcast (q0mode_todo, 300, MPI_INTEGER, root, &
       MPI_COMM_WORLD, errcode)
  call error ('bcast_d3_input', 'at bcast92', errcode)
  call MPI_bcast (istop, 1, MPI_INTEGER, root, MPI_COMM_WORLD, &
       errcode)
  call error ('bcast_d3_input', 'at bcast93', errcode)
  !
  ! real*8
  !
  call MPI_bcast (amass, ntypx, MPI_REAL8, root, MPI_COMM_WORLD, &
       errcode)
  call error ('bcast_ph_input', 'at bcast12', errcode)
  call MPI_bcast (xq, 3, MPI_REAL8, root, MPI_COMM_WORLD, errcode)
  call error ('bcast_ph_input', 'at bcast14', errcode)
  call MPI_bcast (ethr_ph, 1, MPI_REAL8, root, MPI_COMM_WORLD, &
       errcode)
  call error ('bcast_d3_input', 'at bcast16', errcode)
  call MPI_bcast (testreal, 50, MPI_REAL8, root, MPI_COMM_WORLD, &
       errcode)
  call error ('bcast_d3_input', 'at bcast17', errcode)
  !
  ! characters
  !
  if (me.eq.1) write (buffer, '(a75,a14,3a50,a35)') title_ph, &
       filpun, fildyn, fildrho, fild0rho, tmp_dir
  call MPI_bcast (buffer_t3d, 512, MPI_CHARACTER, root, &
       MPI_COMM_WORLD, errcode)
  call error ('bcast_ph_input', 'at bcast character', errcode)
  if (me.ne.1) read (buffer, '(a75,a14,3a50,a35)') title_ph, filpun, &
       fildyn, fildrho, fild0rho, tmp_dir
#endif
  return
end subroutine bcast_d3_input

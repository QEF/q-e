!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine bcast_ph_input
  !-----------------------------------------------------------------------
  !
  !     In this routine the first processor send the input to all
  !     the other processors
  !
  !
#ifdef __PARA
#include "machine.h"

  use pwcom
  use parameters, only : DP
  use phcom
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
  call errore ('bcast_ph_input', 'at barrier', errcode)
  ! logicals
  call MPI_bcast (lgamma, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, &
       errcode)
  call errore ('bcast_ph_input', 'at bcast1', errcode)
  call MPI_bcast (epsil, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, &
       errcode)
  call errore ('bcast_ph_input', 'at bcast2', errcode)
  call MPI_bcast (trans, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, &
       errcode)
  call errore ('bcast_ph_input', 'at bcast3', errcode)
  call MPI_bcast (zue, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, &
       errcode)
  call errore ('bcast_ph_input', 'at bcast4', errcode)
  call MPI_bcast (reduce_io, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, &
       errcode)
  call errore ('bcast_ph_input', 'at bcast5', errcode)
  call MPI_bcast (elph, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, &
       errcode)
  call errore ('bcast_ph_input', 'at bcast5.1', errcode)
  !
  ! integers
  !
  call MPI_bcast (niter_ph, 1, MPI_INTEGER, root, MPI_COMM_WORLD, &
       errcode)
  call errore ('bcast_ph_input', 'at bcast6', errcode)
  call MPI_bcast (nmix_ph, 1, MPI_INTEGER, root, MPI_COMM_WORLD, &
       errcode)
  call errore ('bcast_ph_input', 'at bcast7', errcode)
  call MPI_bcast (maxirr, 1, MPI_INTEGER, root, MPI_COMM_WORLD, &
       errcode)
  call errore ('bcast_ph_input', 'at bcast8', errcode)
  call MPI_bcast (iverbosity, 1, MPI_INTEGER, root, MPI_COMM_WORLD, &
       errcode)
  call errore ('bcast_ph_input', 'at bcast9', errcode)
  !
  ! real*8
  !
  call MPI_bcast (tr2_ph, 1, MPI_REAL8, root, MPI_COMM_WORLD, &
       errcode)
  call errore ('bcast_ph_input', 'at bcast11', errcode)
  call MPI_bcast (amass, ntypx, MPI_REAL8, root, MPI_COMM_WORLD, &
       errcode)
  call errore ('bcast_ph_input', 'at bcast12', errcode)
  call MPI_bcast (alpha_mix, maxter, MPI_REAL8, root, &
       MPI_COMM_WORLD, errcode)
  call errore ('bcast_ph_input', 'at bcast13', errcode)
  call MPI_bcast (xq, 3, MPI_REAL8, root, MPI_COMM_WORLD, errcode)
  call errore ('bcast_ph_input', 'at bcast14', errcode)
  call MPI_bcast (time_max, 1, MPI_REAL8, root, MPI_COMM_WORLD, &
       errcode)
  call errore ('bcast_ph_input', 'at bcast15', errcode)
  !
  ! characters
  !
  if (me.eq.1) write (buffer, '(a75,5a14,2a55)') title_ph, filpun, &
       filelph, fildyn, fildvscf, fildrho, tmp_dir, prefix
  call MPI_bcast (buffer_t3d, 512, MPI_CHARACTER, root, &
       MPI_COMM_WORLD, errcode)
  call errore ('bcast_ph_input', 'at bcast character', errcode)
  if (me.ne.1) read (buffer, '(a75,5a14,2a55)') title_ph, filpun, &
       filelph, fildyn, fildvscf, fildrho, tmp_dir, prefix
#endif
  return
end subroutine bcast_ph_input

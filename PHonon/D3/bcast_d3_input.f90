!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine bcast_d3_input
  !-----------------------------------------------------------------------
  !
  !     In this routine the first processor sends input data to all
  !     the other processors
  !
  !
#ifdef __MPI
  use pwcom
  use qpoint, ONLY: xq
  use phcom
  use d3com
  use mp, only: mp_bcast
  use mp_world, only: world_comm
  use io_files, only: prefix, tmp_dir
  use ions_base, only: amass
  use control_flags, only: iverbosity
  use run_info, only: title

  implicit none
  integer :: root = 0
  !
  ! logicals
  !
  call mp_bcast (lgamma, root, world_comm)
  call mp_bcast (wraux, root, world_comm)
  call mp_bcast (recv, root, world_comm)
  call mp_bcast (testflag,root, world_comm)
  !
  ! integers
  !
  call mp_bcast (iverbosity, root, world_comm)
  call mp_bcast (testint, root, world_comm)
  call mp_bcast (q0mode_todo, root, world_comm)
  call mp_bcast (istop, root, world_comm)
  !
  ! real*8
  !
  call mp_bcast (amass, root, world_comm)
  call mp_bcast (xq, root, world_comm)
  call mp_bcast (ethr_ph, root, world_comm)
  call mp_bcast (testreal, root, world_comm)
  !
  ! characters
  !
  call mp_bcast (title, root, world_comm)
  call mp_bcast (fildyn, root, world_comm)
  call mp_bcast (fildrho, root, world_comm)
  call mp_bcast (fild0rho, root, world_comm)
  call mp_bcast (tmp_dir, root, world_comm)
  call mp_bcast (prefix, root, world_comm)
#endif
  return
end subroutine bcast_d3_input

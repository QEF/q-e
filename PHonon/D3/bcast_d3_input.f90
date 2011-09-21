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
#ifdef __PARA
  use pwcom
  use phcom
  use d3com
  use mp, only: mp_bcast
  use io_files, only: prefix, tmp_dir
  use ions_base, only: amass
  use control_flags, only: iverbosity
  use run_info, only: title

  implicit none
  integer :: root = 0
  !
  ! logicals
  !
  call mp_bcast (lgamma, root)
  call mp_bcast (wraux, root)
  call mp_bcast (recv, root)
  call mp_bcast (testflag,root)
  !
  ! integers
  !
  call mp_bcast (iverbosity, root)
  call mp_bcast (testint, root)
  call mp_bcast (q0mode_todo, root)
  call mp_bcast (istop, root)
  !
  ! real*8
  !
  call mp_bcast (amass, root)
  call mp_bcast (xq, root)
  call mp_bcast (ethr_ph, root)
  call mp_bcast (testreal, root)
  !
  ! characters
  !
  call mp_bcast (title, root)
  call mp_bcast (fildyn, root)
  call mp_bcast (fildrho, root)
  call mp_bcast (fild0rho, root)
  call mp_bcast (tmp_dir, root)
  call mp_bcast (prefix, root)
#endif
  return
end subroutine bcast_d3_input

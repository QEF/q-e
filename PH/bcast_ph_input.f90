!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine bcast_ph_input ( )
  !-----------------------------------------------------------------------
  !
  !     In this routine the first processor sends the input to all
  !     the other processors
  !
  !
#ifdef __PARA
#include "f_defs.h"

  use pwcom
  use phcom
  use mp, only: mp_bcast
  use io_files
  USE control_flags, only: iverbosity, reduce_io
  USE ramanm, ONLY: lraman, elop, dek, eth_rps, eth_ns
  USE input_parameters, ONLY: max_seconds
  USE ions_base,     ONLY : amass

  implicit none
  integer :: root = 0
  !
  ! logicals
  !
  call mp_bcast (lgamma, root)
  call mp_bcast (epsil, root)
  call mp_bcast (trans, root)
  call mp_bcast (zue, root)
  call mp_bcast (reduce_io, root)
  call mp_bcast (elph, root)
  call mp_bcast (lnscf, root)
  call mp_bcast (ldisp, root)
  call mp_bcast (lraman, root)
  call mp_bcast (elop, root)
  call mp_bcast (recover, root)
  !
  ! integers
  !
  call mp_bcast (niter_ph, root)
  call mp_bcast (nmix_ph, root)
  call mp_bcast (maxirr, root)
  call mp_bcast (iverbosity, root)
  !
  ! real*8
  !
  call mp_bcast (tr2_ph, root)
  call mp_bcast (eth_rps, root)
  call mp_bcast (eth_ns, root)
  call mp_bcast (amass, root)
  call mp_bcast (alpha_mix, root)
  call mp_bcast (xq, root)
  call mp_bcast (max_seconds, root)
  call mp_bcast (dek, root)
  !
  ! characters
  !
  call mp_bcast (title, root)
  call mp_bcast (filelph, root)
  call mp_bcast (fildyn, root)
  call mp_bcast (fildvscf, root)
  call mp_bcast (fildrho, root)
  call mp_bcast (tmp_dir, root)
  call mp_bcast (prefix, root)
#endif
  return
end subroutine bcast_ph_input

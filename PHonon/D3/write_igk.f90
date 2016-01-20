!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
subroutine write_igk
  !
  use pwcom
  use phcom
  USE io_files, ONLY : iunigk

  use qpoint,     ONLY : npwq, nksq, igkq
  use control_lr, ONLY : lgamma

  implicit none

  if (nksq.ne.1) return
  rewind (unit = iunigk)
  write (iunigk) npw, igk

  if (.not.lgamma) write (iunigk) npwq, igkq
  return
end subroutine write_igk

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module io
  use parameters, only: ntypx
  character(len=80) :: tmp_dir  ! directory for temporary files
  character(len=80) :: prefix   ! prepended to file names
  character(len=3)  :: nd_nmbr  ! node number (used only in parallel case)
  character(len=80) :: pseudo_dir
  character(len=80) :: pseudop(ntypx)
end module io

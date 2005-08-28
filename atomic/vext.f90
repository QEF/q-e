!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
function vext(r)
  use kinds, only : DP
  implicit none
  real(DP) ::    vext,r 
  vext=0.0_dp
  return
end function vext

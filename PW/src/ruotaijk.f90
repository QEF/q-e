!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine ruotaijk (s, ftau, i, j, k, nr1, nr2, nr3, ri, rj, rk)
  !----------------------------------------------------------------------
  !
  !    This routine computes the rotated of the point i,j,k throught
  !    the symmetry (s,f). Then it computes the equivalent point
  !    on the original mesh
  !
  !
  USE kinds
  implicit none
  !
  !    first the dummy variables
  !
  integer :: s (3, 3), ftau (3), i, j, k, nr1, nr2, nr3, ri, rj, rk
  ! input: the rotation matrix
  ! input: the fractionary translation
  !   !   input: the point to rotate

  ! /
  !   !   input: the dimension of the mesh

  ! /
  !  !  output: the rotated point

  !/
  !
  !  local variable
  !
  ! the rotation matrix in scaled crystallographic
  integer :: ss (3, 3)
  ! axes. Compatibility with the FFT grid must have 
  ! been checked elsewhere (sgam_at)
  !
  ! this is a temporary fix. Much better would be to pass directly the ss 
  ! matrix
  !
  ss (1, 1) = s (1, 1)
  ss (2, 1) = s (2, 1) * nr1 / nr2
  ss (3, 1) = s (3, 1) * nr1 / nr3
  ss (1, 2) = s (1, 2) * nr2 / nr1
  ss (2, 2) = s (2, 2)
  ss (3, 2) = s (3, 2) * nr2 / nr3
  ss (1, 3) = s (1, 3) * nr3 / nr1
  ss (2, 3) = s (2, 3) * nr3 / nr2
  ss (3, 3) = s (3, 3)
  !
  ri = ss (1, 1) * (i - 1) + ss (2, 1) * (j - 1) + ss (3, 1) &
       * (k - 1) - ftau (1)
  ri = mod (ri, nr1) + 1
  if (ri.lt.1) ri = ri + nr1
  rj = ss (1, 2) * (i - 1) + ss (2, 2) * (j - 1) + ss (3, 2) &
       * (k - 1) - ftau (2)
  rj = mod (rj, nr2) + 1
  if (rj.lt.1) rj = rj + nr2
  rk = ss (1, 3) * (i - 1) + ss (2, 3) * (j - 1) + ss (3, 3) &
       * (k - 1) - ftau (3)
  rk = mod (rk, nr3) + 1
  if (rk.lt.1) rk = rk + nr3
  return
end subroutine ruotaijk

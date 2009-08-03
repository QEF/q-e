!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE sym_cart_tensor(tens)
  !-----------------------------------------------------------------------
  !
  ! ... symmetrize a rank-2 tensor in cartesian coordinates
  ! 
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE symme,     ONLY : s, nsym

  IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: tens(3,3)

  ! cartesian to crystal 
  call trntns (tens, at, bg, -1)

  ! symmetrize
  call symtns (tens, nsym, s)
 
  ! crystal to cartesian
  call trntns (tens, at, bg, 1)

END SUBROUTINE sym_cart_tensor

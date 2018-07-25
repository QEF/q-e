!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE export_gstart_2_ppcg(gstart_)
  !----------------------------------------------------------------------------
  !
  USE mp_bands_util, ONLY : gstart
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: gstart_
  !
  gstart = gstart_
  !
END SUBROUTINE export_gstart_2_ppcg

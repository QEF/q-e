!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE scale_fft_3drism(gvec, at_old, laue)
  !---------------------------------------------------------------------------
  !
  ! ... scale G-vectors of 3D-RISM, for Variable Cell
  !
  USE cell_base,   ONLY : bg
  USE gvec_3drism, ONLY : gvec_type, gshells_3drism
  USE kinds,       ONLY : DP
  !
  IMPLICIT NONE
  !
  TYPE(gvec_type), INTENT(INOUT) :: gvec
  REAL(DP),        INTENT(IN)    :: at_old(3, 3)
  LOGICAL,         INTENT(IN)    :: laue
  !
  INTEGER  :: ig
  REAL(DP) :: gx
  REAL(DP) :: gy
  REAL(DP) :: gz
  !
  ! ... scale G-vectors
  CALL cryst_to_cart(gvec%ngm, gvec%g, at_old, -1)
  CALL cryst_to_cart(gvec%ngm, gvec%g, bg,     +1)
  !
  ! ... update G^2
  DO ig = 1, gvec%ngm
    gx = gvec%g(1, ig)
    gy = gvec%g(2, ig)
    gz = gvec%g(3, ig)
    gvec%gg(ig) = gx * gx + gy * gy + gz * gz
  END DO
  !
  ! ... update G-shell
  IF (.NOT. laue) THEN
    CALL gshells_3drism(gvec)
  END IF
  !
END SUBROUTINE scale_fft_3drism

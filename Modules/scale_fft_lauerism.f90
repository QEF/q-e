!
! Copyright (C) 2016 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE scale_fft_lauerism(lfft, at_old)
  !---------------------------------------------------------------------------
  !
  ! ... scale Gxy-vectors of Laue-RISM, for Variable Cell
  !
  USE cell_base, ONLY : bg
  USE kinds,     ONLY : DP
  USE lauefft,   ONLY : lauefft_type
  !
  IMPLICIT NONE
  !
  TYPE(lauefft_type), INTENT(INOUT) :: lfft
  REAL(DP),           INTENT(IN)    :: at_old(3, 3)
  !
  INTEGER  :: ig
  REAL(DP) :: gx
  REAL(DP) :: gy
  REAL(DP) :: gg
  !
  ! ... scale Gxy-vectors
  CALL cryst_to_cart_2d(lfft%ngxy, lfft%gxy, at_old, -1)
  CALL cryst_to_cart_2d(lfft%ngxy, lfft%gxy, bg,     +1)
  !
  ! ... update |Gxy|, Gxy^2
  DO ig = 1, lfft%ngxy
    gx = lfft%gxy(1, ig)
    gy = lfft%gxy(2, ig)
    gg = gx * gx + gy * gy
    lfft%gnxy(ig) = SQRT(gg)
    lfft%ggxy(ig) = gg
  END DO
  !
  ! ... update Gxy-shell
  CALL gxyshells(lfft, .FALSE.) ! lmovecell must be .FALSE.
  !
END SUBROUTINE scale_fft_lauerism

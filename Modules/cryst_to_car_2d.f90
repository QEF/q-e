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
SUBROUTINE cryst_to_cart_2d(nvec, vec, trmat, iflag)
  !---------------------------------------------------------------------------
  !
  ! ... transform 2d-components from crystallographic to cartesian coordinates (iflag = +1)
  ! ...                    , or from cartesian to crystallographic coordinates (iflag = -1).
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: nvec
  REAL(DP), INTENT(INOUT) :: vec(2, nvec)
  REAL(DP), INTENT(IN)    :: trmat(3, 3)
  INTEGER,  INTENT(IN)    :: iflag
  !
  INTEGER  :: ivec
  INTEGER  :: ipol
  REAL(DP) :: tau(2)
  !
  DO ivec = 1, nvec
    !
    IF (iflag > 0) then
      ! ... crystallographic -> cartesian
      DO ipol = 1, 2
        tau(ipol) = trmat(ipol, 1) * vec(1, ivec) &
                & + trmat(ipol, 2) * vec(2, ivec)
      END DO
      !
    ELSE IF (iflag < 0) then
      ! ... cartesian -> crystallographic
      DO ipol = 1, 2
        tau(ipol) = trmat(1, ipol) * vec(1, ivec) &
                & + trmat(2, ipol) * vec(2, ivec)
      END DO
      !
    ELSE
      ! ... not convert
      DO ipol = 1, 2
        tau(ipol) = vec(ipol, ivec)
      END DO
    END IF
    !
    DO ipol = 1, 2
      vec(ipol, ivec) = tau(ipol)
    END DO
    !
  END DO
  !
END SUBROUTINE cryst_to_cart_2d

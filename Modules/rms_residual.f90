!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE rms_residual(mres, nres, vres, rms, comm)
  !---------------------------------------------------------------------------
  !
  ! ... calculate RMS of residual vector
  !
  USE kinds, ONLY : DP
  USE mp,    ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: mres
  INTEGER,  INTENT(IN)  :: nres
  REAL(DP), INTENT(IN)  :: vres(1:*)
  REAL(DP), INTENT(OUT) :: rms
  INTEGER,  INTENT(IN)  :: comm
  !
  REAL(DP) :: rsum
  !
  REAL(DP), EXTERNAL :: ddot
  !
  IF (mres < 1) THEN
    rms = 0.0_DP
    RETURN
  END IF
  !
  IF (nres > 0) THEN
    rsum = ddot(nres, vres, 1, vres, 1)
  ELSE
    rsum = 0.0_DP
  END IF
  !
  CALL mp_sum(rsum, comm)
  !
  rms = SQRT(rsum / DBLE(mres))
  !
END SUBROUTINE rms_residual

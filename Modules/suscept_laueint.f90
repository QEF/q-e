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
SUBROUTINE suscept_laueint(rismt, lhand, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... integrate inter-site susceptibility for Laue-RISM
  !
  USE cell_base, ONLY : alat
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,   ONLY : get_nuniq_in_solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  LOGICAL,         INTENT(IN)    :: lhand  ! if true, right-hand. if false, left-hand.
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER  :: nq
  INTEGER  :: iq1
  INTEGER  :: iq2
  INTEGER  :: iiq2
  INTEGER  :: irz
  REAL(DP) :: rz
  REAL(DP) :: rstep
  REAL(DP) :: x
  REAL(DP) :: x0
  REAL(DP) :: x1
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... set dz (in a.u.)
  rstep = alat * rismt%lfft%zstep
  !
  ! ... integrate xgs
  IF (rismt%nrzl * rismt%nsite * rismt%mp_site%nsite > 0) THEN
    IF (lhand) THEN
      rismt%xgs0 = 0.0_DP
      rismt%xgs1 = 0.0_DP
    ELSE
      rismt%ygs0 = 0.0_DP
      rismt%ygs1 = 0.0_DP
    END IF
  END IF
  !
  ! ... calculation only for Gxy = 0
  IF (rismt%lfft%gxystart > 1) THEN
    !NOP
  ELSE
    GOTO 1
  END IF
  !
  DO iq1 = 1, nq
    DO iq2 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      iiq2 = iq2 - rismt%mp_site%isite_start + 1
      !
      x0 = 0.0_DP
      x1 = 0.0_DP
      !
      IF (lhand) THEN
        DO irz = rismt%lfft%nrz, 1, -1
          rz = DBLE(irz - 1) * rstep
          x  = rismt%xgs(irz, iiq2, iq1)
          x0 = x0 + x * rstep
          x1 = x1 + x * rstep * rz
          rismt%xgs0(irz, iiq2, iq1) = x0
          rismt%xgs1(irz, iiq2, iq1) = x1
        END DO
        !
      ELSE
        DO irz = rismt%lfft%nrz, 1, -1
          rz = DBLE(irz - 1) * rstep
          x  = rismt%xgs(irz, iiq2, iq1)
          x0 = x0 + x * rstep
          x1 = x1 + x * rstep * rz
          rismt%ygs0(irz, iiq2, iq1) = x0
          rismt%ygs1(irz, iiq2, iq1) = x1
        END DO
      END IF
      !
    END DO
  END DO
  !
1 CONTINUE
  IF (rismt%nrzl * rismt%nsite * rismt%mp_site%nsite > 0) THEN
    IF (lhand) THEN
      CALL mp_sum(rismt%xgs0, rismt%mp_site%intra_sitg_comm)
      CALL mp_sum(rismt%xgs1, rismt%mp_site%intra_sitg_comm)
    ELSE
      CALL mp_sum(rismt%ygs0, rismt%mp_site%intra_sitg_comm)
      CALL mp_sum(rismt%ygs1, rismt%mp_site%intra_sitg_comm)
    END IF
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE suscept_laueint

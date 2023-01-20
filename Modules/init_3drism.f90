!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE init_3drism(rism3t, solu, lboth, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... initialize 3D-RISM
  !
  USE err_rism,      ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,         ONLY : DP
  USE rism,          ONLY : rism_type, ITYPE_1DRISM, ITYPE_3DRISM, ITYPE_LAUERISM
  USE rism1d_facade, ONLY : rism1t, bond_width, rism1d_activate_right, rism1d_activate_left
  USE solute,        ONLY : update_solU
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rism3t
  LOGICAL,         INTENT(IN)    :: solu   ! create solute's LJ, or not
  LOGICAL,         INTENT(IN)    :: lboth  ! both-hands calculation, or not
  INTEGER,         INTENT(OUT)   :: ierr
  !
  ! ... check data type
  IF (rism1t%itype /= ITYPE_1DRISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rism3t%itype /= ITYPE_3DRISM .AND. rism3t%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... create Lennard-Jones potentials
  IF (solu) THEN
    CALL update_solU(rism3t, ierr)
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      RETURN
    END IF
  END IF
  !
  ! ... create solvent's susceptibilities
  IF (rism3t%itype == ITYPE_3DRISM) THEN
    CALL suscept_vv(rism1t, rism3t, ierr)
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      RETURN
    END IF
    !
  ELSE !IF (rism3t%itype == ITYPE_LAUERISM) THEN
    CALL rism1d_activate_right()
    CALL suscept_laue(rism1t, rism3t, bond_width, .TRUE., ierr)
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      RETURN
    END IF
    !
    IF (lboth) THEN
      CALL rism1d_activate_left()
      CALL suscept_laue(rism1t, rism3t, bond_width, .FALSE., ierr)
      !
      IF (ierr /= IERR_RISM_NULL) THEN
        RETURN
      END IF
    END IF
  END IF
  !
  ! ... integrate solvent's susceptibilities
  IF (rism3t%itype == ITYPE_LAUERISM) THEN
    CALL rism1d_activate_right()
    CALL suscept_laueint(rism3t, .TRUE., ierr)
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      RETURN
    END IF
    !
    IF (lboth) THEN
      CALL rism1d_activate_left()
      CALL suscept_laueint(rism3t, .FALSE., ierr)
      !
      IF (ierr /= IERR_RISM_NULL) THEN
        RETURN
      END IF
    END IF
  END IF
  !
  ! ... prepare dipole part of Laue-RISM equation
  IF (rism3t%itype == ITYPE_LAUERISM) THEN
    CALL eqn_lauedipole(rism3t, .TRUE., .TRUE., ierr)
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      RETURN
    END IF
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE init_3drism

!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE closure(rismt, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... solve a Closure Equation for RISM
  !
  USE constants, ONLY : K_BOLTZMANN_RY
  USE kinds,     ONLY : DP
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE rism,      ONLY : rism_type, ITYPE_1DRISM, ITYPE_LAUERISM, &
                      & CLOSURE_HNC, CLOSURE_KH
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER  :: iclosure
  INTEGER  :: mr, lr
  REAL(DP) :: beta
  !
  ! ... check data type
  IF (rismt%itype == ITYPE_1DRISM .AND. rismt%nr /= rismt%ng) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... set variables
  iclosure = rismt%closure
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  mr = rismt%nr * rismt%nsite
  !
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    lr = 0
  ELSE
    lr = rismt%nrzl * rismt%nsite
  END IF
  !
  ! ... solve closure equation
  IF (iclosure == CLOSURE_HNC) THEN
    !
    IF (rismt%itype /= ITYPE_LAUERISM) THEN
      IF (mr > 0) THEN
        CALL closure_HNC_x(mr, beta, &
           & rismt%usr(1, 1), rismt%hr(1, 1), rismt%csr(1, 1), rismt%gr(1, 1))
      END IF
      !
    ELSE
      IF (mr > 0) THEN
        CALL closure_HNC_x(mr, beta, &
           & rismt%usr(1, 1), rismt%hr(1, 1), rismt%csdr(1, 1), rismt%gr(1, 1))
      END IF
      !
      IF (lr > 0) THEN
        CALL closure_HNC_x(lr, beta, &
           & rismt%usg0(1, 1), rismt%hg0(1, 1), rismt%csdg0(1, 1), rismt%gg0(1, 1))
      END IF
    END IF
    !
  ELSE IF (iclosure == CLOSURE_KH) THEN
    !
    IF (rismt%itype /= ITYPE_LAUERISM) THEN
      IF (mr > 0) THEN
        CALL closure_KH_x(mr, beta, &
           & rismt%usr(1, 1), rismt%hr(1, 1), rismt%csr(1, 1), rismt%gr(1, 1))
      END IF
      !
    ELSE
      IF (mr > 0) THEN
        CALL closure_KH_x(mr, beta, &
           & rismt%usr(1, 1), rismt%hr(1, 1), rismt%csdr(1, 1), rismt%gr(1, 1))
      END IF
      !
      IF (lr > 0) THEN
        CALL closure_KH_x(lr, beta, &
           & rismt%usg0(1, 1), rismt%hg0(1, 1), rismt%csdg0(1, 1), rismt%gg0(1, 1))
      END IF
    END IF
    !
  ELSE
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... set zero, if R = 0 (for 1D-RISM)
  IF (rismt%itype == ITYPE_1DRISM) THEN
    IF (rismt%mp_task%ivec_start == 1 .AND. rismt%nsite > 0) THEN
      rismt%gr(1, :) = 0.0_DP
    END IF
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE closure
!
!---------------------------------------------------------------------------
SUBROUTINE closure_HNC_x(nr, beta, ur, hr, cr, gr)
  !---------------------------------------------------------------------------
  !
  ! ... HyperNetted-Chain model
  ! ... (J.P.Hansen et al., Theory of simple liquids. Academic Press, London, 1990.)
  ! ...
  ! ...   g(r) = exp(t(r))
  ! ...   t(r) = -beta * u(r) + h(r) - c(r)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nr
  REAL(DP), INTENT(IN)  :: beta
  REAL(DP), INTENT(IN)  :: ur(1:*)
  REAL(DP), INTENT(IN)  :: hr(1:*)
  REAL(DP), INTENT(IN)  :: cr(1:*)
  REAL(DP), INTENT(OUT) :: gr(1:*)
  !
  INTEGER  :: ir
  REAL(DP) :: tr
  !
  REAL(DP), PARAMETER :: MAX_EXP = 100.0_DP
  !
!$omp parallel do default(shared) private(ir, tr)
  DO ir = 1, nr
    tr = -beta * ur(ir) + hr(ir) - cr(ir)
    gr(ir) = EXP(MIN(tr, MAX_EXP))
  END DO
!$omp end parallel do
  !
END SUBROUTINE closure_HNC_x
!
!---------------------------------------------------------------------------
SUBROUTINE closure_KH_x(nr, beta, ur, hr, cr, gr)
  !---------------------------------------------------------------------------
  !
  ! ... Kovalenko and Hirata's model
  ! ... (A.Kovalenko, F.Hirata, J. Chem. Phys. 1999, 110, 10095-10112)
  ! ...
  ! ...   g(r) = exp(t(r)), if t(r) < 0
  ! ...   g(r) = 1 + t(r) , if t(r) > 0
  ! ...   t(r) = -beta * u(r) + h(r) - c(r)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nr
  REAL(DP), INTENT(IN)  :: beta
  REAL(DP), INTENT(IN)  :: ur(1:*)
  REAL(DP), INTENT(IN)  :: hr(1:*)
  REAL(DP), INTENT(IN)  :: cr(1:*)
  REAL(DP), INTENT(OUT) :: gr(1:*)
  !
  INTEGER  :: ir
  REAL(DP) :: tr
  !
!$omp parallel do default(shared) private(ir, tr)
  DO ir = 1, nr
    tr = -beta * ur(ir) + hr(ir) - cr(ir)
    IF (tr < 0.0_DP) THEN
      gr(ir) = EXP(tr)
    ELSE
      gr(ir) = 1.0_DP + tr
    END IF
  END DO
!$omp end parallel do
  !
END SUBROUTINE closure_KH_x


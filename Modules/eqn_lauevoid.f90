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
SUBROUTINE eqn_lauevoid(rismt, lboth, expand, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... solve Laue-RISM equation from void-region, which is defined as
  ! ...
  ! ...                  /
  ! ...   h1(gxy=0,z1) = | dz2 c2(gxy=0,z2) * x21(gxy=0,z2-z1)
  ! ...                  /void-region
  ! ...
  ! ... void-region is in right-hand side, left-hand side or between right- and left-hand side,
  ! ... where solvents does not exist and c2 is linear function.
  ! ... calculated total correlations are added to `hgz' or `hsgz'.
  ! ...
  !
  USE err_rism, ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE rism,     ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,  ONLY : get_nuniq_in_solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  LOGICAL,         INTENT(IN)    :: lboth  ! both-hands calculation, or not
  LOGICAL,         INTENT(IN)    :: expand ! expand-cell(.TRUE.) or unit-cell(.FALSE.)
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER :: nq
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
  IF (rismt%nrzs < rismt%dfft%nr3) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%lfft%xright .AND. rismt%lfft%xleft) THEN
    !
    ! ... void-region is between right- and left-hand side
    CALL eqn_lauevoid_between(rismt, lboth, expand)
    !
  ELSE
    !
    ! ... void-region is right-hand side or left-hand side
    CALL eqn_lauevoid_oneside(rismt, expand)
    !
  END IF
  !
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE eqn_lauevoid
!
!---------------------------------------------------------------------------
SUBROUTINE eqn_lauevoid_oneside(rismt, expand)
  !---------------------------------------------------------------------------
  !
  ! ... Laue-RISM equation from void-region of right-hand side or left-hand side.
  !
  USE cell_base, ONLY : alat
  USE constants, ONLY : K_BOLTZMANN_RY
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type
  USE solvmol,   ONLY : solVs, get_nuniq_in_solVs, &
                      & iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  LOGICAL,         INTENT(IN)    :: expand ! expand-cell(.TRUE.) or unit-cell(.FALSE.)
  !
  INTEGER               :: nq
  INTEGER               :: iq1, iq2
  INTEGER               :: iiq1, iiq2
  INTEGER               :: iv2
  INTEGER               :: isolV2
  INTEGER               :: iatom2
  INTEGER               :: iz
  INTEGER               :: izsta
  INTEGER               :: izend
  INTEGER               :: izsolv
  INTEGER               :: izvoid
  INTEGER               :: nzint
  INTEGER               :: izint
  INTEGER               :: izdelt
  REAL(DP)              :: beta
  REAL(DP)              :: qv2
  REAL(DP)              :: z
  REAL(DP)              :: zstep
  REAL(DP)              :: zoffs
  REAL(DP)              :: zedge
  REAL(DP)              :: voppo
  REAL(DP)              :: vsign
  REAL(DP)              :: cz
  REAL(DP)              :: dz
  REAL(DP), ALLOCATABLE :: c2(:)
  REAL(DP), ALLOCATABLE :: d2(:)
  REAL(DP), ALLOCATABLE :: h1(:)
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... beta = 1 / (kB * T)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... set integral regions as index of long Z-stick (i.e. expanded cell)
  IF (rismt%lfft%xright) THEN
    IF (expand) THEN
      izsta = rismt%lfft%izright_gedge
      izend = rismt%lfft%nrz
    ELSE
      izsta = rismt%lfft%izright_start0
      izend = rismt%lfft%izright_end0
    END IF
    !
    izsolv = rismt%lfft%izright_start0
    izvoid = izsolv - 1
    !
    IF (rismt%lfft%gxystart > 1) THEN
      voppo = DBLE(rismt%vleft(1)) / alat
    ELSE
      voppo = 0.0_DP
    END IF
    !
    vsign = -1.0_DP
    !
  ELSE !IF (rismt%lfft%xleft) THEN
    IF (expand) THEN
      izsta = 1
      izend = rismt%lfft%izleft_gedge
    ELSE
      izsta = rismt%lfft%izleft_start0
      izend = rismt%lfft%izleft_end0
    END IF
    !
    izsolv = rismt%lfft%izleft_end0
    izvoid = izsolv + 1
    !
    IF (rismt%lfft%gxystart > 1) THEN
      voppo = DBLE(rismt%vright(1)) / alat
    ELSE
      voppo = 0.0_DP
    END IF
    !
    vsign = +1.0_DP
  END IF
  !
  ! ... count integral points along Z
  nzint = izend - izsta + 1
  !
  ! ... properties about length (in a.u.)
  zstep = alat * rismt%lfft%zstep
  zoffs = alat * (rismt%lfft%zleft + rismt%lfft%zoffset)
  zedge = zoffs + zstep * DBLE(izsolv - 1)
  !
  ! ... allocate working memory
  IF (rismt%nsite > 0) THEN
    ALLOCATE(c2(rismt%nsite))
    ALLOCATE(d2(rismt%nsite))
  END IF
  IF (nzint > 0) THEN
    ALLOCATE(h1(nzint))
  END IF
  !
  ! ... calculate c2, d2
  DO iq2 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq2   = iq2 - rismt%mp_site%isite_start + 1
    iv2    = iuniq_to_isite(1, iq2)
    isolV2 = isite_to_isolV(iv2)
    iatom2 = isite_to_iatom(iv2)
    qv2    = solVs(isolV2)%charge(iatom2)
    !
    IF (rismt%lfft%gxystart > 1) THEN
      c2(iiq2) = rismt%csdg0(izsolv, iiq2) &
             & - beta * qv2 * DBLE(rismt%vlgz(izsolv))
      d2(iiq2) = -beta * qv2 * voppo
    ELSE
      c2(iiq2) = 0.0_DP
      d2(iiq2) = 0.0_DP
    END IF
  END DO
  !
  IF (rismt%nsite > 0) THEN
    CALL mp_sum(c2, rismt%mp_site%intra_sitg_comm)
    CALL mp_sum(d2, rismt%mp_site%intra_sitg_comm)
  END IF
  !
  ! ... Laue-RISM equation of void-region
  DO iq1 = 1, nq
    IF (rismt%mp_site%isite_start <= iq1 .AND. iq1 <= rismt%mp_site%isite_end) THEN
      iiq1 = iq1 - rismt%mp_site%isite_start + 1
    ELSE
      iiq1 = 0
    END IF
    !
    IF (nzint > 0) THEN
      h1 = 0.0_DP
    END IF
    !
    DO iq2 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      iiq2 = iq2 - rismt%mp_site%isite_start + 1
      !
      ! ... h1(z1)
!$omp parallel do default(shared) private(iz, izint, izdelt, z, cz, dz)
      DO iz = izsta, izend
        izint  = iz - izsta + 1
        izdelt = ABS(iz - izvoid) + 1
        !
        IF (izdelt <= rismt%lfft%nrz) THEN
          z  = zoffs + zstep * DBLE(iz - 1)
          cz = c2(iiq2) + d2(iiq2) * (z - zedge)
          dz = d2(iiq2) * vsign
          h1(izint) = h1(izint) &
          & + cz * rismt%xgs0(izdelt, iiq2, iq1) &
          & + dz * rismt%xgs1(izdelt, iiq2, iq1)
        END IF
      END DO
!$omp end parallel do
    END DO
    !
    IF (nzint > 0) THEN
      CALL mp_sum(h1, rismt%mp_site%inter_sitg_comm)
    END IF
    !
    IF (iiq1 > 0) THEN
      IF (expand) THEN
        ! ... add h1 -> hsgz
        IF (rismt%lfft%gxystart > 1) THEN
!$omp parallel do default(shared) private(iz, izint)
          DO iz = izsta, izend
            izint = iz - izsta + 1
            rismt%hsgz(iz, iiq1) = rismt%hsgz(iz, iiq1) + CMPLX(h1(izint), 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
        END IF
        !
      ELSE
        ! ... add h1 -> hg0
!$omp parallel do default(shared) private(iz, izint)
        DO iz = izsta, izend
          izint = iz - izsta + 1
          rismt%hg0(iz, iiq1) = rismt%hg0(iz, iiq1) + h1(izint)
        END DO
!$omp end parallel do
      END IF
    END IF
    !
  END DO
  !
  ! ... deallocate working memory
  IF (rismt%nsite > 0) THEN
    DEALLOCATE(c2)
    DEALLOCATE(d2)
  END IF
  IF (nzint > 0) THEN
    DEALLOCATE(h1)
  END IF
  !
END SUBROUTINE eqn_lauevoid_oneside
!
!---------------------------------------------------------------------------
SUBROUTINE eqn_lauevoid_between(rismt, lboth, expand)
  !---------------------------------------------------------------------------
  !
  ! ... Laue-RISM equation from void-region between right- and left-hand side.
  !
  USE cell_base, ONLY : alat
  USE constants, ONLY : K_BOLTZMANN_RY
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type
  USE solvmol,   ONLY : solVs, get_nuniq_in_solVs, &
                      & iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  LOGICAL,         INTENT(IN)    :: lboth  ! both-hands calculation, or not
  LOGICAL,         INTENT(IN)    :: expand ! expand-cell(.TRUE.) or unit-cell(.FALSE.)
  !
  INTEGER               :: nq
  INTEGER               :: iq1, iq2
  INTEGER               :: iiq1, iiq2
  INTEGER               :: iv2
  INTEGER               :: isolV2
  INTEGER               :: iatom2
  INTEGER               :: iz
  INTEGER               :: nzright
  INTEGER               :: izright_sta
  INTEGER               :: izright_end
  INTEGER               :: izright_solv
  INTEGER               :: izright_void
  INTEGER               :: nzleft
  INTEGER               :: izleft_sta
  INTEGER               :: izleft_end
  INTEGER               :: izleft_solv
  INTEGER               :: izleft_void
  INTEGER               :: izint
  INTEGER               :: izdelt1
  INTEGER               :: izdelt2
  REAL(DP)              :: beta
  REAL(DP)              :: qv2
  REAL(DP)              :: z
  REAL(DP)              :: zstep
  REAL(DP)              :: zoffs
  REAL(DP)              :: zright_edge
  REAL(DP)              :: zleft_edge
  REAL(DP)              :: c2, cz
  REAL(DP)              :: d2, dz
  REAL(DP), ALLOCATABLE :: xg0(:)
  REAL(DP), ALLOCATABLE :: xg1(:)
  REAL(DP), ALLOCATABLE :: cright(:)
  REAL(DP), ALLOCATABLE :: cleft(:)
  REAL(DP), ALLOCATABLE :: hright(:)
  REAL(DP), ALLOCATABLE :: hleft(:)
  !
  ! ... has void-region ?
  IF ((rismt%lfft%izright_start0 - rismt%lfft%izleft_end0) <= 1) THEN
    RETURN
  END IF
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... beta = 1 / (kB * T)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... set integral regions as index of long Z-stick (i.e. expanded cell)
  IF (expand) THEN
    izright_sta = rismt%lfft%izright_gedge
    izright_end = rismt%lfft%nrz
    izleft_sta  = 1
    izleft_end  = rismt%lfft%izleft_gedge
  ELSE
    izright_sta = rismt%lfft%izright_start0
    izright_end = rismt%lfft%izright_end0
    izleft_sta  = rismt%lfft%izleft_start0
    izleft_end  = rismt%lfft%izleft_end0
  END IF
  !
  izright_solv = rismt%lfft%izright_start0
  izright_void = izright_solv - 1
  izleft_solv  = rismt%lfft%izleft_end0
  izleft_void  = izleft_solv  + 1
  !
  ! ... count integral points along Z
  nzright = izright_end - izright_sta + 1
  nzleft  = izleft_end  - izleft_sta  + 1
  !
  ! ... properties about length (in a.u.)
  zstep = alat * rismt%lfft%zstep
  zoffs = alat * (rismt%lfft%zleft + rismt%lfft%zoffset)
  !
  zright_edge = zoffs + zstep * DBLE(izright_solv - 1)
  zleft_edge  = zoffs + zstep * DBLE(izleft_solv  - 1)
  !
  ! ... allocate working memory
  IF (rismt%nrzl > 0) THEN
    ALLOCATE(xg0(rismt%nrzl))
    ALLOCATE(xg1(rismt%nrzl))
  END IF
  IF (rismt%nsite > 0) THEN
    ALLOCATE(cright(rismt%nsite))
    ALLOCATE(cleft( rismt%nsite))
  END IF
  IF (nzright > 0) THEN
    ALLOCATE(hright(nzright))
  END IF
  IF (nzleft > 0) THEN
    ALLOCATE(hleft(nzleft))
  END IF
  !
  ! ... calculate cright(z2), cleft(z2)
  DO iq2 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq2   = iq2 - rismt%mp_site%isite_start + 1
    iv2    = iuniq_to_isite(1, iq2)
    isolV2 = isite_to_isolV(iv2)
    iatom2 = isite_to_iatom(iv2)
    qv2    = solVs(isolV2)%charge(iatom2)
    !
    IF (rismt%lfft%gxystart > 1) THEN
      cright(iiq2) = rismt%csdg0(izright_solv, iiq2) &
                   & - beta * qv2 * DBLE(rismt%vlgz(izright_solv))
      cleft( iiq2) = rismt%csdg0(izleft_solv, iiq2) &
                   & - beta * qv2 * DBLE(rismt%vlgz(izleft_solv))
    ELSE
      cright(iiq2) = 0.0_DP
      cleft( iiq2) = 0.0_DP
    END IF
  END DO
  !
  IF (rismt%nsite > 0) THEN
    CALL mp_sum(cright, rismt%mp_site%intra_sitg_comm)
    CALL mp_sum(cleft,  rismt%mp_site%intra_sitg_comm)
  END IF
  !
  ! ... Laue-RISM equation of void-region
  DO iq1 = 1, nq
    IF (rismt%mp_site%isite_start <= iq1 .AND. iq1 <= rismt%mp_site%isite_end) THEN
      iiq1 = iq1 - rismt%mp_site%isite_start + 1
    ELSE
      iiq1 = 0
    END IF
    !
    IF (nzright > 0) THEN
      hright = 0.0_DP
    END IF
    IF (nzleft > 0) THEN
      hleft = 0.0_DP
    END IF
    !
    DO iq2 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      iiq2 = iq2 - rismt%mp_site%isite_start + 1
      !
      d2 = (cright(iiq2) - cleft(iiq2)) / (zright_edge - zleft_edge)
      !
      ! ... hleft(z1)
      c2 = cleft(iiq2)
      !
      IF (.NOT. lboth) THEN
        xg0 = rismt%xgs0(1:rismt%nrzl, iiq2, iq1)
        xg1 = rismt%xgs1(1:rismt%nrzl, iiq2, iq1)
      ELSE
        xg0 = rismt%ygs0(1:rismt%nrzl, iiq2, iq1)
        xg1 = rismt%ygs1(1:rismt%nrzl, iiq2, iq1)
      END IF
      !
!$omp parallel do default(shared) private(iz, izint, izdelt1, izdelt2, z, cz, dz)
      DO iz = izleft_sta, izleft_end
        izint   = iz - izleft_sta + 1
        izdelt1 = ABS(iz - izleft_void ) + 1
        izdelt2 = ABS(iz - izright_solv) + 1
        !
        IF (izdelt1 <= rismt%lfft%nrz) THEN
          z  = zoffs + zstep * DBLE(iz - 1)
          cz = c2 + d2 * (z - zleft_edge)
          dz = d2
          hleft(izint) = hleft(izint) &
          & + cz * xg0(izdelt1) &
          & + dz * xg1(izdelt1)
        END IF
        !
        IF (izdelt2 <= rismt%lfft%nrz) THEN
          z  = zoffs + zstep * DBLE(iz - 1)
          cz = c2 + d2 * (z - zleft_edge)
          dz = d2
          hleft(izint) = hleft(izint) &
          & - cz * xg0(izdelt2) &
          & - dz * xg1(izdelt2)
        END IF
      END DO
!$omp end parallel do
      !
      ! ... hright(z1)
      c2 = cright(iiq2)
      !
      xg0 = rismt%xgs0(1:rismt%nrzl, iiq2, iq1)
      xg1 = rismt%xgs1(1:rismt%nrzl, iiq2, iq1)
      !
!$omp parallel do default(shared) private(iz, izint, izdelt1, izdelt2, z, cz, dz)
      DO iz = izright_sta, izright_end
        izint   = iz - izright_sta + 1
        izdelt1 = ABS(iz - izright_void) + 1
        izdelt2 = ABS(iz - izleft_solv ) + 1
        !
        IF (izdelt1 <= rismt%lfft%nrz) THEN
          z  = zoffs + zstep * DBLE(iz - 1)
          cz = c2 + d2 * (z - zright_edge)
          dz = -d2
          hright(izint) = hright(izint) &
          & + cz * xg0(izdelt1) &
          & + dz * xg1(izdelt1)
        END IF
        !
        IF (izdelt2 <= rismt%lfft%nrz) THEN
          z  = zoffs + zstep * DBLE(iz - 1)
          cz = c2 + d2 * (z - zright_edge)
          dz = -d2
          hright(izint) = hright(izint) &
          & - cz * xg0(izdelt2) &
          & - dz * xg1(izdelt2)
        END IF
      END DO
!$omp end parallel do
    END DO
    !
    IF (nzright > 0) THEN
      CALL mp_sum(hright, rismt%mp_site%inter_sitg_comm)
    END IF
    IF (nzleft > 0) THEN
      CALL mp_sum(hleft, rismt%mp_site%inter_sitg_comm)
    END IF
    !
    IF (iiq1 > 0) THEN
      IF (expand) THEN
        IF (rismt%lfft%gxystart > 1) THEN
          ! ... add hleft -> hsgz
!$omp parallel do default(shared) private(iz, izint)
          DO iz = izleft_sta, izleft_end
            izint = iz - izleft_sta + 1
            rismt%hsgz(iz, iiq1) = rismt%hsgz(iz, iiq1) + CMPLX(hleft(izint), 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
          ! ... add hright -> hsgz
!$omp parallel do default(shared) private(iz, izint)
          DO iz = izright_sta, izright_end
            izint = iz - izright_sta + 1
            rismt%hsgz(iz, iiq1) = rismt%hsgz(iz, iiq1) + CMPLX(hright(izint), 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
        END IF
        !
      ELSE
        ! ... add hleft -> hg0
!$omp parallel do default(shared) private(iz, izint)
        DO iz = izleft_sta, izleft_end
          izint = iz - izleft_sta + 1
          rismt%hg0(iz, iiq1) = rismt%hg0(iz, iiq1) + hleft(izint)
        END DO
!$omp end parallel do
        !
        ! ... add hright -> hg0
!$omp parallel do default(shared) private(iz, izint)
        DO iz = izright_sta, izright_end
          izint = iz - izright_sta + 1
          rismt%hg0(iz, iiq1) = rismt%hg0(iz, iiq1) + hright(izint)
        END DO
!$omp end parallel do
      END IF
    END IF
    !
  END DO
  !
  ! ... deallocate working memory
  IF (rismt%nrzl > 0) THEN
    DEALLOCATE(xg0)
    DEALLOCATE(xg1)
  END IF
  IF (rismt%nsite > 0) THEN
    DEALLOCATE(cright)
    DEALLOCATE(cleft)
  END IF
  IF (nzright > 0) THEN
    DEALLOCATE(hright)
  END IF
  IF (nzleft > 0) THEN
    DEALLOCATE(hleft)
  END IF
  !
END SUBROUTINE eqn_lauevoid_between

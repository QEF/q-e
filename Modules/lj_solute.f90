!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define RSMIN__  1.0E-6_DP
!
!--------------------------------------------------------------------------
SUBROUTINE lj_setup_solU_tau(rismt, rsmax, count_only, ierr)
  !--------------------------------------------------------------------------
  !
  ! ... setup coordinate of solute's atoms,
  ! ... which can contribute to Lennard-Jones potentials
  !
  USE cell_base, ONLY : at, bg, alat
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE ions_base, ONLY : nat, tau
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  USE solute,    ONLY : solU_nat, solU_tau, solU_ljsig, isup_to_iuni
  USE solvmol,   ONLY : nsolV, solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(IN)  :: rsmax
  LOGICAL,         INTENT(IN)  :: count_only
  INTEGER,         INTENT(OUT) :: ierr
  !
  LOGICAL               :: laue
  INTEGER               :: isolV
  INTEGER               :: iatom
  INTEGER               :: ia
  INTEGER               :: im1, im2, im3
  INTEGER               :: nm1, nm2, nm3
  REAL(DP)              :: rmax
  REAL(DP)              :: sv, su, suv
  REAL(DP)              :: rm1, rm2, rm3
  REAL(DP)              :: bgnrm(3)
  REAL(DP)              :: tau_tmp(3)
  REAL(DP), ALLOCATABLE :: tau_uni(:,:)
  !
  REAL(DP), EXTERNAL :: dnrm2
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... alloc memory
  ALLOCATE(tau_uni(3, nat))
  !
  ! ... set variables
  laue = .FALSE.
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    laue = .TRUE.
  END IF
  !
  bgnrm(1) = dnrm2(3, bg (1, 1), 1)
  bgnrm(2) = dnrm2(3, bg (1, 2), 1)
  bgnrm(3) = dnrm2(3, bg (1, 3), 1)
  !
  sv = 0.0_DP
  DO isolV = 1, nsolV
    DO iatom = 1, solVs(isolV)%natom
      sv = MAX(sv, solVs(isolV)%ljsig(iatom))
    END DO
  END DO
  !
  ! ... count maximum of cell size
  su = 0.0_DP
  DO ia = 1, nat
    su = MAX(su, solU_ljsig(ia))
  END DO
  !
  suv  = 0.5_DP * (sv + su)
  rmax = rsmax * suv / alat
  !
  nm1 = CEILING(bgnrm(1) * rmax)
  nm2 = CEILING(bgnrm(2) * rmax)
  IF (.NOT. laue) THEN
    nm3 = CEILING(bgnrm(3) * rmax)
  ELSE
    nm3 = 0
  END IF
  !
  ! ... wrap coordinates back into cell
  tau_uni = tau
  CALL cryst_to_cart(nat, tau_uni, bg, -1)
  IF (.NOT. laue) THEN
    tau_uni = tau_uni - FLOOR(tau_uni)
  ELSE
    tau_uni(1:2, :) = tau_uni(1:2, :) - FLOOR(tau_uni(1:2, :))
  END IF
  !
  ! ... set unit cell
  solU_nat = nat
  IF (.NOT. count_only) THEN
    DO ia = 1, nat
      solU_tau(:,  ia) = tau_uni(:, ia)
      isup_to_iuni(ia) = ia
    END DO
  END IF
  !
  ! ... set neighbor cells
  DO im1 = -nm1, nm1
    DO im2 = -nm2, nm2
      DO im3 = -nm3, nm3
        !
        IF (im1 == 0 .AND. im2 == 0 .AND. im3 == 0) THEN
          CYCLE
        END IF
        !
        DO ia = 1, nat
          su   = solU_ljsig(ia)
          suv  = 0.5_DP * (sv + su)
          rmax = rsmax * suv / alat
          rm1  = bgnrm(1) * rmax
          rm2  = bgnrm(2) * rmax
          rm3  = bgnrm(3) * rmax
          !
          tau_tmp(1) = tau_uni(1, ia) + DBLE(im1)
          tau_tmp(2) = tau_uni(2, ia) + DBLE(im2)
          tau_tmp(3) = tau_uni(3, ia) + DBLE(im3)
          !
          IF (tau_tmp(1) < -rm1 .OR. (rm1 + 1.0_DP) < tau_tmp(1)) THEN
            CYCLE
          END IF
          !
          IF (tau_tmp(2) < -rm2 .OR. (rm2 + 1.0_DP) < tau_tmp(2)) THEN
            CYCLE
          END IF
          !
          IF (.NOT. laue) THEN
            IF (tau_tmp(3) < -rm3 .OR. (rm3 + 1.0_DP) < tau_tmp(3)) THEN
              CYCLE
            END IF
          END IF
          !
          solU_nat = solU_nat + 1
          IF (.NOT. count_only) THEN
            solU_tau(:,  solU_nat) = tau_tmp(:)
            isup_to_iuni(solU_nat) = ia
          END IF
        END DO
        !
      END DO
    END DO
  END DO
  !
  IF (.NOT. count_only) THEN
    CALL cryst_to_cart(solU_nat, solU_tau, at, +1)
  END IF
  !
  ! ... dealloc memory
  DEALLOCATE(tau_uni)
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE lj_setup_solU_tau
!
!--------------------------------------------------------------------------
SUBROUTINE lj_setup_solU_vlj(rismt, rsmax, ierr)
  !--------------------------------------------------------------------------
  !
  ! ... calculate solute-solvent's Lennard-Jones potential
  ! ...
  ! ...   ----           [ [  sig  ]12    [  sig  ]6 ]
  ! ...   >    4 * esp * [ [-------]   -  [-------]  ]
  ! ...   ----           [ [|r - R|]      [|r - R|]  ]
  ! ...    R
  !
  USE err_rism, ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,    ONLY : DP
  USE rism,     ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  USE solvmol,  ONLY : get_nuniq_in_solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: rsmax
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER :: nq
  INTEGER :: iq
  LOGICAL :: laue
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%dfft%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... Laue-RISM or not
  laue = .FALSE.
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    laue = .TRUE.
  END IF
  !
  ! ... calculate Lennard-Jones
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    CALL lj_setup_solU_vlj_x(iq, rismt, rsmax, laue)
  END DO
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE lj_setup_solU_vlj
!
!--------------------------------------------------------------------------
SUBROUTINE lj_setup_solU_vlj_x(iq, rismt, rsmax, laue)
  !--------------------------------------------------------------------------
  !
  ! ... calculate solute-solvent's Lennard-Jones potential
  ! ... for a solvent's site
  !
  USE cell_base, ONLY : alat, at
  USE fft_types, ONLY : fft_index_to_3d
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type
  USE solute,    ONLY : solU_nat, solU_tau, solU_ljeps, solU_ljsig, isup_to_iuni
  USE solvmol,   ONLY : iuniq_to_isite, isite_to_isolV, isite_to_iatom, solVs
  !
  IMPLICIT NONE
  !
  INTEGER,         INTENT(IN)    :: iq
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: rsmax
  LOGICAL,         INTENT(IN)    :: laue
  !
  REAL(DP), PARAMETER :: RSMIN = RSMIN__
  !
  INTEGER  :: iiq
  INTEGER  :: iv
  INTEGER  :: isolV
  INTEGER  :: iatom
  INTEGER  :: ir, nr, mr
  INTEGER  :: i1, i2, i3
  INTEGER  :: n1, n2, n3
  INTEGER  :: ia, iia
  LOGICAL  :: offrange
  REAL(DP) :: tau_r(3)
  REAL(DP) :: r1, r2, r3
  REAL(DP) :: r3offset
  REAL(DP) :: ev, eu, euv
  REAL(DP) :: sv, su, suv
  REAL(DP) :: rmax
  REAL(DP) :: rmin
  REAL(DP) :: ruv2
  REAL(DP) :: xuv, yuv, zuv
  REAL(DP) :: sr2, sr6, sr12
  REAL(DP) :: vlj
  !
  ! ... FFT box
  n1 = rismt%dfft%nr1
  n2 = rismt%dfft%nr2
  n3 = rismt%dfft%nr3
  nr = rismt%dfft%nr1x * rismt%dfft%my_nr2p * rismt%dfft%my_nr3p
  mr = rismt%dfft%nnr
  !
  ! ... solvent properties
  iiq   = iq - rismt%mp_site%isite_start + 1
  iv    = iuniq_to_isite(1, iq)
  isolV = isite_to_isolV(iv)
  iatom = isite_to_iatom(iv)
  sv    = solVs(isolV)%ljsig(iatom)
  ev    = solVs(isolV)%ljeps(iatom)
  !
  ! ... offset for Laue-RISM
  IF (laue) THEN
#if defined (__ESM_NOT_SYMMETRIC)
    r3offset = 0.0_DP
#else
    IF (MOD(n3, 2) == 0) THEN
      r3offset = 0.5_DP / DBLE(n3)
    ELSE
      r3offset = 0.0_DP
    END IF
#endif
  END IF
  !
  ! ... calculate potential on each FFT grid
!$omp parallel do default(shared) private(ir, i1, i2, i3, offrange, r1, r2, r3, tau_r, vlj, &
!$omp             ia, iia, su, suv, rmax, rmin, xuv, yuv, zuv, ruv2, eu, euv, sr2, sr6, sr12)
  DO ir = 1, mr
    !
    ! ... create coordinate of a FFT grid
    IF (ir <= nr) THEN
      CALL fft_index_to_3d(ir, rismt%dfft, i1, i2, i3, offrange)
    ELSE
      offrange = .TRUE.
    END IF
    !
    IF (offrange) THEN
      rismt%uljr(ir, iiq) = 0.0_DP
      CYCLE
    END IF
    !
    r1 = DBLE(i1) / DBLE(n1)
    r2 = DBLE(i2) / DBLE(n2)
    r3 = DBLE(i3) / DBLE(n3)
    !
    IF (laue) THEN
      r3 = r3 + r3offset
      IF (i3 >= (n3 - (n3 / 2))) THEN
        r3 =  r3 - 1.0_DP
      END IF
    END IF
    !
    tau_r(:) = r1 * at(:, 1) + r2 * at(:, 2) + r3 * at(:, 3)
    !
    ! ... contribution from each solute's atom
    vlj = 0.0_DP
    !
    DO ia = 1, solU_nat
      iia  = isup_to_iuni(ia)
      su   = solU_ljsig(iia)
      suv  = 0.5_DP * (sv + su)
      rmax = rsmax * suv / alat
      rmin = RSMIN * suv / alat
      xuv  = tau_r(1) - solU_tau(1, ia)
      yuv  = tau_r(2) - solU_tau(2, ia)
      zuv  = tau_r(3) - solU_tau(3, ia)
      ruv2 = xuv * xuv + yuv * yuv + zuv * zuv
      IF (ruv2 > (rmax * rmax)) THEN
        CYCLE
      END IF
      IF (ruv2 < (rmin * rmin)) THEN
        ruv2 = rmin * rmin
      END IF
      !
      eu   = solU_ljeps(iia)
      euv  = SQRT(ev * eu)
      sr2  = suv * suv / ruv2 / alat / alat
      sr6  = sr2 * sr2 * sr2
      sr12 = sr6 * sr6
      vlj  = vlj + 4.0_DP * euv * (sr12 - sr6)
    END DO
    !
    rismt%uljr(ir, iiq) = vlj
    !
  END DO
!$omp end parallel do
  !
END SUBROUTINE lj_setup_solU_vlj_x
!
!--------------------------------------------------------------------------
SUBROUTINE lj_setup_wall(rismt, rsmax, ierr)
  !--------------------------------------------------------------------------
  !
  ! ... calculate wall-solvent's Lennard-Jones repulsive potential
  ! ...
  ! ...                                  sig^12
  ! ...   (+-) 2pi * rho * 4 * esp * --------------
  ! ...                              90 * |z - Z|^9
  ! ...
  ! ... optionally, attractive potential is added
  ! ...
  ! ...                                  sig^6
  ! ...   (-+) 2pi * rho * 4 * esp * --------------
  ! ...                              12 * |z - Z|^3
  !
  USE err_rism, ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,    ONLY : DP
  USE rism,     ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,  ONLY : get_nuniq_in_solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: rsmax
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER :: nq
  INTEGER :: iq
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
  IF (rismt%nr < rismt%dfft%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... calculate Lennard-Jones wall
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    CALL lj_setup_wall_x(iq, rismt, rsmax)
  END DO
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE lj_setup_wall
!
!--------------------------------------------------------------------------
SUBROUTINE lj_setup_wall_x(iq, rismt, rsmax)
  !--------------------------------------------------------------------------
  !
  ! ... calculate wall-solvent's Lennard-Jones repulsive potential
  ! ... for a solvent's site
  !
  USE cell_base, ONLY : alat, at
  USE constants, ONLY : tpi
  USE fft_types, ONLY : fft_index_to_3d
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type
  USE solute,    ONLY : iwall, wall_tau, wall_rho, wall_ljeps, wall_ljsig, &
                        wall_lj6, IWALL_RIGHT, IWALL_LEFT
  USE solvmol,   ONLY : iuniq_to_isite, isite_to_isolV, isite_to_iatom, solVs
  !
  IMPLICIT NONE
  !
  INTEGER,         INTENT(IN)    :: iq
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: rsmax
  !
  REAL(DP), PARAMETER :: RSMIN = RSMIN__
  !
  INTEGER  :: iiq
  INTEGER  :: iv
  INTEGER  :: isolV
  INTEGER  :: iatom
  INTEGER  :: ir, nr, mr
  INTEGER  :: i1, i2, i3
  INTEGER  :: n1, n2, n3
  LOGICAL  :: offrange
  REAL(DP) :: tau_z
  REAL(DP) :: r3
  REAL(DP) :: r3offset
  REAL(DP) :: rho
  REAL(DP) :: ev, eu, euv
  REAL(DP) :: sv, su, suv
  REAL(DP) :: rmax
  REAL(DP) :: rmin
  REAL(DP) :: zuv
  REAL(DP) :: sr, sr2, sr3, sr6, sr9
  REAL(DP) :: rsign
  REAL(DP) :: vw
  !
  ! ... which type of wall ?
  IF (iwall == IWALL_RIGHT) THEN
    rsign = -1.0_DP
    !
  ELSE IF (iwall == IWALL_LEFT) THEN
    rsign = +1.0_DP
    !
  ELSE !IF (iwall == IWALL_NULL) THEN
    IF (rismt%dfft%nnr > 0) THEN
      iiq = iq - rismt%mp_site%isite_start + 1
      rismt%uwr(1:rismt%dfft%nnr, iiq) = 0.0_DP
    END IF
    RETURN
  END IF
  !
  ! ... FFT box
  n1  = rismt%dfft%nr1
  n2  = rismt%dfft%nr2
  n3  = rismt%dfft%nr3
  nr = rismt%dfft%nr1x * rismt%dfft%my_nr2p * rismt%dfft%my_nr3p
  mr = rismt%dfft%nnr
  !
  ! ... solvent properties
  iiq   = iq - rismt%mp_site%isite_start + 1
  iv    = iuniq_to_isite(1, iq)
  isolV = isite_to_isolV(iv)
  iatom = isite_to_iatom(iv)
  sv    = solVs(isolV)%ljsig(iatom)
  ev    = solVs(isolV)%ljeps(iatom)
  !
  ! ... wall properties
  rho = wall_rho
  su  = wall_ljsig
  eu  = wall_ljeps
  !
  ! ... wall-solvent properties
  suv  = 0.5_DP * (sv + su)
  euv  = SQRT(ev * eu)
  rmax = rsmax * suv / alat
  rmin = RSMIN * suv / alat
  !
  ! ... offset for Laue-RISM
#if defined (__ESM_NOT_SYMMETRIC)
  r3offset = 0.0_DP
#else
  IF (MOD(n3, 2) == 0) THEN
    r3offset = 0.5_DP / DBLE(n3)
  ELSE
    r3offset = 0.0_DP
  END IF
#endif
  !
  ! ... calculate potential on each FFT grid
!$omp parallel do default(shared) private(ir, i1, i2, i3, offrange, r3, tau_z, &
!$omp             vw, zuv, sr, sr2, sr3, sr6, sr9)
  DO ir = 1, mr
    !
    ! ... create coordinate of a FFT grid
    IF (ir <= nr) THEN
      CALL fft_index_to_3d(ir, rismt%dfft, i1, i2, i3, offrange)
    ELSE
      offrange = .TRUE.
    END IF
    !
    IF (offrange) THEN
      rismt%uwr(ir, iiq) = 0.0_DP
      CYCLE
    END IF
    !
    r3 = r3offset + DBLE(i3) / DBLE(n3)
    IF (i3 >= (n3 - (n3 / 2))) THEN
      r3 =  r3 - 1.0_DP
    END IF
    !
    tau_z = r3 * at(3, 3)
    !
    ! ... potential from the wall
    zuv = rsign * (tau_z - wall_tau)
    IF (zuv < rmin) THEN
      zuv = rmin
    END IF
    !
    IF (zuv > rmax) THEN
      vw = 0.0_DP
      !
    ELSE
      sr   = suv / zuv / alat
      sr2  = sr  * sr
      sr3  = sr2 * sr
      sr6  = sr3 * sr3
      sr9  = sr6 * sr3
      IF (wall_lj6) THEN
        vw = tpi * rho * 4.0_DP * euv * suv * suv * suv * (sr9 / 90.0_DP - sr3 / 12.0_DP)
      ELSE
        vw = tpi * rho * 4.0_DP * euv * suv * suv * suv * sr9 / 90.0_DP
      END IF
      !
    END IF
    !
    rismt%uwr(ir, iiq) = vw
    !
  END DO
!$omp end parallel do
  !
END SUBROUTINE lj_setup_wall_x
!
!--------------------------------------------------------------------------
SUBROUTINE lj_get_wall_edge(tau0, v0)
  !--------------------------------------------------------------------------
  !
  ! ... calculate edge position of repulsive wall
  !
  USE kinds,   ONLY : DP
  USE solvmol, ONLY : get_nuniq_in_solVs
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: tau0 ! edge position of wall
  REAL(DP), INTENT(IN)  :: v0   ! threshold of LJ-potential
  !
  INTEGER :: nq
  INTEGER :: iq
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... calculate edge position of wall
  tau0 = 1.0E+99_DP
  DO iq = 1, nq
    CALl lj_get_wall_edge_x(iq, tau0, v0)
  END DO
  !
END SUBROUTINE lj_get_wall_edge
!
!--------------------------------------------------------------------------
SUBROUTINE lj_get_wall_edge_x(iq, tau0, v0)
  !--------------------------------------------------------------------------
  !
  ! ... calculate edge position of repulsive wall
  ! ... for a solvent's site
  !
  USE cell_base, ONLY : alat
  USE constants, ONLY : tpi
  USE kinds,     ONLY : DP
  USE solute,    ONLY : iwall, wall_rho, wall_ljeps, wall_ljsig
  USE solvmol,   ONLY : iuniq_to_isite, isite_to_isolV, isite_to_iatom, solVs
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: iq
  REAL(DP), INTENT(INOUT) :: tau0
  REAL(DP), INTENT(IN)    :: v0
  !
  INTEGER  :: iv
  INTEGER  :: isolV
  INTEGER  :: iatom
  REAL(DP) :: rho
  REAL(DP) :: ev, eu, euv
  REAL(DP) :: sv, su, suv
  REAL(DP) :: suv2, suv4, suv8, suv12
  REAL(DP) :: z0
  !
  IF(v0 <= 0.0_DP) THEN
    RETURN
  END IF
  !
  ! ... solvent properties
  iv    = iuniq_to_isite(1, iq)
  isolV = isite_to_isolV(iv)
  iatom = isite_to_iatom(iv)
  sv    = solVs(isolV)%ljsig(iatom)
  ev    = solVs(isolV)%ljeps(iatom)
  !
  ! ... wall properties
  rho = wall_rho
  su  = wall_ljsig
  eu  = wall_ljeps
  !
  ! ... wall-solvent properties
  suv = 0.5_DP * (sv + su)
  euv = SQRT(ev * eu)
  !
  ! ... calculate edge position of wall
  suv2  = suv  * suv
  suv4  = suv2 * suv2
  suv8  = suv4 * suv4
  suv12 = suv8 * suv4
  z0    = tpi * rho * 4.0_DP * euv * suv12 / 90.0_DP / v0
  !
  IF (z0 > 0.0_DP) THEN
    z0 = z0 ** (1.0_DP / 9.0_DP)
    tau0 = MIN(tau0, z0 / alat)
  END IF
  !
END SUBROUTINE lj_get_wall_edge_x
!
!--------------------------------------------------------------------------
SUBROUTINE lj_get_force(rismt, force, rsmax, ierr)
  !--------------------------------------------------------------------------
  !
  ! ... calculate solute-solvent's Lennard-Jones force
  ! ...
  ! ...           /                     [ 12*(x - X) [  sig  ]12    6*(x - X) [  sig  ]6 ]
  ! ...   - rho * | dr g(r) * 4 * esp * [ ---------- [-------]   -  --------- [-------]  ]
  ! ...           /                     [ |r - R|^2  [|r - R|]      |r - R|^2 [|r - R|]  ]
  ! ...
  !
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE ions_base, ONLY : nat
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  USE solvmol,   ONLY : get_nuniq_in_solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(OUT) :: force(3, nat)
  REAL(DP),        INTENT(IN)  :: rsmax
  INTEGER,         INTENT(OUT) :: ierr
  !
  INTEGER :: nq
  INTEGER :: iq
  LOGICAL :: laue
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%dfft%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... Laue-RISM or not
  laue = .FALSE.
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    laue = .TRUE.
  END IF
  !
  ! ... calculate Lennard-Jones force
  force = 0.0_DP
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    CALL lj_get_force_x(iq, rismt, force, rsmax, laue)
  END DO
  !
  CALL mp_sum(force, rismt%mp_site%inter_sitg_comm)
  CALL mp_sum(force, rismt%mp_site%intra_sitg_comm)
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE lj_get_force
!
!--------------------------------------------------------------------------
SUBROUTINE lj_get_force_x(iq, rismt, force, rsmax, laue)
  !--------------------------------------------------------------------------
  !
  ! ... calculate solute-solvent's Lennard-Jones force
  ! ... for a solvent's site
  !
  USE cell_base, ONLY : alat, at, omega
  USE constants, ONLY : eps32
  USE fft_types, ONLY : fft_index_to_3d
#if defined(_OPENMP)
  USE ions_base, ONLY : nat
#endif
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type
  USE solute,    ONLY : solU_nat, solU_tau, solU_ljeps, solU_ljsig, isup_to_iuni
  USE solvmol,   ONLY : solVs, iuniq_to_isite, iuniq_to_nsite, &
                      & isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  INTEGER,         INTENT(IN)    :: iq
  TYPE(rism_type), INTENT(IN)    :: rismt
  REAL(DP),        INTENT(INOUT) :: force(1:3, 1:*)
  REAL(DP),        INTENT(IN)    :: rsmax
  LOGICAL,         INTENT(IN)    :: laue
  !
  REAL(DP), PARAMETER :: RSMIN = RSMIN__
  !
  INTEGER  :: iiq
  INTEGER  :: iv
  INTEGER  :: nv
  INTEGER  :: isolV
  INTEGER  :: iatom
  INTEGER  :: ir, nr
  INTEGER  :: i1, i2, i3
  INTEGER  :: n1, n2, n3
  INTEGER  :: nx1, nx2, nx3
  INTEGER  :: iz
  INTEGER  :: ia, iia
  LOGICAL  :: offrange
  REAL(DP) :: fac
  REAL(DP) :: weight
  REAL(DP) :: rho_right
  REAL(DP) :: rho_left
  REAL(DP) :: rhog
  REAL(DP) :: tau_r(3)
  REAL(DP) :: r1, r2, r3
  REAL(DP) :: r3offset
  REAL(DP) :: ev, eu, euv
  REAL(DP) :: sv, su, suv
  REAL(DP) :: rmax
  REAL(DP) :: rmin
  REAL(DP) :: ruv2
  REAL(DP) :: xuv, yuv, zuv
  REAL(DP) :: sr2, sr6, sr12
#if defined(_OPENMP)
  REAL(DP), ALLOCATABLE :: fromp(:,:)
#endif
  !
  ! ... FFT box
  n1 = rismt%dfft%nr1
  n2 = rismt%dfft%nr2
  n3 = rismt%dfft%nr3
  nr = rismt%dfft%nr1x * rismt%dfft%my_nr2p * rismt%dfft%my_nr3p
  !
  weight = omega / DBLE(n1 * n2 * n3)
  !
  ! ... solvent properties
  iiq       = iq - rismt%mp_site%isite_start + 1
  iv        = iuniq_to_isite(1, iq)
  nv        = iuniq_to_nsite(iq)
  isolV     = isite_to_isolV(iv)
  iatom     = isite_to_iatom(iv)
  sv        = solVs(isolV)%ljsig(iatom)
  ev        = solVs(isolV)%ljeps(iatom)
  rho_right = DBLE(nv) * solVs(isolV)%density
  rho_left  = DBLE(nv) * solVs(isolV)%subdensity
  !
  ! ... offset for Laue-RISM
  IF (laue) THEN
#if defined (__ESM_NOT_SYMMETRIC)
    r3offset = 0.0_DP
#else
    IF (MOD(n3, 2) == 0) THEN
      r3offset = 0.5_DP / DBLE(n3)
    ELSE
      r3offset = 0.0_DP
    END IF
#endif
  END IF
  !
  ! ... calculate potential on each FFT grid
!$omp parallel default(shared) private(ir, i1, i2, i3, offrange, r1, r2, r3, tau_r, rhog, iz, &
!$omp          ia, iia, su, suv, rmax, rmin, xuv, yuv, zuv, ruv2, eu, euv, sr2, sr6, sr12, &
!$omp          fac, fromp)
#if defined(_OPENMP)
  ALLOCATE(fromp(3, nat))
  fromp = 0.0_DP
#endif
!$omp do
  DO ir = 1, nr
    !
    ! ... create coordinate of a FFT grid
    CALL fft_index_to_3d(ir, rismt%dfft, i1, i2, i3, offrange)
    IF (offrange) THEN
      CYCLE
    END IF
    !
    r1 = DBLE(i1) / DBLE(n1)
    r2 = DBLE(i2) / DBLE(n2)
    r3 = DBLE(i3) / DBLE(n3)
    !
    IF (laue) THEN
      r3 = r3 + r3offset
      IF (i3 >= (n3 - (n3 / 2))) THEN
        r3 =  r3 - 1.0_DP
      END IF
    END IF
    !
    tau_r(:) = r1 * at(:, 1) + r2 * at(:, 2) + r3 * at(:, 3)
    !
    IF (.NOT. laue) THEN
      rhog = rho_right * rismt%gr(ir, iiq)
      !
    ELSE
      IF (i3 < (n3 - (n3 / 2))) THEN
        iz = i3 + (n3 / 2)
      ELSE
        iz = i3 - n3 + (n3 / 2)
      END IF
      iz = iz + rismt%lfft%izcell_start
      !
      IF (iz <= rismt%lfft%izleft_gedge) THEN
        rhog = rho_left  * rismt%gr(ir, iiq)
      ELSE IF (iz >= rismt%lfft%izright_gedge) THEN
        rhog = rho_right * rismt%gr(ir, iiq)
      ELSE
        rhog = 0.0_DP
      END IF
    END IF
    !
    IF (ABS(rhog) < eps32) THEN
      CYCLE
    END IF
    !
    ! ... contribution from each solute's atom
    DO ia = 1, solU_nat
      iia  = isup_to_iuni(ia)
      su   = solU_ljsig(iia)
      suv  = 0.5_DP * (sv + su)
      rmax = rsmax * suv / alat
      rmin = RSMIN * suv / alat
      xuv  = tau_r(1) - solU_tau(1, ia)
      yuv  = tau_r(2) - solU_tau(2, ia)
      zuv  = tau_r(3) - solU_tau(3, ia)
      ruv2 = xuv * xuv + yuv * yuv + zuv * zuv
      IF (ruv2 > (rmax * rmax)) THEN
        CYCLE
      END IF
      IF (ruv2 < (rmin * rmin)) THEN
        ruv2 = rmin * rmin
      END IF
      !
      eu   = solU_ljeps(iia)
      euv  = SQRT(ev * eu)
      sr2  = suv * suv / ruv2 / alat / alat
      sr6  = sr2 * sr2 * sr2
      sr12 = sr6 * sr6
      fac  = 4.0_DP * euv * (12.0_DP * sr12 - 6.0_DP * sr6) / ruv2 / alat
#if defined(_OPENMP)
      fromp(1, iia) = fromp(1, iia) - weight * rhog * fac * xuv
      fromp(2, iia) = fromp(2, iia) - weight * rhog * fac * yuv
      fromp(3, iia) = fromp(3, iia) - weight * rhog * fac * zuv
#else
      force(1, iia) = force(1, iia) - weight * rhog * fac * xuv
      force(2, iia) = force(2, iia) - weight * rhog * fac * yuv
      force(3, iia) = force(3, iia) - weight * rhog * fac * zuv
#endif
    END DO
    !
  END DO
!$omp end do
#if defined(_OPENMP)
!$omp critical
  force(1:3, 1:nat) = force(1:3, 1:nat) + fromp(1:3, 1:nat)
!$omp end critical
  DEALLOCATE(fromp)
#endif
!$omp end parallel
  !
END SUBROUTINE lj_get_force_x
!
!--------------------------------------------------------------------------
SUBROUTINE lj_get_stress(rismt, sigma, rsmax, ierr)
  !--------------------------------------------------------------------------
  !
  ! ... calculate solute-solvent's Lennard-Jones stress
  !
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  USE solvmol,   ONLY : get_nuniq_in_solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(OUT) :: sigma(3, 3)
  REAL(DP),        INTENT(IN)  :: rsmax
  INTEGER,         INTENT(OUT) :: ierr
  !
  INTEGER :: nq
  INTEGER :: iq
  LOGICAL :: laue
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%dfft%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... Laue-RISM or not
  laue = .FALSE.
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    laue = .TRUE.
  END IF
  !
  ! ... calculate Lennard-Jones stress
  sigma = 0.0_DP
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    CALL lj_get_stress_x(iq, rismt, sigma, rsmax, laue)
  END DO
  !
  CALL mp_sum(sigma, rismt%mp_site%inter_sitg_comm)
  CALL mp_sum(sigma, rismt%mp_site%intra_sitg_comm)
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE lj_get_stress
!
!--------------------------------------------------------------------------
SUBROUTINE lj_get_stress_x(iq, rismt, sigma, rsmax, laue)
  !--------------------------------------------------------------------------
  !
  ! ... calculate solute-solvent's Lennard-Jones force
  ! ... for a solvent's site
  !
  USE cell_base, ONLY : alat, at, omega
  USE constants, ONLY : eps32
  USE fft_types, ONLY : fft_index_to_3d
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type
  USE solute,    ONLY : solU_nat, solU_tau, solU_ljeps, solU_ljsig, isup_to_iuni
  USE solvmol,   ONLY : solVs, iuniq_to_isite, iuniq_to_nsite, &
                      & isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  INTEGER,         INTENT(IN)    :: iq
  TYPE(rism_type), INTENT(IN)    :: rismt
  REAL(DP),        INTENT(INOUT) :: sigma(3, 3)
  REAL(DP),        INTENT(IN)    :: rsmax
  LOGICAL,         INTENT(IN)    :: laue
  !
  REAL(DP), PARAMETER :: RSMIN = RSMIN__
  !
  INTEGER  :: iiq
  INTEGER  :: iv
  INTEGER  :: nv
  INTEGER  :: isolV
  INTEGER  :: iatom
  INTEGER  :: ir, nr
  INTEGER  :: i1, i2, i3
  INTEGER  :: n1, n2, n3
  INTEGER  :: iz
  INTEGER  :: ia, iia
  LOGICAL  :: offrange
  REAL(DP) :: fac
  REAL(DP) :: weight
  REAL(DP) :: rho_right
  REAL(DP) :: rho_left
  REAL(DP) :: rhog
  REAL(DP) :: tau_r(3)
  REAL(DP) :: r1, r2, r3
  REAL(DP) :: r3offset
  REAL(DP) :: ev, eu, euv
  REAL(DP) :: sv, su, suv
  REAL(DP) :: rmax
  REAL(DP) :: rmin
  REAL(DP) :: ruv2
  REAL(DP) :: xuv, yuv, zuv
  REAL(DP) :: sr2, sr6, sr12
#if defined(_OPENMP)
  REAL(DP) :: sgomp(3, 3)
#endif
  !
  ! ... FFT box
  n1 = rismt%dfft%nr1
  n2 = rismt%dfft%nr2
  n3 = rismt%dfft%nr3
  nr = rismt%dfft%nr1x * rismt%dfft%my_nr2p * rismt%dfft%my_nr3p
  !
  weight = omega / DBLE(n1 * n2 * n3)
  !
  ! ... solvent properties
  iiq       = iq - rismt%mp_site%isite_start + 1
  iv        = iuniq_to_isite(1, iq)
  nv        = iuniq_to_nsite(iq)
  isolV     = isite_to_isolV(iv)
  iatom     = isite_to_iatom(iv)
  sv        = solVs(isolV)%ljsig(iatom)
  ev        = solVs(isolV)%ljeps(iatom)
  rho_right = DBLE(nv) * solVs(isolV)%density
  rho_left  = DBLE(nv) * solVs(isolV)%subdensity
  !
  ! ... offset for Laue-RISM
  IF (laue) THEN
#if defined (__ESM_NOT_SYMMETRIC)
    r3offset = 0.0_DP
#else
    IF (MOD(n3, 2) == 0) THEN
      r3offset = 0.5_DP / DBLE(n3)
    ELSE
      r3offset = 0.0_DP
    END IF
#endif
  END IF
  !
  ! ... calculate potential on each FFT grid
!$omp parallel default(shared) private(ir, i1, i2, i3, offrange, r1, r2, r3, tau_r, rhog, iz, &
!$omp          ia, iia, su, suv, rmax, rmin, xuv, yuv, zuv, ruv2, eu, euv, sr2, sr6, sr12, &
!$omp          fac, sgomp)
#if defined(_OPENMP)
  sgomp = 0.0_DP
#endif
!$omp do
  DO ir = 1, nr
    !
    ! ... create coordinate of a FFT grid
    CALL fft_index_to_3d(ir, rismt%dfft, i1, i2, i3, offrange)
    IF (offrange) THEN
      CYCLE
    END IF
    !
    r1 = DBLE(i1) / DBLE(n1)
    r2 = DBLE(i2) / DBLE(n2)
    r3 = DBLE(i3) / DBLE(n3)
    !
    IF (laue) THEN
      r3 = r3 + r3offset
      IF (i3 >= (n3 - (n3 / 2))) THEN
        r3 =  r3 - 1.0_DP
      END IF
    END IF
    !
    tau_r(:) = r1 * at(:, 1) + r2 * at(:, 2) + r3 * at(:, 3)
    !
    IF (.NOT. laue) THEN
      rhog = rho_right * rismt%gr(ir, iiq)
      !
    ELSE
      IF (i3 < (n3 - (n3 / 2))) THEN
        iz = i3 + (n3 / 2)
      ELSE
        iz = i3 - n3 + (n3 / 2)
      END IF
      iz = iz + rismt%lfft%izcell_start
      !
      IF (iz <= rismt%lfft%izleft_gedge) THEN
        rhog = rho_left  * rismt%gr(ir, iiq)
      ELSE IF (iz >= rismt%lfft%izright_gedge) THEN
        rhog = rho_right * rismt%gr(ir, iiq)
      ELSE
        rhog = 0.0_DP
      END IF
    END IF
    !
    IF (ABS(rhog) < eps32) THEN
      CYCLE
    END IF
    !
    ! ... contribution from each solute's atom
    DO ia = 1, solU_nat
      iia  = isup_to_iuni(ia)
      su   = solU_ljsig(iia)
      suv  = 0.5_DP * (sv + su)
      rmax = rsmax * suv / alat
      rmin = RSMIN * suv / alat
      xuv  = tau_r(1) - solU_tau(1, ia)
      yuv  = tau_r(2) - solU_tau(2, ia)
      zuv  = tau_r(3) - solU_tau(3, ia)
      ruv2 = xuv * xuv + yuv * yuv + zuv * zuv
      IF (ruv2 > (rmax * rmax)) THEN
        CYCLE
      END IF
      IF (ruv2 < (rmin * rmin)) THEN
        ruv2 = rmin * rmin
      END IF
      !
      eu   = solU_ljeps(iia)
      euv  = SQRT(ev * eu)
      sr2  = suv * suv / ruv2 / alat / alat
      sr6  = sr2 * sr2 * sr2
      sr12 = sr6 * sr6

#if defined(_OPENMP)
      ! TODO
      ! TODO set sgomp
      ! TODO
#else
      ! TODO
      ! TODO set sigma
      ! TODO
#endif
    END DO
    !
  END DO
!$omp end do
#if defined(_OPENMP)
!$omp critical
  sigma = sigma + sgomp
!$omp end critical
#endif
!$omp end parallel
  !
END SUBROUTINE lj_get_stress_x

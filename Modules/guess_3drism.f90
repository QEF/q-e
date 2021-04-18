!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE guess_3drism(rismt, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... create initial guess of 3D-RISM or Laue-RISM.
  !
  USE cell_base, ONLY : at, alat
  USE constants, ONLY : eps4, K_BOLTZMANN_RY
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE fft_types, ONLY : fft_index_to_3d
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_max
  USE rism,      ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  USE solvmol,   ONLY : get_nuniq_in_solVs, solVs, &
                      & iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER  :: nq
  INTEGER  :: iq
  INTEGER  :: iiq
  INTEGER  :: iv
  INTEGER  :: isolV
  INTEGER  :: iatom
  INTEGER  :: ir
  INTEGER  :: i1, i2, i3
  LOGICAL  :: offrange
  LOGICAL  :: laue
  REAL(DP) :: beta
  REAL(DP) :: qv
  REAL(DP) :: vlj
  REAL(DP) :: cs0
  REAL(DP) :: csmax
  REAL(DP) :: erf0
  !
  REAL(DP), PARAMETER :: VLJ_MAX  = eps4 ! Ry
  REAL(DP), PARAMETER :: CS_SCALE = 0.1_DP
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
  ! ... if no data, return as normally done
  IF (rismt%nsite < 1) THEN
    GOTO 1
  END IF
  !
  ! ... Laue-RISM or not
  laue = .FALSE.
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    laue = .TRUE.
  END IF
  !
  ! ... beta = 1 / (kB * T)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... create guess for each solvent's site
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    isolV = isite_to_isolV(iv)
    iatom = isite_to_iatom(iv)
    qv    = solVs(isolV)%charge(iatom)
    !
    ! ... set csr to be zero
    csmax = 0.0_DP
    rismt%csr(:, iiq) = 0.0_DP
    !
    ! ... set csr initially
    DO ir = 1, rismt%dfft%nr1x * rismt%dfft%my_nr2p * rismt%dfft%my_nr3p
      !
      CALL fft_index_to_3d(ir, rismt%dfft, i1, i2, i3, offrange)
      IF (offrange) THEN
        CYCLE
      END IF
      !
      vlj = rismt%uljr(ir, iiq)
      IF (laue) THEN  ! add LJ-wall
        vlj =  vlj + rismt%uwr(ir, iiq)
      END IF
      !
      IF (vlj >= VLJ_MAX) THEN
        cs0 = beta * qv * rismt%vlr(ir)
        csmax = MAX(csmax, ABS(cs0))
        rismt%csr(ir, iiq) = cs0
      END IF
      !
    END DO
    !
    CALL mp_max(csmax, rismt%mp_site%intra_sitg_comm)
    !
    ! ... correct csr to be smooth
    DO ir = 1, rismt%dfft%nr1x * rismt%dfft%my_nr2p * rismt%dfft%my_nr3p
      !
      CALL fft_index_to_3d(ir, rismt%dfft, i1, i2, i3, offrange)
      IF (offrange) THEN
        CYCLE
      END IF
      !
      IF (csmax > 0.0_DP) THEN
        cs0  = rismt%csr(ir, iiq)
        erf0 = erf(ABS(cs0) / (CS_SCALE * csmax))
        rismt%csr(ir, iiq) = cs0 * erf0 * erf0
      END IF
      !
    END DO
    !
  END DO
  !
  ! ... correction for Laue-RISM
  IF (laue) THEN
    CALL correct_edge()
  END IF
  !
  ! ... set dipole part for Laue-RISM
  IF (laue) THEN
    IF (rismt%nsite > 0) THEN
      rismt%cda = 0.0_DP
    END IF
  END IF
  !
  ! ... set Gxy=0 term for Laue-RISM
  IF (laue) THEN
    IF (rismt%nrzl * rismt%nsite > 0) THEN
      rismt%csg0 = 0.0_DP
    END IF
    !
    CALL corrgxy0_laue(rismt, .TRUE., rismt%csr, rismt%csg0, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      RETURN
    END IF
  END IF
  !
  ! ... normally done
1 CONTINUE
  ierr = IERR_RISM_NULL
  !
CONTAINS
  !
  SUBROUTINE correct_edge()
    IMPLICIT NONE
    !
    INTEGER  :: iz
    REAL(DP) :: z
    REAL(DP) :: z0
    !
    REAL(DP), PARAMETER :: Z_SCALE = 5.0_DP ! bohr
    !
    IF (rismt%nsite < 1) THEN
      RETURN
    END IF
    !
    z0 = 0.5_DP * at(3, 3)
    !
    DO ir = 1, rismt%dfft%nr1x * rismt%dfft%my_nr2p * rismt%dfft%my_nr3p
      !
      CALL fft_index_to_3d(ir, rismt%dfft, i1, i2, i3, offrange)
      IF (offrange) THEN
        CYCLE
      END IF
      !
      IF (i3 < (rismt%dfft%nr3 - (rismt%dfft%nr3 / 2))) THEN
        iz = i3 + (rismt%dfft%nr3 / 2)
      ELSE
        iz = i3 - rismt%dfft%nr3 + (rismt%dfft%nr3 / 2)
      END IF
      iz = iz + rismt%lfft%izcell_start
      !
      z = rismt%lfft%zleft + rismt%lfft%zoffset + rismt%lfft%zstep * DBLE(iz - 1)
      !
      IF (rismt%lfft%xright) THEN
        erf0 = erf(alat * (z0 - z) / Z_SCALE)
        rismt%csr(ir, :) = rismt%csr(ir, :) * (erf0 * erf0)
      END IF
      !
      IF (rismt%lfft%xleft) THEN
        erf0 = erf(alat * (z + z0) / Z_SCALE)
        rismt%csr(ir, :) = rismt%csr(ir, :) * (erf0 * erf0)
      END IF
      !
    END DO
    !
  END SUBROUTINE correct_edge
  !
END SUBROUTINE guess_3drism

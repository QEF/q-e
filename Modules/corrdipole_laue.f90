!
! Copyright (C) 2017 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE corrdipole_laue(rismt, lextract, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... sum (and extract) dipole parts of direct correlations, in Laue-RISM calculation.
  ! ...
  ! ... Variables:
  ! ...   lextract: if .TRUE.  extract dipole parts of direct correlations
  ! ...                        and sum Cs(r) + Cd(z),
  ! ...             if .FALSE. not extract dipole parts of direct correlations
  ! ...                        but sum Cs(r) + Cd(z)
  !
  USE cell_base, ONLY : alat
  USE constants, ONLY : K_BOLTZMANN_RY
  USE fft_types, ONLY : fft_index_to_3d
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,   ONLY : solVs, iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  LOGICAL,         INTENT(IN)    :: lextract
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER               :: iq
  INTEGER               :: iiq
  INTEGER               :: iv
  INTEGER               :: isolV
  INTEGER               :: iatom
  INTEGER               :: izsolv
  REAL(DP)              :: qv
  REAL(DP)              :: beta
  REAL(DP)              :: zstep
  REAL(DP)              :: zoffs
  REAL(DP)              :: zedge
  REAL(DP)              :: vline0
  REAL(DP)              :: vline1
  REAL(DP), ALLOCATABLE :: cd0(:)
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%dfft%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... is one-hand ?
  IF (rismt%lfft%xright .AND. rismt%lfft%xleft) THEN
    !
    IF (lextract) THEN
      IF (rismt%nsite > 0) THEN
        rismt%cda = 0.0_DP
      END IF
    END IF
    !
    IF (rismt%nr * rismt%nsite > 0) THEN
      rismt%csdr  = rismt%csr
      rismt%csdg0 = rismt%csg0
    END IF
    !
    ierr = IERR_RISM_NULL
    RETURN
  END IF
  !
  ! ... sum Cs(r) + Cd(z), then return
  IF (.NOT. lextract) THEN
    !
    IF (rismt%nr * rismt%nsite > 0) THEN
      rismt%csdr = 0.0_DP
    END IF
    !
    IF (rismt%nrzl * rismt%nsite > 0) THEN
      rismt%csdg0 = 0.0_DP
    END IF
    !
    CALL update_short_range(.FALSE.)
    !
    ierr = IERR_RISM_NULL
    RETURN
  END IF
  !
  ! ... beta = 1 / (kB * T)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... set some variables
  IF (rismt%lfft%xright) THEN
    izsolv = rismt%lfft%izright_start0
    !
    IF (rismt%lfft%gxystart > 1) THEN
      vline0 = AIMAG(rismt%vleft(1))
      vline1 = DBLE( rismt%vleft(1)) / alat
    END IF
    !
  ELSE !IF (rismt%lfft%xleft) THEN
    izsolv = rismt%lfft%izleft_end0
    !
    IF (rismt%lfft%gxystart > 1) THEN
      vline0 = AIMAG(rismt%vright(1))
      vline1 = DBLE( rismt%vright(1)) / alat
    END IF
  END IF
  !
  ! ... properties about length (in a.u.)
  zstep = alat * rismt%lfft%zstep
  zoffs = alat * (rismt%lfft%zleft + rismt%lfft%zoffset)
  zedge = zoffs + zstep * DBLE(izsolv - 1)
  !
  ! ... allocate memory
  IF (rismt%nsite > 0) THEN
    ALLOCATE(cd0(rismt%nsite))
    cd0 = 0.0_DP
  END IF
  !
  ! ... calculate amplitudes of dipole parts,
  ! ... which is included in short-range direct correlations
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    isolV = isite_to_isolV(iv)
    iatom = isite_to_iatom(iv)
    qv    = solVs(isolV)%charge(iatom)
    !
    IF (rismt%lfft%gxystart > 1) THEN
      cd0(iiq) = rismt%csg0(izsolv, iiq)              &  ! short-range
             & - beta * qv * DBLE(rismt%vlgz(izsolv)) &  ! long-range
             & + beta * qv * (vline1 * zedge + vline0)   ! from void-region
    END IF
  END DO
  !
  IF (rismt%nsite > 0) THEN
    CALL mp_sum(cd0, rismt%mp_site%intra_sitg_comm)
  END IF
  !
  ! ... update amplitudes of dipole parts of direct correlations
  IF (rismt%nsite > 0) THEN
    rismt%cda = cd0
  END IF
  !
  ! ... update short-range direct correlations
  IF (rismt%nr * rismt%nsite > 0) THEN
    rismt%csdr = 0.0_DP
  END IF
  !
  IF (rismt%nrzl * rismt%nsite > 0) THEN
    rismt%csdg0 = 0.0_DP
  END IF
  !
  CALL update_short_range(.TRUE.)
  !
  ! ... deallocate memoery
  IF (rismt%nsite > 0) THEN
    DEALLOCATE(cd0)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
CONTAINS
  !
  SUBROUTINE update_short_range(modify_csr)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: modify_csr
    !
    INTEGER :: ir
    INTEGER :: i1, i2, i3
    INTEGER :: iz
    INTEGER :: isite
    LOGICAL :: offrange
    !
    IF (rismt%nsite < 1) THEN
      RETURN
    END IF
    !
    ! ... update R-space
!$omp parallel do default(shared) private(ir, i1, i2, i3, iz, isite, offrange)
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
      IF (iz > rismt%lfft%izright_end .OR. iz < rismt%lfft%izleft_start) THEN
        CYCLE
      END IF
      IF (iz < rismt%lfft%izright_start .AND. iz > rismt%lfft%izleft_end) THEN
        CYCLE
      END IF
      !
      IF (modify_csr) THEN
        DO isite = 1, rismt%nsite
          rismt%csr(ir, isite) = rismt%csr(ir, isite) - rismt%cda(isite) * rismt%cdz(iz)
        END DO
      END IF
      !
      DO isite = 1, rismt%nsite
        rismt%csdr(ir, isite) = rismt%csr(ir, isite) + rismt%cda(isite) * rismt%cdz(iz)
      END DO
      !
    END DO
!$omp end parallel do
    !
    ! ... update Gxy=0 term
!$omp parallel do default(shared) private(iz, isite)
    DO iz = 1, rismt%nrzl
      !
      IF (modify_csr) THEN
        DO isite = 1, rismt%nsite
          rismt%csg0(iz, isite) = rismt%csg0(iz, isite) - rismt%cda(isite) * rismt%cdz(iz)
        END DO
      END IF
      !
      DO isite = 1, rismt%nsite
        rismt%csdg0(iz, isite) = rismt%csg0(iz, isite) + rismt%cda(isite) * rismt%cdz(iz)
      END DO
      !
    END DO
!$omp end parallel do
    !
  END SUBROUTINE update_short_range
  !
END SUBROUTINE corrdipole_laue

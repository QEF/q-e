!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE print_corr_vv(rismt, label, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... print 1D-RISM's correlation functions
  !
  USE constants, ONLY : K_BOLTZMANN_RY, BOHR_RADIUS_ANGS
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE io_files,  ONLY : tmp_dir, prefix
  USE io_global, ONLY : ionode
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum, mp_get
  USE rism,      ONLY : rism_type, ITYPE_1DRISM
  USE solvmol,   ONLY : nsolV, solVs, get_nuniq_in_solVs, get_nsite_in_solVs, &
                      & iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type),  INTENT(IN)  :: rismt
  CHARACTER(LEN=*), INTENT(IN)  :: label
  INTEGER,          INTENT(OUT) :: ierr
  !
  INTEGER               :: nv
  INTEGER               :: iun
  CHARACTER(LEN=256)    :: filrism
  REAL(DP)              :: beta
  REAL(DP), ALLOCATABLE :: zvv(:,:)
  !
  ! ... number of sites in solvents
  nv = get_nsite_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_1DRISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr /= rismt%ng) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nsite < (nv * (nv + 1) / 2)) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (.NOT. rismt%is_intra) THEN
    ierr = IERR_RISM_NULL
    RETURN
  END IF
  !
  ! ... beta = 1 / (kB * T)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ALLOCATE(zvv(rismt%nr, rismt%nsite))
  !
  ! ... write Gvv(r)
#if defined (__DEBUG_RISM)
  filrism = TRIM(tmp_dir) // TRIM(prefix) // '.1drism_gvv' // TRIM(ADJUSTL(label))
#else
  filrism = TRIM(tmp_dir) // TRIM(prefix) // '.1drism' // TRIM(ADJUSTL(label))
#endif
  CALL open_filrism(iun, TRIM(filrism))
  CALL write_comment(iun, 'Pair distribution function')
  CALL write_comment(iun, 'Gvv(r)')
  CALL write_comment(iun, '')
  CALL write_zvv(iun, rismt%gr)
  CALL close_filrism(iun)
  !
#if defined (__DEBUG_RISM)
  ! ... write Hvv(r)
  filrism = TRIM(tmp_dir) // TRIM(prefix) // '.1drism_hvv' // TRIM(ADJUSTL(label))
  CALL open_filrism(iun, TRIM(filrism))
  CALL write_comment(iun, 'Total correlation function')
  CALL write_comment(iun, 'Hvv(r)')
  CALL write_comment(iun, '')
  CALL write_zvv(iun, rismt%hr)
  CALL close_filrism(iun)
  !
  ! ... write Cvv(r)
  filrism = TRIM(tmp_dir) // TRIM(prefix) // '.1drism_cvv' // TRIM(ADJUSTL(label))
  CALL open_filrism(iun, TRIM(filrism))
  CALL write_comment(iun, 'Direct correlation function')
  CALL write_comment(iun, 'Cvv(r)')
  CALL write_comment(iun, '')
  zvv = rismt%csr - beta * rismt%ulr
  CALL write_zvv(iun, zvv)
  CALL close_filrism(iun)
  !
  ! ... write Uvv(r) / kB * T
  filrism = TRIM(tmp_dir) // TRIM(prefix) // '.1drism_uvv' // TRIM(ADJUSTL(label))
  CALL open_filrism(iun, TRIM(filrism))
  CALL write_comment(iun, 'Inter-site potential')
  CALL write_comment(iun, 'Uvv(r) / kB * T')
  CALL write_comment(iun, '')
  zvv = beta * (rismt%usr + rismt%ulr)
  CALL write_zvv(iun, zvv)
  CALL close_filrism(iun)
  !
#endif
  DEALLOCATE(zvv)
  !
  ierr = IERR_RISM_NULL
  !
CONTAINS
  !
  SUBROUTINE open_filrism(iun, filrism)
    IMPLICIT NONE
    INTEGER,          INTENT(OUT) :: iun
    CHARACTER(LEN=*), INTENT(IN)  :: filrism
    !
    INTEGER :: ios
    !
    INTEGER, EXTERNAL :: find_free_unit
    !
    iun = find_free_unit()
    IF (ionode) THEN
      OPEN(unit=iun, file=filrism, status='unknown', form='formatted', action='write', iostat=ios)
      ios = ABS(ios)
    ELSE
      ios = 0
    END IF
    !
    CALL mp_sum(ios, rismt%mp_task%itask_comm)
    IF (ios > 0) THEN
      CALL errore('print_corr_vv', 'cannot open file' // TRIM(filrism), ios)
    END IF
  END SUBROUTINE open_filrism
  !
  SUBROUTINE close_filrism(iun)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iun
    !
    IF (ionode) THEN
      CLOSE(unit=iun)
    END IF
  END SUBROUTINE close_filrism
  !
  SUBROUTINE write_comment(iun, str)
    IMPLICIT NONE
    INTEGER,          INTENT(IN) :: iun
    CHARACTER(LEN=*), INTENT(IN) :: str
    !
    IF (ionode) THEN
      WRITE(iun, '("# ", A)') str
    END IF
  END SUBROUTINE write_comment
  !
  SUBROUTINE write_zvv(iun, zvv)
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: iun
    REAL(DP), INTENT(IN) :: zvv(:,:)
    !
    INTEGER,  PARAMETER :: LEN_STR = 6
    INTEGER,  PARAMETER :: LEN_COL = 15
    REAL(DP), PARAMETER :: RMAX = 100.0_DP ! angs.
    !
    INTEGER                             :: io_proc
    INTEGER                             :: iproc
    INTEGER                             :: ir
    INTEGER                             :: mr
    INTEGER                             :: nr
    INTEGER                             :: nq
    INTEGER                             :: nqq
    INTEGER                             :: iqq
    INTEGER                             :: iq1, iq2
    INTEGER                             :: iv1, iv2
    INTEGER                             :: iw1, iw2
    INTEGER                             :: isolV1
    INTEGER                             :: isolV2
    INTEGER                             :: iatom1
    INTEGER                             :: iatom2
    INTEGER                             :: ncount
    REAL(DP)                            :: r
    REAL(DP),               ALLOCATABLE :: yvv(:,:)
    INTEGER,                ALLOCATABLE :: ivv(:)
    INTEGER                             :: ivv_
    INTEGER,                ALLOCATABLE :: qqhash(:)
    INTEGER                             :: qqhash_
    CHARACTER(LEN=LEN_STR)              :: ssolV1
    CHARACTER(LEN=LEN_STR)              :: ssolV2
    CHARACTER(LEN=LEN_STR)              :: satom1
    CHARACTER(LEN=LEN_STR)              :: satom2
    CHARACTER(LEN=LEN_COL), ALLOCATABLE :: title1(:)
    CHARACTER(LEN=LEN_COL)              :: title1_
    CHARACTER(LEN=LEN_COL), ALLOCATABLE :: title2(:)
    CHARACTER(LEN=LEN_COL)              :: title2_
    !
    ! ... find the index of the ionode
    io_proc = 0
    IF (ionode) THEN
      io_proc = rismt%mp_task%me_task
    END IF
    CALL mp_sum(io_proc, rismt%mp_task%itask_comm)
    !
    ! ... set index
    nq = get_nuniq_in_solVs()
    nqq = nq * (nq + 1) / 2
    ALLOCATE(ivv(nqq))
    ALLOCATE(qqhash(nqq))
    ALLOCATE(title1(nqq))
    ALLOCATE(title2(nqq))
    !
    iqq = 0
    DO iq1 = 1, nq
      iv1 = iuniq_to_isite(1, iq1)
      DO iq2 = iq1, nq
        iv2 = iuniq_to_isite(1, iq2)
        iw1 = MIN(iv1, iv2)
        iw2 = MAX(iv1, iv2)
        isolV1 = isite_to_isolV(iw1)
        isolV2 = isite_to_isolV(iw2)
        iatom1 = isite_to_iatom(iw1)
        iatom2 = isite_to_iatom(iw2)
        ssolV1 = ADJUSTL(solVs(isolV1)%name)
        ssolV2 = ADJUSTL(solVs(isolV2)%name)
        satom1 = ADJUSTL(solVs(isolV1)%aname(iatom1))
        satom2 = ADJUSTL(solVs(isolV2)%aname(iatom2))
        iqq = iqq + 1
        ivv(iqq) = iw2 * (iw2 - 1) / 2 + iw1
        qqhash(iqq) = nsolV * (isolV1 - 1) + isolV2
        title1(iqq) = '  ' // TRIM(ssolV1) // ':' // TRIM(ssolV2)
        title2(iqq) = '  ' // TRIM(satom1) // ':' // TRIM(satom2)
      END DO
    END DO
    !
    ! ... sort index
    ncount = 1
    DO WHILE (ncount > 0)
      ncount = 0
      DO iqq = 1, (nqq - 1)
        IF (qqhash(iqq) > qqhash(iqq + 1)) THEN
          ncount = ncount + 1
          ivv_    = ivv(   iqq)
          qqhash_ = qqhash(iqq)
          title1_ = title1(iqq)
          title2_ = title2(iqq)
          ivv(   iqq) = ivv(   iqq + 1)
          qqhash(iqq) = qqhash(iqq + 1)
          title1(iqq) = title1(iqq + 1)
          title2(iqq) = title2(iqq + 1)
          ivv(   iqq + 1) = ivv_
          qqhash(iqq + 1) = qqhash_
          title1(iqq + 1) = title1_
          title2(iqq + 1) = title2_
        END IF
      END DO
    END DO
    !
    title1_ = ''
    DO iqq = 1, nqq
      IF (title1(iqq) == title1_) THEN
        title1(iqq) = ''
      ELSE
        title1_ = title1(iqq)
      END IF
    END DO
    !
    ! ... write header
    IF (rismt%mp_task%me_task == io_proc) THEN
      WRITE(iun, '(15X)', advance='no')
      DO iqq = 1, nqq
        WRITE(iun, '(A15)', advance='no') title1(iqq)
      END DO
      WRITE(iun, '()', advance='yes')
      !
      WRITE(iun, '("  r/Angstrom   ")', advance='no')
      DO iqq = 1, nqq
        WRITE(iun, '(A15)', advance='no') title2(iqq)
      END DO
      WRITE(iun, '()', advance='yes')
    END IF
    !
    ! ... write zvv
    DO iproc = 0, (rismt%mp_task%nproc_task - 1)
      mr = rismt%mp_task%idis_vecs(iproc + 1)
      nr = rismt%mp_task%ilen_vecs(iproc + 1)
      !
      IF (rismt%mp_task%me_task == io_proc) THEN
        ALLOCATE(yvv(nr, rismt%nsite))
      ELSE
        ALLOCATE(yvv(1, 1))
      END IF
      !
      CALL mp_get(yvv, zvv, rismt%mp_task%me_task, &
      & io_proc, iproc, iproc + 1, rismt%mp_task%itask_comm)
      !
      IF (rismt%mp_task%me_task == io_proc) THEN
        DO ir = 1, nr
          r = rismt%rfft%rgrid(ir + mr) * BOHR_RADIUS_ANGS
#if !defined (__DEBUG_RISM)
          IF (r > RMAX) EXIT
#endif
          !
          WRITE(iun, '(E15.6e3)', advance='no') r
          DO iqq = 1, nqq
            WRITE(iun, '(E15.6e3)', advance='no') yvv(ir, ivv(iqq))
          END DO
          WRITE(iun, '()', advance='yes')
        END DO
      END IF
      !
      DEALLOCATE(yvv)
    END DO
    !
    DEALLOCATE(ivv)
    DEALLOCATE(qqhash)
    DEALLOCATE(title1)
    DEALLOCATE(title2)
  END SUBROUTINE write_zvv
  !
END SUBROUTINE print_corr_vv

!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE io_rism_xml
  !----------------------------------------------------------------------------
  !
  USE err_rism,    ONLY : stop_by_err_rism, IERR_RISM_INCORRECT_DATA_TYPE
  USE io_files,    ONLY : restart_dir, create_directory
  USE io_global,   ONLY : ionode, ionode_id
  USE kinds,       ONLY : DP
  USE rism,        ONLY : rism_type, ITYPE_1DRISM, ITYPE_3DRISM, ITYPE_LAUERISM
  USE xml_io_rism, ONLY : write_1drism_xml, read_1drism_xml, &
                        & write_3drism_xml, read_3drism_xml, &
                        & write_lauerism_xml, read_lauerism_xml, &
                        & write_lauedipole_xml, read_lauedipole_xml, &
                        & write_lauegxy0_xml, read_lauegxy0_xml
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... public components
  PUBLIC :: write_1drism
  PUBLIC :: read_1drism
  PUBLIC :: write_3drism
  PUBLIC :: read_3drism
  !
CONTAINS
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_1drism(rismt, extension)
    !------------------------------------------------------------------------
    !
    ! ... this routine writes the 1D-RISM's data in xml format into the
    ! ... '.save' directory
    ! ... the '.save' directory is created if not already present
    !
    IMPLICIT NONE
    !
    TYPE(rism_type),            INTENT(IN) :: rismt
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: extension
    !
    CHARACTER(LEN=256) :: dirname
    CHARACTER(LEN=256) :: file_base
    CHARACTER(LEN=256) :: ext
    !
    ! ... check data type
    IF (rismt%itype /= ITYPE_1DRISM) THEN
      CALL stop_by_err_rism('write_1drism', IERR_RISM_INCORRECT_DATA_TYPE)
    END IF
    !
    IF (rismt%nr /= rismt%ng) THEN
      CALL stop_by_err_rism('write_1drism', IERR_RISM_INCORRECT_DATA_TYPE)
    END IF
    !
    ! ... set file name
    dirname = restart_dir()
    CALL create_directory(dirname)
    !
    ext = ' '
    IF (PRESENT(extension)) THEN
      ext = '.' // TRIM(extension)
    END IF
    !
    ! ... return if out of intra group
    IF (.NOT. rismt%is_intra) THEN
      RETURN
    END IF
    !
    ! ... write csr
    file_base = TRIM(dirname) // '/1d-rism_csvv_r' // TRIM(ext)
    !
    CALL write_1drism_xml(file_base, rismt%csr, 'Csvv(r)', rismt%mp_task%nvec, &
       & rismt%nsite, rismt%mp_task%idis_vecs, rismt%mp_task%ilen_vecs, &
       & ionode, rismt%mp_task%itask_comm)
    !
    ! ... write hr
    file_base = TRIM(dirname) // '/1d-rism_hvv_r' // TRIM(ext)
    !
    CALL write_1drism_xml(file_base, rismt%hr, 'Hvv(r)', rismt%mp_task%nvec, &
       & rismt%nsite, rismt%mp_task%idis_vecs, rismt%mp_task%ilen_vecs, &
       & ionode, rismt%mp_task%itask_comm)
    !
    ! ... write gr
    file_base = TRIM(dirname) // '/1d-rism_gvv_r' // TRIM(ext)
    !
    CALL write_1drism_xml(file_base, rismt%gr, 'Gvv(r)', rismt%mp_task%nvec, &
       & rismt%nsite, rismt%mp_task%idis_vecs, rismt%mp_task%ilen_vecs, &
       & ionode, rismt%mp_task%itask_comm)
    !
    ! ... write csg
    file_base = TRIM(dirname) // '/1d-rism_csvv_g' // TRIM(ext)
    !
    CALL write_1drism_xml(file_base, rismt%csg, 'Csvv(g)', rismt%mp_task%nvec, &
       & rismt%nsite, rismt%mp_task%idis_vecs, rismt%mp_task%ilen_vecs, &
       & ionode, rismt%mp_task%itask_comm)
    !
    ! ... write hg
    file_base = TRIM(dirname) // '/1d-rism_hvv_g' // TRIM(ext)
    !
    CALL write_1drism_xml(file_base, rismt%hg, 'Hvv(g)', rismt%mp_task%nvec, &
       & rismt%nsite, rismt%mp_task%idis_vecs, rismt%mp_task%ilen_vecs, &
       & ionode, rismt%mp_task%itask_comm)
    !
  END SUBROUTINE write_1drism
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_1drism(rismt, extension)
    !------------------------------------------------------------------------
    !
    ! ... this routine reads the 1D-RISM's data in xml format from the
    ! ... files saved into the '.save' directory
    !
    IMPLICIT NONE
    !
    TYPE(rism_type),            INTENT(INOUT) :: rismt
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: extension
    !
    CHARACTER(LEN=256) :: dirname
    CHARACTER(LEN=256) :: file_base
    CHARACTER(LEN=256) :: ext
    !
    ! ... check data type
    IF (rismt%itype /= ITYPE_1DRISM) THEN
      CALL stop_by_err_rism('read_1drism', IERR_RISM_INCORRECT_DATA_TYPE)
    END IF
    !
    IF (rismt%nr /= rismt%ng) THEN
      CALL stop_by_err_rism('read_1drism', IERR_RISM_INCORRECT_DATA_TYPE)
    END IF
    !
    ! ... set file name
    dirname = restart_dir()
    ext = ' '
    IF (PRESENT(extension)) THEN
      ext = '.' // TRIM(extension)
    END IF
    !
    ! ... return if out of intra group
    IF (.NOT. rismt%is_intra) THEN
      RETURN
    END IF
    !
    ! ... read csr
    file_base = TRIM(dirname) // '/1d-rism_csvv_r' // TRIM(ext)
    !
    CALL read_1drism_xml(file_base, rismt%csr, rismt%mp_task%nvec, &
       & rismt%nsite, rismt%mp_task%idis_vecs, rismt%mp_task%ilen_vecs, &
       & ionode, ionode_id, rismt%mp_task%itask_comm)
    !
    ! ... read hr
    file_base = TRIM(dirname) // '/1d-rism_hvv_r' // TRIM(ext)
    !
    CALL read_1drism_xml(file_base, rismt%hr, rismt%mp_task%nvec, &
       & rismt%nsite, rismt%mp_task%idis_vecs, rismt%mp_task%ilen_vecs, &
       & ionode, ionode_id, rismt%mp_task%itask_comm)
    !
    ! ... read gr
    file_base = TRIM(dirname) // '/1d-rism_gvv_r' // TRIM(ext)
    !
    CALL read_1drism_xml(file_base, rismt%gr, rismt%mp_task%nvec, &
       & rismt%nsite, rismt%mp_task%idis_vecs, rismt%mp_task%ilen_vecs, &
       & ionode, ionode_id, rismt%mp_task%itask_comm)
    !
    ! ... read csg
    file_base = TRIM(dirname) // '/1d-rism_csvv_g' // TRIM(ext)
    !
    CALL read_1drism_xml(file_base, rismt%csg, rismt%mp_task%nvec, &
       & rismt%nsite, rismt%mp_task%idis_vecs, rismt%mp_task%ilen_vecs, &
       & ionode, ionode_id, rismt%mp_task%itask_comm)
    !
    ! ... read hg
    file_base = TRIM(dirname) // '/1d-rism_hvv_g' // TRIM(ext)
    !
    CALL read_1drism_xml(file_base, rismt%hg, rismt%mp_task%nvec, &
       & rismt%nsite, rismt%mp_task%idis_vecs, rismt%mp_task%ilen_vecs, &
       & ionode, ionode_id, rismt%mp_task%itask_comm)
    !
  END SUBROUTINE read_1drism
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_3drism(rismt, ecut, gamma_only, extension)
    !------------------------------------------------------------------------
    !
    ! ... this routine writes the 3D-RISM's data in xml format into the
    ! ... '.save' directory
    ! ... the '.save' directory is created if not already present
    !
    IMPLICIT NONE
    !
    TYPE(rism_type),            INTENT(IN) :: rismt
    REAL(DP),                   INTENT(IN) :: ecut
    LOGICAL,                    INTENT(IN) :: gamma_only
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: extension
    !
    CHARACTER(LEN=256) :: dirname
    CHARACTER(LEN=256) :: file_base
    CHARACTER(LEN=256) :: datname
    CHARACTER(LEN=256) :: ext
    !
    REAL(DP)    :: rdummy(1, 1)
    REAL(DP)    :: ddummy(1)
    COMPLEX(DP) :: cdummy(1, 1)
    !
    ! ... check data type
    IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
      CALL stop_by_err_rism('write_3drism', IERR_RISM_INCORRECT_DATA_TYPE)
    END IF
    !
    IF (rismt%nr < rismt%dfft%nnr) THEN
      CALL stop_by_err_rism('write_3drism', IERR_RISM_INCORRECT_DATA_TYPE)
    END IF
    !
    IF (rismt%itype == ITYPE_LAUERISM) THEN
      IF (rismt%nrzl < rismt%lfft%nrz) THEN
        CALL stop_by_err_rism('write_3drism', IERR_RISM_INCORRECT_DATA_TYPE)
      END IF
      !
      IF (rismt%ngxy < rismt%lfft%ngxy) THEN
        CALL stop_by_err_rism('write_3drism', IERR_RISM_INCORRECT_DATA_TYPE)
      END IF
      !
    END IF
    !
    ! ... set file name
    dirname = restart_dir()
    CALL create_directory(dirname)
    !
    ext = ' '
    IF (PRESENT(extension)) THEN
      ext = '.' // TRIM(extension)
    END IF
    !
    ! ... write csr
    file_base = TRIM(dirname) // '/3d-rism_csuv_r' // TRIM(ext)
    datname   = 'Csuv(r)'
    !
    IF (rismt%dfft%nnr * rismt%nsite > 0) THEN
      CALL write_3drism_x(rismt, rismt%dfft%nnr, rismt%nsite, &
                        & rismt%csr, ecut, file_base, datname)
    ELSE
      CALL write_3drism_x(rismt, rismt%dfft%nnr, rismt%nsite, &
                        & rdummy, ecut, file_base, datname)
    END IF
    !
    ! ... write hr
    file_base = TRIM(dirname) // '/3d-rism_huv_r' // TRIM(ext)
    datname   = 'Huv(r)'
    !
    IF (rismt%dfft%nnr * rismt%nsite > 0) THEN
      CALL write_3drism_x(rismt, rismt%dfft%nnr, rismt%nsite, &
                        & rismt%hr, ecut, file_base, datname)
    ELSE
      CALL write_3drism_x(rismt, rismt%dfft%nnr, rismt%nsite, &
                        & rdummy, ecut, file_base, datname)
    END IF
    !
    ! ... write gr
    file_base = TRIM(dirname) // '/3d-rism_guv_r' // TRIM(ext)
    datname   = 'Guv(r)'
    !
    IF (rismt%dfft%nnr * rismt%nsite > 0) THEN
      CALL write_3drism_x(rismt, rismt%dfft%nnr, rismt%nsite, &
                        & rismt%gr, ecut, file_base, datname)
    ELSE
      CALL write_3drism_x(rismt, rismt%dfft%nnr, rismt%nsite, &
                        & rdummy, ecut, file_base, datname)
    END IF
    !
    IF (rismt%itype == ITYPE_LAUERISM) THEN
      !
      ! ... write cd(dipole)
      file_base = TRIM(dirname) // '/3d-rism_cduv' // TRIM(ext)
      datname   = 'Cduv'
      !
      IF (rismt%nsite > 0) THEN
        CALL write_lauedipole_x(rismt, rismt%nsite, rismt%cda, file_base, datname)
      ELSE
        CALL write_lauedipole_x(rismt, rismt%nsite, ddummy, file_base, datname)
      END IF
      !
      ! ... write cs(Gxy=0)
      file_base = TRIM(dirname) // '/3d-rism_csuv_0' // TRIM(ext)
      datname   = 'Csuv(0,z)'
      !
      IF (rismt%nsite > 0) THEN
        CALL write_lauegxy0_x(rismt, rismt%nrzl, rismt%nsite, &
                            & rismt%csg0, file_base, datname)
      ELSE
        CALL write_lauegxy0_x(rismt, rismt%nrzl, rismt%nsite, &
                            & rdummy, file_base, datname)
      END IF
      !
      ! ... write hs(laue)
      file_base = TRIM(dirname) // '/3d-rism_hsuv_l' // TRIM(ext)
      datname   = 'Hsuv(gxy,z)'
      !
      IF (rismt%nrzl * rismt%lfft%ngxy * rismt%nsite > 0) THEN
        CALL write_lauerism_x(rismt, (rismt%nrzl * rismt%lfft%ngxy), rismt%nsite, &
                            & rismt%hsgz, ecut, gamma_only, file_base, datname)
      ELSE
        CALL write_lauerism_x(rismt, (rismt%nrzl * rismt%lfft%ngxy), rismt%nsite, &
                            & cdummy, ecut, gamma_only, file_base, datname)
      END IF
      !
      ! ... write hl(laue)
      file_base = TRIM(dirname) // '/3d-rism_hluv_l' // TRIM(ext)
      datname   = 'Hluv(gxy,z)'
      !
      IF (rismt%nrzl * rismt%lfft%ngxy * rismt%nsite > 0) THEN
        CALL write_lauerism_x(rismt, (rismt%nrzl * rismt%lfft%ngxy), rismt%nsite, &
                            & rismt%hlgz, ecut, gamma_only, file_base, datname)
      ELSE
        CALL write_lauerism_x(rismt, (rismt%nrzl * rismt%lfft%ngxy), rismt%nsite, &
                            & cdummy, ecut, gamma_only, file_base, datname)
      END IF
      !
    END IF
    !
  END SUBROUTINE write_3drism
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_3drism_x(rismt, nr, nsite, zr, ecut, file_base, name)
    !------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE(rism_type),  INTENT(IN) :: rismt
    INTEGER,          INTENT(IN) :: nr
    INTEGER,          INTENT(IN) :: nsite
    REAL(DP),         INTENT(IN) :: zr(nr, nsite)
    REAL(DP),         INTENT(IN) :: ecut
    CHARACTER(LEN=*), INTENT(IN) :: file_base
    CHARACTER(LEN=*), INTENT(IN) :: name
    !
    CALL write_3drism_xml(file_base, zr(1:nr, 1:nsite), name, &
                        & rismt%mp_site%nsite, rismt%mp_site%isite_start, rismt%mp_site%isite_end, &
                        & ecut, rismt%dfft, ionode, rismt%mp_site%inter_sitg_comm)
    !
  END SUBROUTINE write_3drism_x
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_lauerism_x(rismt, nl, nsite, zl, ecut, gamma_only, file_base, name)
    !------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE(rism_type),  INTENT(IN) :: rismt
    INTEGER,          INTENT(IN) :: nl
    INTEGER,          INTENT(IN) :: nsite
    COMPLEX(DP),      INTENT(IN) :: zl(nl, nsite)
    REAL(DP),         INTENT(IN) :: ecut
    LOGICAL,          INTENT(IN) :: gamma_only
    CHARACTER(LEN=*), INTENT(IN) :: file_base
    CHARACTER(LEN=*), INTENT(IN) :: name
    !
    CALL write_lauerism_xml(file_base, zl(1:nl, 1:nsite), name, &
                          & rismt%mp_site%nsite, rismt%mp_site%isite_start, rismt%mp_site%isite_end, &
                          & ecut, gamma_only, rismt%lfft, ionode, &
                          & rismt%mp_site%intra_sitg_comm, rismt%mp_site%inter_sitg_comm)
    !
  END SUBROUTINE write_lauerism_x
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_lauedipole_x(rismt, nsite, zd, file_base, name)
    !------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE(rism_type),  INTENT(IN) :: rismt
    INTEGER,          INTENT(IN) :: nsite
    REAL(DP),         INTENT(IN) :: zd(nsite)
    CHARACTER(LEN=*), INTENT(IN) :: file_base
    CHARACTER(LEN=*), INTENT(IN) :: name
    !
    CALL write_lauedipole_xml(file_base, zd(1:nsite), name, &
                            & rismt%mp_site%nsite, rismt%mp_site%isite_start, rismt%mp_site%isite_end, &
                            & ionode, rismt%mp_site%intra_sitg_comm, rismt%mp_site%inter_sitg_comm)
    !
  END SUBROUTINE write_lauedipole_x
  !
  !------------------------------------------------------------------------
  SUBROUTINE write_lauegxy0_x(rismt, nl, nsite, zl, file_base, name)
    !------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE(rism_type),  INTENT(IN) :: rismt
    INTEGER,          INTENT(IN) :: nl
    INTEGER,          INTENT(IN) :: nsite
    REAL(DP),         INTENT(IN) :: zl(nl, nsite)
    CHARACTER(LEN=*), INTENT(IN) :: file_base
    CHARACTER(LEN=*), INTENT(IN) :: name
    !
    CALL write_lauegxy0_xml(file_base, zl(1:nl, 1:nsite), name, &
                          & rismt%mp_site%nsite, rismt%mp_site%isite_start, rismt%mp_site%isite_end, &
                          & rismt%lfft, ionode, &
                          & rismt%mp_site%intra_sitg_comm, rismt%mp_site%inter_sitg_comm)
    !
  END SUBROUTINE write_lauegxy0_x
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_3drism(rismt, ecut, extension)
    !------------------------------------------------------------------------
    !
    ! ... this routine reads the 3D-RISM's data in xml format from the
    ! ... files saved into the '.save' directory
    !
    IMPLICIT NONE
    !
    TYPE(rism_type),            INTENT(INOUT) :: rismt
    REAL(DP),                   INTENT(IN)    :: ecut
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: extension
    !
    CHARACTER(LEN=256) :: dirname
    CHARACTER(LEN=256) :: file_base
    CHARACTER(LEN=256) :: ext
    !
    REAL(DP)    :: rdummy(1, 1)
    REAL(DP)    :: ddummy(1)
    COMPLEX(DP) :: cdummy(1, 1)
    !
    ! ... check data type
    IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
      CALL stop_by_err_rism('read_3drism', IERR_RISM_INCORRECT_DATA_TYPE)
    END IF
    !
    IF (rismt%nr < rismt%dfft%nnr) THEN
      CALL stop_by_err_rism('read_3drism', IERR_RISM_INCORRECT_DATA_TYPE)
    END IF
    !
    ! ... set file name
    dirname = restart_dir()
    ext = ' '
    IF (PRESENT(extension)) THEN
      ext = '.' // TRIM(extension)
    END IF
    !
    ! ... read csr
    file_base = TRIM(dirname) // '/3d-rism_csuv_r' // TRIM(ext)
    !
    IF (rismt%dfft%nnr * rismt%nsite > 0) THEN
      CALL read_3drism_x(rismt, rismt%dfft%nnr, rismt%nsite, &
                       & rismt%csr, ecut, file_base)
    ELSE
      CALL read_3drism_x(rismt, rismt%dfft%nnr, rismt%nsite, &
                       & rdummy, ecut, file_base)
    END IF
    !
    ! ... read hr
    file_base = TRIM(dirname) // '/3d-rism_huv_r' // TRIM(ext)
    !
    IF (rismt%dfft%nnr * rismt%nsite > 0) THEN
      CALL read_3drism_x(rismt, rismt%dfft%nnr, rismt%nsite, &
                       & rismt%hr, ecut, file_base)
    ELSE
      CALL read_3drism_x(rismt, rismt%dfft%nnr, rismt%nsite, &
                       & rdummy, ecut, file_base)
    END IF
    !
    ! ... read gr
    file_base = TRIM(dirname) // '/3d-rism_guv_r' // TRIM(ext)
    !
    IF (rismt%dfft%nnr * rismt%nsite > 0) THEN
      CALL read_3drism_x(rismt, rismt%dfft%nnr, rismt%nsite, &
                       & rismt%gr, ecut, file_base)
    ELSE
      CALL read_3drism_x(rismt, rismt%dfft%nnr, rismt%nsite, &
                       & rdummy, ecut, file_base)
    END IF
    !
    IF (rismt%itype == ITYPE_LAUERISM) THEN
      !
      ! ... read cd(dipole)
      file_base = TRIM(dirname) // '/3d-rism_cduv' // TRIM(ext)
      !
      IF (rismt%nsite > 0) THEN
        CALL read_lauedipole_x(rismt, rismt%nsite, rismt%cda, file_base)
      ELSE
        CALL read_lauedipole_x(rismt, rismt%nsite, ddummy, file_base)
      END IF
      !
      ! ... read cs(Gxy=0)
      file_base = TRIM(dirname) // '/3d-rism_csuv_0' // TRIM(ext)
      !
      IF (rismt%nsite > 0) THEN
        CALL read_lauegxy0_x(rismt, rismt%nrzl, rismt%nsite, &
                           & rismt%csg0, file_base)
      ELSE
        CALL read_lauegxy0_x(rismt, rismt%nrzl, rismt%nsite, &
                           & rdummy, file_base)
      END IF
      !
      ! ... read hs(laue)
      file_base = TRIM(dirname) // '/3d-rism_hsuv_l' // TRIM(ext)
      !
      IF (rismt%nrzl * rismt%lfft%ngxy * rismt%nsite > 0) THEN
        CALL read_lauerism_x(rismt, (rismt%nrzl * rismt%lfft%ngxy), rismt%nsite, &
                           & rismt%hsgz, ecut, file_base)
      ELSE
        CALL read_lauerism_x(rismt, (rismt%nrzl * rismt%lfft%ngxy), rismt%nsite, &
                           & cdummy, ecut, file_base)
      END IF
      !
      ! ... read hl(laue)
      file_base = TRIM(dirname) // '/3d-rism_hluv_l' // TRIM(ext)
      !
      IF (rismt%nrzl * rismt%lfft%ngxy * rismt%nsite > 0) THEN
        CALL read_lauerism_x(rismt, (rismt%nrzl * rismt%lfft%ngxy), rismt%nsite, &
                           & rismt%hlgz, ecut, file_base)
      ELSE
        CALL read_lauerism_x(rismt, (rismt%nrzl * rismt%lfft%ngxy), rismt%nsite, &
                           & cdummy, ecut, file_base)
      END IF
      !
    END IF
    !
  END SUBROUTINE read_3drism
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_3drism_x(rismt, nr, nsite, zr, ecut, file_base)
    !------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE(rism_type),  INTENT(IN)  :: rismt
    INTEGER,          INTENT(IN)  :: nr
    INTEGER,          INTENT(IN)  :: nsite
    REAL(DP),         INTENT(OUT) :: zr(nr, nsite)
    REAL(DP),         INTENT(IN)  :: ecut
    CHARACTER(LEN=*), INTENT(IN)  :: file_base
    !
    CALL read_3drism_xml(file_base, zr(1:nr, 1:nsite), &
                       & rismt%mp_site%nsite, rismt%mp_site%isite_start, rismt%mp_site%isite_end, &
                       & ecut, rismt%dfft, ionode, rismt%mp_site%inter_sitg_comm)
    !
  END SUBROUTINE read_3drism_x
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_lauerism_x(rismt, nl, nsite, zl, ecut, file_base)
    !------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE(rism_type),  INTENT(IN)  :: rismt
    INTEGER,          INTENT(IN)  :: nl
    INTEGER,          INTENT(IN)  :: nsite
    COMPLEX(DP),      INTENT(OUT) :: zl(nl, nsite)
    REAL(DP),         INTENT(IN)  :: ecut
    CHARACTER(LEN=*), INTENT(IN)  :: file_base
    !
    CALL read_lauerism_xml(file_base, zl(1:nl, 1:nsite), &
                         & rismt%mp_site%nsite, rismt%mp_site%isite_start, rismt%mp_site%isite_end, &
                         & ecut, rismt%lfft, ionode, &
                         & rismt%mp_site%intra_sitg_comm, rismt%mp_site%inter_sitg_comm)
    !
  END SUBROUTINE read_lauerism_x
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_lauedipole_x(rismt, nsite, zd, file_base)
    !------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE(rism_type),  INTENT(IN)  :: rismt
    INTEGER,          INTENT(IN)  :: nsite
    REAL(DP),         INTENT(OUT) :: zd(nsite)
    CHARACTER(LEN=*), INTENT(IN)  :: file_base
    !
    CALL read_lauedipole_xml(file_base, zd(1:nsite), &
                           & rismt%mp_site%nsite, rismt%mp_site%isite_start, rismt%mp_site%isite_end, &
                           & ionode, rismt%mp_site%intra_sitg_comm, rismt%mp_site%inter_sitg_comm)
    !
  END SUBROUTINE read_lauedipole_x
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_lauegxy0_x(rismt, nl, nsite, zl, file_base)
    !------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE(rism_type),  INTENT(IN)  :: rismt
    INTEGER,          INTENT(IN)  :: nl
    INTEGER,          INTENT(IN)  :: nsite
    REAL(DP),         INTENT(OUT) :: zl(nl, nsite)
    CHARACTER(LEN=*), INTENT(IN)  :: file_base
    !
    CALL read_lauegxy0_xml(file_base, zl(1:nl, 1:nsite), &
                         & rismt%mp_site%nsite, rismt%mp_site%isite_start, rismt%mp_site%isite_end, &
                         & rismt%lfft, ionode, &
                         & rismt%mp_site%intra_sitg_comm, rismt%mp_site%inter_sitg_comm)
    !
  END SUBROUTINE read_lauegxy0_x
  !
END MODULE io_rism_xml

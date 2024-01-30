  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Roxana Margine
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE io_eliashberg
  !----------------------------------------------------------------------
  !!
  !! This module contains all the IO part of the superconductivity part of EPW
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_read_aniso_iaxis(itemp)
    !-----------------------------------------------------------------------
    !!
    !! This routine reads from file the anisotropic delta and znorm on the imaginary-axis
    !!
    !! input
    !!
    !! itemp  - temperature point
    !!
    !! SH: Modified to read the full-bandwidth runs' files (Nov 2021).
    !----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nstemp, fsthick, fbw
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nsiw, gap0, agap, wsi, znormi, deltai, &
                              aznormi, adeltai, nkfs, nbndfs, ef0, ekfs, &
                              dosef, wkfs, w0g, ashifti, shifti
    USE constants_epw, ONLY : kelvin2eV, eps6, zero
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE low_lvl,       ONLY : mem_size_eliashberg
    USE division,      ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    !
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! output file name
    CHARACTER(LEN = 256) :: word
    !! character read from file
    !
    INTEGER :: iw
    !! Counter on frequency
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ibnd
    !! Counter on band
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: nks
    !! Number of k points per pool
    INTEGER(8) :: imelt
    !! Required allocation of memory
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: temp
    !! Temperature in K
    REAL(KIND = DP) :: eband
    !! Temporary variable for eigenvalue
    REAL(KIND = DP) :: omega
    !! Temporary variable for frequency
    REAL(KIND = DP) :: weight
    !! Temporaty variable 
    REAL(KIND = DP), ALLOCATABLE :: adeltai_tmp(:, :, :)
    !! Temporary array to collect adeltai from all pools
    REAL(KIND = DP), ALLOCATABLE :: aznormi_tmp(:, :, :)
    !! Temporary array to collect aznormi from all pools
    REAL(KIND = DP), ALLOCATABLE :: ashifti_tmp(:, :, :)
    !! Temporary array to collect ashifti from all pools
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    nks = upper_bnd - lower_bnd + 1
    !
    ! RM - adeltai, aznormi, ashifti are defined per k-points per pool
    !
    IF (fbw) THEN
      ! get the size of required memory for deltai, znormi, shifti,
      ! adeltai, aznormi, ashifti, adeltai_tmp, aznormi_tmp, ashifti_tmp
      imelt = (3 + 3 * nbndfs * nks + 3 * nbndfs * nkfs) * nsiw(itemp)
    ELSE
      ! get the size of required memory for deltai, znormi,
      ! adeltai, aznormi, adeltai_tmp, aznormi_tmp
      imelt = (2 + 2 * nbndfs * nks + 2 * nbndfs * nkfs) * nsiw(itemp)
    ENDIF
    CALL mem_size_eliashberg(2, imelt)
    !
    ALLOCATE(deltai(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating deltai', 1)
    ALLOCATE(znormi(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating znormi', 1)
    !
    ALLOCATE(adeltai(nsiw(itemp), nbndfs, lower_bnd:upper_bnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating adeltai', 1)
    ALLOCATE(aznormi(nsiw(itemp), nbndfs, lower_bnd:upper_bnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating aznormi', 1)
    ALLOCATE(adeltai_tmp(nsiw(itemp), nbndfs, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating adeltai_tmp', 1)
    ALLOCATE(aznormi_tmp(nsiw(itemp), nbndfs, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating aznormi_tmp', 1)    
    deltai(:)  = zero
    znormi(:)  = zero
    adeltai(:, :, :) = zero
    aznormi(:, :, :) = zero
    adeltai_tmp(:, :, :) = zero
    aznormi_tmp(:, :, :) = zero    
    !
    ! SH: allocate extra arrays needed for fbw
    IF (fbw) THEN
      ALLOCATE(shifti(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating shifti', 1)
      ALLOCATE(ashifti(nsiw(itemp), nbndfs, lower_bnd:upper_bnd), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating ashifti', 1)
      ALLOCATE(ashifti_tmp(nsiw(itemp), nbndfs, nkfs), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating ashifti_tmp', 1)
      shifti(:) = zero
      ashifti(:, :, :) = zero
      ashifti_tmp(:, :, :) = zero
    ENDIF
    !
    IF (mpime == ionode_id) THEN
      !
      temp = gtemp(itemp) / kelvin2eV
      ! anisotropic case
      IF (temp < 10.d0) THEN
        WRITE(name1, '(a, a14, f4.2)') TRIM(prefix), '.imag_aniso_00', temp
      ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
        WRITE(name1, '(a, a13, f5.2)') TRIM(prefix), '.imag_aniso_0', temp
      ELSEIF (temp >= 100.d0) THEN
        WRITE(name1, '(a, a12, f6.2)') TRIM(prefix), '.imag_aniso_', temp
      ENDIF
      !
      OPEN(UNIT = iufilgap, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'error opening file ' // name1, iufilgap)
      READ(iufilgap, '(a)') word
      DO iw = 1, nsiw(itemp) ! loop over omega
        DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              ! SH: read fbw entries
              IF (fbw) THEN
                READ(iufilgap, '(5ES20.10)') omega, eband, &
                  aznormi_tmp(iw, ibnd, ik), adeltai_tmp(iw, ibnd, ik), ashifti_tmp(iw, ibnd, ik)
              ELSE
                READ(iufilgap, '(4ES20.10)') omega, eband, &
                  aznormi_tmp(iw, ibnd, ik), adeltai_tmp(iw, ibnd, ik)
              ENDIF
              IF (iw == 1) &
                agap(ibnd, ik) = adeltai_tmp(1, ibnd, ik)
            ENDIF
          ENDDO ! ibnd
        ENDDO ! ik
        IF (ABS(wsi(iw) - omega) > eps6) &
          CALL errore('eliashberg_read_aniso_iaxis', 'temperature not the same with the input', 1)
      ENDDO ! iw
      CLOSE(iufilgap)
      !
      DO iw = 1, nsiw(itemp) ! loop over omega
        DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              weight = 0.5d0 * wkfs(ik) * w0g(ibnd, ik) / dosef
              znormi(iw) = znormi(iw) + weight * aznormi_tmp(iw, ibnd, ik)
              deltai(iw) = deltai(iw) + weight * adeltai_tmp(iw, ibnd, ik)
              ! SH: shifti is present only for fbw runs
              IF (fbw) shifti(iw) = shifti(iw) + weight * ashifti_tmp(iw, ibnd, ik)
            ENDIF
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! iw
      gap0 = deltai(1)
      !
      CALL gap_FS(itemp)
      !
    ENDIF
    CALL mp_bcast(gap0, ionode_id, inter_pool_comm)
    CALL mp_bcast(agap, ionode_id, inter_pool_comm)
    CALL mp_bcast(deltai, ionode_id, inter_pool_comm)
    CALL mp_bcast(znormi, ionode_id, inter_pool_comm)
    CALL mp_bcast(adeltai_tmp, ionode_id, inter_pool_comm)
    CALL mp_bcast(aznormi_tmp, ionode_id, inter_pool_comm)
    ! SH: extra arrays for fbw runs
    IF (fbw) THEN
      CALL mp_bcast(shifti, ionode_id, inter_pool_comm)
      CALL mp_bcast(ashifti_tmp, ionode_id, inter_pool_comm)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    adeltai(:, :, lower_bnd:upper_bnd)  = adeltai_tmp(:, :, lower_bnd:upper_bnd)
    aznormi(:, :, lower_bnd:upper_bnd)  = aznormi_tmp(:, :, lower_bnd:upper_bnd)
    !
    IF (fbw) THEN
      ashifti(:, :, lower_bnd:upper_bnd) = ashifti_tmp(:, :, lower_bnd:upper_bnd) 
    ENDIF
    !
    DEALLOCATE(deltai, STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error deallocating deltai', 1)
    DEALLOCATE(znormi, STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error deallocating znormi', 1)    
    DEALLOCATE(adeltai_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error deallocating adeltai_tmp', 1)
    DEALLOCATE(aznormi_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error deallocating aznormi_tmp', 1)
    IF (fbw) THEN
      DEALLOCATE(shifti, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error deallocating shifti', 1)      
      DEALLOCATE(ashifti_tmp, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error deallocating ashifti_tmp', 1)
    ENDIF
    !
    IF (fbw) THEN
      ! remove memory allocated for deltai, znormi, shifti, adeltai_tmp, aznormi_tmp, ashifti_tmp
      imelt = (3 + 3 * nbndfs * nkfs) * nsiw(itemp)
    ELSE
      ! remove memory allocated for deltai, znormi, adeltai_tmp, aznormi_tmp
      imelt = (2 + 2 * nbndfs * nkfs) * nsiw(itemp)      
    ENDIF
    CALL mem_size_eliashberg(2, -imelt)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE eliashberg_read_aniso_iaxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_write_iaxis(itemp)
    !-----------------------------------------------------------------------
    !!
    !! This routine writes to files results from the solutions of the Eliashberg equations
    !! on the imaginary-axis
    !!
    !! SH: Modified to write the full-bandwidth runs' files, and re-ordered
    !!       "deltai, ..." arrays' indices for efficiency (Nov 2021).
    !!
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : fsthick, laniso, liso, fbw
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nkfs, nbndfs, ef0, ekfs, nsiw, wsi, agap, &
                              deltai, znormi, shifti, adeltai, aznormi, ashifti
    USE constants_epw, ONLY : kelvin2eV, zero
    USE division,      ONLY : fkbounds
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE io_global,     ONLY : ionode_id
    USE mp_world,      ONLY : mpime
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter for temperature
    !
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! output file name
    CHARACTER(LEN = 256) :: cname
    !! character in output file name
    !
    INTEGER :: iw
    !! Counter on frequency imag-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: nks
    !! Number of k points per pool
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios
    !! IO error message
    !
    REAL(KIND = DP) :: temp
    !! Temperature in K
    REAL(KIND = DP), ALLOCATABLE :: adeltai_tmp(:, :, :)
    !! Temporary array to collect adeltai from all pools
    REAL(KIND = DP), ALLOCATABLE :: aznormi_tmp(:, :, :)
    !! Temporary array to collect aznormi from all pools
    REAL(KIND = DP), ALLOCATABLE :: ashifti_tmp(:, :, :)
    !! Temporary array to collect ashifti from all pools
    !
    temp = gtemp(itemp) / kelvin2eV
    !
    cname = 'imag'
    !
    IF (laniso) THEN
      !
      CALL fkbounds(nkfs, lower_bnd, upper_bnd)
      !
      nks = upper_bnd - lower_bnd + 1
      !
      ALLOCATE(adeltai_tmp(nsiw(itemp), nbndfs, nkfs), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_write_iaxis', 'Error allocating adeltai_tmp', 1)
      ALLOCATE(aznormi_tmp(nsiw(itemp), nbndfs, nkfs), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_write_iaxis', 'Error allocating aznormi_tmp', 1)
      adeltai_tmp(:, :, :) = zero
      aznormi_tmp(:, :, :) = zero
      !
      adeltai_tmp(:, :, lower_bnd:upper_bnd) = adeltai(:, :, lower_bnd:upper_bnd)
      aznormi_tmp(:, :, lower_bnd:upper_bnd) = aznormi(:, :, lower_bnd:upper_bnd)
      !
      IF (fbw) THEN 
        ALLOCATE(ashifti_tmp(nsiw(itemp), nbndfs, nkfs), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_write_iaxis', 'Error allocating ashifti_tmp', 1)
        ashifti_tmp(:, :, :) = zero
        ashifti_tmp(:, :, lower_bnd:upper_bnd) = ashifti(:, :, lower_bnd:upper_bnd)
      ENDIF     
      !   
      ! collect contributions from all pools
      CALL mp_sum(adeltai_tmp, inter_pool_comm)
      CALL mp_sum(aznormi_tmp, inter_pool_comm)
      IF (fbw) CALL mp_sum(ashifti_tmp, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
      IF (mpime == ionode_id) THEN
        ! 
        IF (temp < 10.d0) THEN
          WRITE(name1, '(a, a1, a4, a9, f4.2)') TRIM(prefix), '.', cname, '_aniso_00', temp
        ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
          WRITE(name1, '(a, a1, a4, a8, f5.2)') TRIM(prefix), '.', cname, '_aniso_0', temp
        ELSEIF (temp >= 100.d0) THEN
          WRITE(name1, '(a, a1, a4, a7, f6.2)') TRIM(prefix), '.', cname, '_aniso_', temp
        ENDIF
        OPEN(UNIT = iufilgap, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('eliashberg_write_iaxis', 'error opening file ' // name1, iufilgap)
        !
        ! SH: write the header for fbw runs
        IF (fbw) THEN
          WRITE(iufilgap,'(6a20)') '#        w [eV]', 'Enk-Ef [eV]', 'znorm(w)', 'delta(w) [eV]', 'shift(w) [eV]'
        ELSE
          WRITE(iufilgap,'(5a20)') '#        w [eV]', 'Enk-Ef [eV]', 'znorm(w)', 'delta(w) [eV]'
        ENDIF
        DO iw = 1, nsiw(itemp) ! loop over omega
          DO ik = 1, nkfs
            DO ibnd = 1, nbndfs
              IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
                ! SH: write the entries for fbw runs
                IF (fbw) THEN
                  WRITE(iufilgap, '(5ES20.10)') wsi(iw), ekfs(ibnd, ik) - ef0, &
                    aznormi_tmp(iw, ibnd, ik), adeltai_tmp(iw, ibnd, ik), ashifti_tmp(iw, ibnd, ik)
                ELSE
                  WRITE(iufilgap, '(4ES20.10)') wsi(iw), ekfs(ibnd, ik) - ef0, &
                    aznormi_tmp(iw, ibnd, ik), adeltai_tmp(iw, ibnd, ik)
                ENDIF
                IF (iw == 1) &
                  agap(ibnd, ik) = adeltai_tmp(iw, ibnd, ik)
              ENDIF
            ENDDO ! ibnd
          ENDDO ! ik
        ENDDO ! iw
        CLOSE(iufilgap)
        !
        CALL gap_distribution_FS(itemp, cname)
        !
        CALL gap_FS(itemp)
        !
      ENDIF ! mpime
      !
      CALL mp_bcast(agap, ionode_id, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
      DEALLOCATE(adeltai_tmp, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_write_iaxis', 'Error deallocating adeltai_tmp', 1)
      DEALLOCATE(aznormi_tmp, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_write_iaxis', 'Error deallocating aznormi_tmp', 1)
      IF (fbw) THEN 
        DEALLOCATE(ashifti_tmp, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_write_iaxis', 'Error deallocating ashifti_tmp', 1)
      ENDIF
      !      
    ENDIF ! laniso
    !
    ! isotropic case
    IF (liso) THEN
      IF (temp < 10.d0) THEN
        WRITE(name1, '(a, a1, a4, a7, f4.2)') TRIM(prefix), '.', cname, '_iso_00', temp
      ELSEIF (temp >= 10.d0 .AND. temp < 100.d0 ) THEN
        WRITE(name1, '(a, a1, a4, a6, f5.2)') TRIM(prefix), '.', cname, '_iso_0', temp
      ELSEIF (temp >= 100.d0) THEN
        WRITE(name1, '(a, a1, a4, a5, f6.2)') TRIM(prefix), '.', cname, '_iso_', temp
      ENDIF
      OPEN(UNIT = iufilgap, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('eliashberg_write_iaxis', 'error opening file ' // name1, iufilgap)
      !
      ! SH: write the header for fbw runs
      IF (fbw) THEN
        WRITE(iufilgap, '(4a20)') 'w [eV]', 'znorm(w)', 'delta(w) [eV]', 'shift(w) [eV]'
      ELSE
        WRITE(iufilgap, '(3a20)') 'w [eV]', 'znorm(w)', 'delta(w) [eV]'
      ENDIF
      DO iw = 1, nsiw(itemp) ! loop over omega
        ! SH: write the entries for fbw runs
        IF (fbw) THEN
          WRITE(iufilgap, '(4ES20.10)') wsi(iw), znormi(iw), deltai(iw), shifti(iw)
        ELSE
          WRITE(iufilgap, '(3ES20.10)') wsi(iw), znormi(iw), deltai(iw)
        ENDIF
      ENDDO
      CLOSE(iufilgap)
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE eliashberg_write_iaxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_write_raxis(itemp, cname)
    !-----------------------------------------------------------------------
    !!
    !! This routine writes to files results from the solutions of the Eliashberg
    !! equations on the real-axis
    !!
    !! SH: Modified to write the full-bandwidth runs' files (Nov 2021).
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nqstep, fsthick, laniso, liso, fbw
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nsw, ws, gap0, agap, delta, znorm, adelta, aznorm, &
                              nkfs, nbndfs, ef0, ekfs, shift, ashift
    USE constants_epw, ONLY : kelvin2eV, zero, czero
    USE division,      ONLY : fkbounds
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE io_global,     ONLY : ionode_id
    USE low_lvl,       ONLY : mem_size_eliashberg    
    USE mp_world,      ONLY : mpime
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter for temperature
    CHARACTER(LEN = 256), INTENT(in) :: cname
    !! character in output file name
    !
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! output file name
    !
    LOGICAL :: lgap
    !! True if gap found
    !
    INTEGER :: iw
    !! Counter on frequency real-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: nks
    !! Number of k points per pool    
    INTEGER(8) :: imelt
    !! Required allocation of memory    
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status    
    !
    REAL(KIND = DP) :: temp
    !! Temperature in K
    REAL(KIND = DP) :: var1, var2, var3, var4
    !! Temporary working variables
    !
    COMPLEX(KIND = DP), ALLOCATABLE :: adelta_tmp(:, :, :)
    !! Temporary array to collect adelta from all pools
    COMPLEX(KIND = DP), ALLOCATABLE :: aznorm_tmp(:, :, :)
    !! Temporary array to collect aznorm from all pools
    COMPLEX(KIND = DP), ALLOCATABLE :: ashift_tmp(:, :, :)
    !! Temporary array to collect ashift from all pools
    !
    temp = gtemp(itemp) / kelvin2eV
    !
    IF (laniso) THEN
      !
      CALL fkbounds(nkfs, lower_bnd, upper_bnd)
      !
      nks = upper_bnd - lower_bnd + 1
      !
      agap(:, :) = zero
      DO ik = lower_bnd, upper_bnd
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
            lgap = .TRUE.
            ! DO iw = 1, nsw
            DO iw = 1, nsw - 1   ! FG: this change is to prevent segfault in ws(iw+1) and adelta(*,*,iw+1)
              var1 = REAL(adelta(iw, ibnd, ik))
              var2 = REAL(adelta(iw + 1, ibnd, ik))
              var3 = var1 - ws(iw)
              var4 = var2 - ws(iw + 1)
              IF (lgap .AND. iw < nqstep .AND. var1 > 0.d0 .AND. var2 > 0.d0 .AND. var3 * var4 < 0.d0) THEN
                agap(ibnd, ik) = (var3 * ws(iw + 1) - var4 * ws(iw)) / (var3 - var4)
                lgap = .FALSE.
              ENDIF
            ENDDO ! iw
            IF (lgap) &
              agap(ibnd, ik) = REAL(adelta(1, ibnd, ik))
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
      !
      ! collect contributions from all pools
      CALL mp_sum(agap, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
      CALL gap_distribution_FS(itemp, cname)
      !
      IF (iverbosity == 2) THEN 
        IF (fbw) THEN
          ! get the size of required memory for adelta_tmp, aznorm_tmp, ashift_tmp
          imelt = 2 * 3 * nbndfs * nkfs * nsw
        ELSE
          ! get the size of required memory for adelta_tmp, aznorm_tmp, ashift_tmp
          imelt = 2 * 2 * nbndfs * nkfs * nsw                
        ENDIF
        CALL mem_size_eliashberg(2, imelt)
        !      
        ALLOCATE(adelta_tmp(nsw, nbndfs, nkfs), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_write_raxis', 'Error allocating adelta_tmp', 1)      
        ALLOCATE(aznorm_tmp(nsw, nbndfs, nkfs), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_write_raxis', 'Error allocating aznorm_tmp', 1)
        adelta_tmp(:, :, :) = czero
        adelta_tmp(:, :, lower_bnd:upper_bnd) = adelta(:, :, lower_bnd:upper_bnd)
        aznorm_tmp(:, :, :) = czero
        aznorm_tmp(:, :, lower_bnd:upper_bnd) = aznorm(:, :, lower_bnd:upper_bnd)
        !
        IF (fbw) THEN
          ALLOCATE(ashift_tmp(nsw, nbndfs, nkfs), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_write_raxis', 'Error allocating ashift_tmp', 1)
          ashift_tmp(:, :, :) = czero
          ashift_tmp(:, :, lower_bnd:upper_bnd) = ashift(:, :, lower_bnd:upper_bnd)
        ENDIF
        !
        ! collect contributions from all pools
        CALL mp_sum(adelta_tmp, inter_pool_comm)
        CALL mp_sum(aznorm_tmp, inter_pool_comm)
        IF (fbw) CALL mp_sum(ashift_tmp, inter_pool_comm)
        CALL mp_barrier(inter_pool_comm)
        ! 
        IF (mpime == ionode_id) THEN
          IF (temp < 10.d0) THEN
            WRITE(name1, '(a, a1, a4, a9, f4.2)') TRIM(prefix), '.', cname, '_aniso_00', temp
          ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
            WRITE(name1, '(a, a1, a4, a8, f5.2)') TRIM(prefix), '.', cname, '_aniso_0', temp
          ELSEIF (temp >= 100.d0) THEN
            WRITE(name1, '(a, a1, a4, a7, f6.2)') TRIM(prefix), '.', cname, '_aniso_', temp
          ENDIF
          OPEN(UNIT = iufilgap, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
          IF (ios /= 0) CALL errore('eliashberg_write_raxis', 'error opening file ' // name1, iufilgap)
          ! SH: write the header for fbw runs
          IF (fbw) THEN
            WRITE(iufilgap, '(8a20)') '#        w [eV]', 'Enk-Ef [eV]', 'Re[znorm(w)]', 'Im[znorm(w)]', &
              'Re[delta(w)] [eV]', 'Im[delta(w)] [eV]', 'Re[shift(w)] [eV]', 'Im[shift(w)] [eV]'
          ELSE
            WRITE(iufilgap, '(6a20)') '#        w [eV]', 'Enk-Ef [eV]', 'Re[znorm(w)]', 'Im[znorm(w)]', &
              'Re[delta(w)] [eV]', 'Im[delta(w)] [eV]'
          ENDIF
          !
          DO ik = 1, nkfs
            DO ibnd = 1, nbndfs
              IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
                DO iw = 1, nsw
                  ! SH: write the entries for fbw runs
                  IF (fbw) THEN
                    WRITE(iufilgap, '(8ES20.10)') ws(iw), ekfs(ibnd, ik) - ef0, &
                      REAL(aznorm_tmp(iw, ibnd, ik)), AIMAG(aznorm_tmp(iw, ibnd, ik)), &
                      REAL(adelta_tmp(iw, ibnd, ik)), AIMAG(adelta_tmp(iw, ibnd, ik)), &
                      REAL(ashift_tmp(iw, ibnd, ik)), AIMAG(ashift_tmp(iw, ibnd, ik))
                  ELSE
                    WRITE(iufilgap, '(6ES20.10)') ws(iw), ekfs(ibnd, ik) - ef0, &
                      REAL(aznorm_tmp(iw, ibnd, ik)), AIMAG(aznorm_tmp(iw, ibnd, ik)), &
                      REAL(adelta_tmp(iw, ibnd, ik)), AIMAG(adelta_tmp(iw, ibnd, ik))
                  ENDIF
                ENDDO ! iw
              ENDIF
            ENDDO ! ibnd
          ENDDO ! ik
          CLOSE(iufilgap)
          !
        ENDIF ! mpime
        ! 
        DEALLOCATE(adelta_tmp, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_write_raxis', 'Error deallocating adelta_tmp', 1)
        DEALLOCATE(aznorm_tmp, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_write_raxis', 'Error deallocating aznorm_tmp', 1)
        IF (fbw) THEN
          DEALLOCATE(ashift_tmp, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_write_raxis', 'Error deallocating ashift_tmp', 1)
        ENDIF
        !
        IF (fbw) THEN
          ! get the size of required memory for adelta_tmp, aznorm_tmp, ashift_tmp
          imelt = 2 * 3 * nbndfs * nkfs * nsw
        ELSE
          ! get the size of required memory for adelta_tmp, aznorm_tmp, ashift_tmp
          imelt = 2 * 2 * nbndfs * nkfs * nsw
        ENDIF
        CALL mem_size_eliashberg(2, -imelt)
        !  
      ENDIF ! iverbosity 
      !
    ENDIF ! laniso
    !
    ! isotropic case
    IF (liso) THEN
      IF (temp < 10.d0) THEN
        WRITE(name1, '(a, a1, a4, a7, f4.2)') TRIM(prefix), '.', cname, '_iso_00', temp
      ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
        WRITE(name1, '(a, a1, a4, a6, f5.2)') TRIM(prefix), '.', cname, '_iso_0', temp
      ELSEIF (temp >= 100.d0) THEN
        WRITE(name1, '(a, a1, a4, a5, f6.2)') TRIM(prefix), '.', cname, '_iso_', temp
      ENDIF
      OPEN(UNIT = iufilgap, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('eliashberg_write_raxis', 'error opening file ' // name1, iufilgap)
      ! SH: write the header for fbw runs
      IF (fbw) THEN
        WRITE(iufilgap, '(7a20)') 'w [eV]', 'Re[znorm(w)]', 'Im[znorm(w)]', 'Re[delta(w)] [eV]', 'Im[delta(w)] [eV]', &
          'Re[shift(w)] [eV]', 'Im[shift(w)] [eV]'
      ELSE
        WRITE(iufilgap, '(5a20)') 'w [eV]', 'Re[znorm(w)]', 'Im[znorm(w)]', 'Re[delta(w)] [eV]', 'Im[delta(w)] [eV]'
      ENDIF
      !
      lgap = .TRUE.
      ! DO iw = 1, nsw
      DO iw = 1, nsw - 1   ! this change is to prevent segfault in delta(iw+1) and ws(iw+1)
        var1 = REAL(delta(iw))
        var2 = REAL(delta(iw + 1))
        var3 = var1 - ws(iw)
        var4 = var2 - ws(iw + 1)
        IF (lgap .AND. iw < nqstep .AND. var1 > 0.d0 .AND. var2 > 0.d0 .AND. var3 * var4 < 0.d0) THEN
            gap0 = (var3 * ws(iw + 1) - var4 * ws(iw)) / (var3 - var4)
          lgap = .FALSE.
        ENDIF
        ! SH: write the entries for fbw runs
        IF (fbw) THEN
          WRITE(iufilgap, '(7ES20.10)') ws(iw), REAL(znorm(iw)), AIMAG(znorm(iw)), &
             REAL(delta(iw)), AIMAG(delta(iw)), REAL(shift(iw)), AIMAG(shift(iw))
        ELSE
          WRITE(iufilgap, '(5ES20.10)') ws(iw), REAL(znorm(iw)), AIMAG(znorm(iw)), &
             REAL(delta(iw)), AIMAG(delta(iw))
        ENDIF
      ENDDO ! iw
      CLOSE(iufilgap)
      IF (lgap) &
        gap0 = REAL(delta(1))
    ENDIF ! liso
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE eliashberg_write_raxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE write_phdos(nrr_q, irvec_q, ndegen_q, nrws, rws)
    !-----------------------------------------------------------------------
    !!
    !! SH: Writes "phdos" files for phonon density of states (Nov 2021).
    !! RM: Updated (Nov 2021).
    !!
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iudosfil, iunselecq, iufilfreq
    USE io_files,      ONLY : prefix, tmp_dir
    USE ions_base,     ONLY : nat
    USE epwcom,        ONLY : nqstep, nqsmear, degaussq, delta_qsmear, eps_acustic, &
                              lifc
    USE constants_epw, ONLY : zero, ryd2mev
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id, npool
    USE mp_world,      ONLY : mpime, world_comm
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE io_global,     ONLY : stdout, ionode_id
    !
    USE modes,         ONLY : nmodes
    USE elph2,         ONLY : wf, wqf, nqtotf, xqf
    USE wan2bloch,     ONLY : dynwan2bloch, dynifc2blochf
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrr_q
    !! number of WS points
    INTEGER, INTENT(in) :: irvec_q(3, nrr_q)
    !! INTEGER components of the ir-th Wigner-Seitz grid point for phonons
    INTEGER, INTENT(in) :: ndegen_q(nrr_q, nat, nat)
    !! Wigner-Seitz weights for the phonon grid that depend on atomic positions $R + \tau(nb) - \tau(na)$
    INTEGER, INTENT(in) :: nrws
    !! Number of Wigner-Size real space vectors
    REAL(KIND = DP), INTENT(in) :: rws(0:3, nrws)
    !! Real space Wigner-Seitz vector
    !
    ! Local variables
    CHARACTER(LEN = 256) :: name1, filfreq, dirname
    !! file names
    !
    INTEGER :: ismear
    !! Counter on smearings
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios
    !! IO error message
    INTEGER :: iq
    !! Counter on q-points
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: iwph
    !! Counter over frequncy
    !
    REAL(KIND = DP) :: xxq(3)
    !! Current q-point
    REAL(KIND = DP) :: weightq
    !! factors in lambda_eph and a2f
    REAL(KIND = DP) :: sigma
    !! smearing in delta function
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP) :: w2_tmp(nmodes)
    !! Interpolated phonon frequency
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss: an approximation to the delta function
    REAL(KIND = DP), ALLOCATABLE :: ww(:)
    !! Current frequency 
    REAL(KIND = DP), ALLOCATABLE :: phdos(:, :)
    !! Phonon density of states for different ismear
    REAL(KIND = DP), ALLOCATABLE :: phdos_modeproj(:, :)
    !! Phonon density of states  projected over modes for different ismear
    COMPLEX(KIND = DP) :: uf_tmp(nmodes, nmodes)
    !! Rotation matrix for phonons
    !
    ! HM: The phonon frequencies should have been calculated only for some points: From iqq = iq_restart to iqq = totq.
    !     Here we calculate the phonon frequencies for all the q points to calculate phdos and a2F.
    !
    DO iq = 1, nqtotf
      xxq = xqf(:, iq)
      IF (.NOT. lifc) THEN
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq, uf_tmp, w2_tmp)
      ELSE
        CALL dynifc2blochf(nmodes, rws, nrws, xxq, uf_tmp, w2_tmp)
      ENDIF
      !
      DO imode = 1, nmodes
        !
        ! wf are the interpolated eigenfrequencies (omega on fine grid)
        IF (w2_tmp(imode) > 0.d0) THEN
          wf(imode, iq) =  SQRT(ABS(w2_tmp(imode)))
        ELSE
          wf(imode, iq) = -SQRT(ABS(w2_tmp(imode)))
        ENDIF
      ENDDO
    ENDDO
    !
    ! Starting the calculation of phdos
    !
    IF (mpime == ionode_id) THEN
      !
      ALLOCATE(ww(nqstep), STAT = ierr)
      IF (ierr /= 0) CALL errore('write_phdos', 'Error allocating ww', 1)
      ww(:) = zero
      !
      dw = 1.1d0 * MAXVAL(wf(:, :)) /  DBLE(nqstep) ! increase by 10%
      DO iwph = 1, nqstep  
        ww(iwph) = DBLE(iwph) * dw
      ENDDO
      !
      ALLOCATE(phdos(nqstep, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('write_phdos', 'Error allocating phdos', 1)
      ALLOCATE(phdos_modeproj(nmodes, nqstep), STAT = ierr)
      IF (ierr /= 0) CALL errore('write_phdos', 'Error allocating phdos_modeproj', 1)
      phdos(:, :) = zero
      phdos_modeproj(:, :) = zero
      !
      DO ismear = 1, nqsmear
        sigma = degaussq + (ismear - 1) * delta_qsmear
        DO iq = 1, nqtotf
          DO imode = 1, nmodes
            IF (wf(imode, iq) > eps_acustic) THEN
              DO iwph = 1, nqstep
                weightq  = w0gauss((ww(iwph) - wf(imode, iq)) / sigma, 0) / sigma
                phdos(iwph, ismear) = phdos(iwph, ismear) + wqf(iq) * weightq
                IF (ismear == 1) &
                  phdos_modeproj(imode, iwph) = phdos_modeproj(imode, iwph) + wqf(iq) * weightq
              ENDDO ! iwph
            ENDIF ! wf
          ENDDO ! imode
        ENDDO ! iq
      ENDDO ! ismear
      !
      name1 = TRIM(prefix) // '.phdos'
      OPEN(UNIT = iudosfil, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('write_phdos', 'error opening file ' // name1, iudosfil)
      !
      DO ismear = 1, nqsmear
        IF (ismear == nqsmear) &
          WRITE(iudosfil, '(" w[meV] phdos[states/meV] for ", i4, " smearing values")') ismear
        DO iwph = 1, nqstep
          IF (ismear == nqsmear) &
            WRITE(iudosfil, '(f12.7, 15f15.7)') ww(iwph) * ryd2mev, phdos(iwph, :) / ryd2mev
        ENDDO
      ENDDO
      !
      CLOSE(iudosfil)
      !
      name1 = TRIM(prefix) // '.phdos_proj'
      OPEN(UNIT = iudosfil, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('write_phdos', 'error opening file ' // name1, iudosfil)
      WRITE(iudosfil, '("w[meV] phdos[states/meV] phdos_modeproj[states/meV]")')
      DO iwph = 1, nqstep
        ! wsph in meV (from eV) and phdos in states/meV (from states/eV)
        WRITE(iudosfil, '(f12.7, 100f15.7)') ww(iwph) * ryd2mev, & 
              phdos(iwph, 1) / ryd2mev, phdos_modeproj(:, iwph) / ryd2mev
      ENDDO
      CLOSE(iudosfil)
      !
      DEALLOCATE(ww, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_phdos', 'Error deallocating ww', 1)
      DEALLOCATE(phdos, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_phdos', 'Error deallocating phdos', 1)
      DEALLOCATE(phdos_modeproj, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_phdos', 'Error deallocating phdos_modeproj', 1)
    ENDIF
    !
    WRITE(stdout, '(/5x, a/)') 'Finish writing phdos files ' // TRIM(prefix) // '.phdos' &
      // ' and ' // TRIM(prefix) // '.phdos_proj'
    ! 
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE write_phdos
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE write_dos(eferm, nele)
    !-----------------------------------------------------------------------
    !!
    !! SH: This routine writes electronic "dos" file; it will be 
    !!       used for full-bandwidth iso run (Nov 2021).
    !!
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufildos
    USE io_files,      ONLY : prefix, tmp_dir
    USE epwcom,        ONLY : nbndsub, ngaussw, degaussw, fsthick, dos_del
    USE constants_epw, ONLY : ryd2ev, zero
    USE elph2,         ONLY : etf, wkf, nkqf
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP) :: eferm
    !! Fermi energy
    REAL(KIND = DP) :: nele
    !! pre-calculated number of electrons in Fermi window
    !
    ! Local variables
    CHARACTER(LEN = 256) :: fildos
    !! Name of Fermi ene. window dos file
    !
    INTEGER :: i
    !! counter over energy intervals
    INTEGER :: wnum
    !! number of energy intervals in Fermi window
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: wene
    !! energies in Fermi window
    REAL(KIND = DP) :: wval
    !! dos value in fermi window
    REAL(KIND = DP) :: wele
    !! no of electrons in fermi window
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! Function to compute the density of states at the Fermi level
    !
    fildos = TRIM(tmp_dir) // '/' // TRIM(prefix) // '.dos'
    !
    wnum = NINT(2.d0 * fsthick * ryd2ev / dos_del)
    wele = zero
    !
    OPEN(UNIT = iufildos, FILE = fildos, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
    IF (ios /= 0) CALL errore('write_dos', 'error opening file ' // fildos, iufildos)
    WRITE(iufildos, '(a, ES20.10, a, ES20.10, a)') '# EFermi[eV]    ', eferm * ryd2ev ,'   dos_EFermi[eV^-1] ', &
      dos_ef(ngaussw, degaussw, eferm, etf, wkf, nkqf, nbndsub) / ryd2ev
    WRITE(iufildos, '(a, ES20.10, a, ES20.10)') '# FermiWind[eV] ', fsthick * ryd2ev, '   Nr_electrons      ', nele
    WRITE(iufildos, '(a, ES20.10, a, i8)')     '# dos_del[eV]   ', dos_del,          '   Nr_bins           ', wnum
    WRITE(iufildos, '(a)') '#            E [eV]           dos[state/eV]       Int dos[#]'
    !
    DO i = 1, wnum
      wene = (eferm - fsthick) + DBLE(i - 1) * dos_del / ryd2ev
      wval = dos_ef(ngaussw, degaussw, wene, etf, wkf, nkqf, nbndsub)
      wele = wele + wval / ryd2ev * dos_del
      ! total DOS is being written (not per spin!) for compatibility with QE format
      WRITE(iufildos, '(5x, 3ES20.10)') wene * ryd2ev, wval / ryd2ev , wele
    ENDDO
    !
    CLOSE(iufildos)
    WRITE(stdout, '(/5x, a/)') 'Finish writing dos file ' // TRIM(prefix) // '.dos'
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE write_dos
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_dos()
    !-----------------------------------------------------------------------
    !!
    !! SH: Read the electronic dos in Fermi window from "dos" file (Nov 2021).
    !!
    USE eliashbergcom, ONLY : en, dosen, ndos, ef0, dosef
    USE constants_epw, ONLY : zero, two
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufildos
    USE io_global,     ONLY : stdout
    USE io_files,      ONLY : prefix, tmp_dir
    USE epwcom,        ONLY : dos_del, fila2f
    !
    IMPLICIT NONE
    !
    ! Local variables
    CHARACTER(LEN = 256) :: fildos
    !! Output file name
    CHARACTER(LEN = 256) :: tmpst0, tmpst1
    !! Temporary string variable
    !
    LOGICAL :: file_exists
    !! Status of file
    !
    INTEGER :: idos
    !! Counter over ndos
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: neband
    !! Integrated dos (Nr of electrons)
    !
    fildos = TRIM(tmp_dir) // '/' // TRIM(prefix) // '.dos'
    !
    WRITE(stdout, '(/5x, a/)') 'Start reading dos file ' // TRIM(prefix) // '.dos'
    WRITE(stdout, '(a)') ' '
    !
    INQUIRE(FILE=fildos, EXIST=file_exists)
    IF (.NOT. file_exists) CALL errore('read_dos', 'error opening file ' // fildos, 1)
    !
    OPEN(UNIT = iufildos, FILE = fildos, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
    !
    READ(iufildos, '(a)') tmpst0
    READ(iufildos, '(a)') tmpst1
    READ(iufildos, *) tmpst1, tmpst1, dos_del, tmpst1, ndos
    READ(iufildos, '(a)') tmpst1
    !! if it's a run from a2f file, Ef and dos(Ef) are needed!
    IF (fila2f /= ' ') THEN
      READ(tmpst0, *) tmpst1, tmpst1, ef0, tmpst1, dosef
      dosef = dosef / two
    ENDIF
    !
    ALLOCATE(en(ndos), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_dos', 'Error allocating en', 1)
    ALLOCATE(dosen(ndos), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_dos', 'Error allocating dosen', 1)
    en(:)    = zero
    dosen(:) = zero
    neband   = zero
    !
    DO idos = 1, ndos
      READ(iufildos,  *) en(idos), dosen(idos), neband
    ENDDO
    dosen(:) = dosen(:) / two ! switch to per spin
    !
    CLOSE(iufildos)
    !
    WRITE(stdout, '(/5x, a/)') 'Finish reading dos file ' // TRIM(prefix) // '.dos'
    WRITE(stdout, '(a)') ' '
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE read_dos
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_a2f()
    !-----------------------------------------------------------------------
    !!
    !! Read the eliashberg spectral function from fila2f
    !!
    !! SH: Modified with removing "a2f_iso file" (Nov 2021).
    !!
    USE epwcom,        ONLY : nqstep, fila2f, laniso
    USE eliashbergcom, ONLY : wsphmax, wsph, dwsph, a2f_iso
    USE constants_epw, ONLY : zero
    USE io_var,        ONLY : iua2ffil
    USE io_global,     ONLY : ionode_id, stdout
    USE io_files,      ONLY : prefix
    USE mp_global,     ONLY : npool, inter_pool_comm
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,      ONLY : mpime
    !
    IMPLICIT NONE
    !
    INTEGER :: iwph
    !! Counter over frequncy
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    !
    ALLOCATE(a2f_iso(nqstep), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_a2f', 'Error allocating a2f_iso', 1)
    ALLOCATE(wsph(nqstep), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_a2f', 'Error allocating wsph', 1)
    a2f_iso(:) = zero
    wsph(:) = zero
    !
    IF (fila2f == ' ') WRITE(fila2f, '(a, a4)') TRIM(prefix), '.a2f'
    IF (mpime == ionode_id) THEN
      OPEN(UNIT = iua2ffil, FILE = fila2f, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_a2f', 'error opening file ' // fila2f, iua2ffil)
      !
      READ(iua2ffil, *)
      DO iwph = 1, nqstep
        READ(iua2ffil, *) wsph(iwph), a2f_iso(iwph) ! freq from meV to eV
        wsph(iwph) = wsph(iwph) / 1000.d0
      ENDDO
      wsphmax = wsph(nqstep)
      !dwsph = wsphmax / DBLE(nqstep - 1)
      dwsph = wsphmax / DBLE(nqstep)
      CLOSE(iua2ffil)
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    CALL mp_bcast(a2f_iso, ionode_id, inter_pool_comm)
    CALL mp_bcast(wsph, ionode_id, inter_pool_comm)
    CALL mp_bcast(wsphmax, ionode_id, inter_pool_comm)
    CALL mp_bcast(dwsph, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout, '(/5x, a/)') 'Finish reading a2f file'
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE read_a2f
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_frequencies()
    !-----------------------------------------------------------------------
    !
    ! read the frequencies obtained from a previous epw run
    !
    USE io_global, ONLY : stdout, ionode_id
    USE io_var,    ONLY : iufilfreq, iunselecq
    USE io_files,  ONLY : prefix, tmp_dir
    USE modes,     ONLY : nmodes
    USE elph2,     ONLY : nqtotf, wf, wqf, xqf, totq, selecq
    USE epwcom,    ONLY : nqf1, nqf2, nqf3, nqstep
    USE eliashbergcom, ONLY : wsphmax, dwsph, wsph
    USE constants_epw, ONLY : ryd2ev, zero
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256) :: filfreq
    !! file name
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory where ikmap/egnv/freq/ephmat files are saved
    !
    INTEGER :: iq
    !! Counter on q points
    INTEGER :: imode
    !! Counter on modes
    INTEGER :: iwph
    !! Counter over frequncy
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    INTEGER :: iqq
    !! Q-point index from selecq.fmt window
    INTEGER :: nqf1_, nqf2_, nqf3_
    !! Temporary variable for number of q-points along each direction
    !
    IF (mpime == ionode_id) THEN
      ! read 'selecq.fmt' file
      OPEN(UNIT = iunselecq, FILE = 'selecq.fmt', STATUS = 'old', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_frequencies', 'error opening selecq.fmt', 1)
      READ(iunselecq, *) totq
      ALLOCATE(selecq(totq), STAT = ierr)
      IF (ierr /= 0) CALL errore('read_frequencies', 'Error allocating selecq', 1)
      selecq(:) = 0
      READ(iunselecq, *) nqtotf
      IF (nqtotf /= nqf1 * nqf2 * nqf3) &
        CALL errore('read_frequencies', 'selecq.fmt is not calculated on the nqf1, nqf2, nqf3 mesh', 1)
      READ(iunselecq, *) selecq(:)
      CLOSE(iunselecq)
      !
      ! read frequencies from file
      dirname = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat'
      filfreq = TRIM(dirname) // '/' // 'freq'
      !OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', &
           FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_frequencies', 'error opening file ' // filfreq, iufilfreq)
      !READ(iufilfreq, '(5i7)') nqtotf, nqf1_, nqf2_, nqf3_, nmodes
      READ(iufilfreq) nqtotf, nqf1_, nqf2_, nqf3_, nmodes
      IF (nqtotf /= nqf1 * nqf2 * nqf3) &
        CALL errore('read_frequencies', 'e-ph mat elements were not calculated on the nqf1, nqf2, nqf3 mesh', 1)
    !
    ENDIF
    CALL mp_bcast(totq, ionode_id, inter_pool_comm)
    IF (mpime /= ionode_id) ALLOCATE(selecq(totq))
    CALL mp_bcast(selecq, ionode_id, inter_pool_comm)
    CALL mp_bcast(nqtotf, ionode_id, inter_pool_comm)
    CALL mp_bcast(nmodes, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    ALLOCATE(wf(nmodes, nqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_frequencies', 'Error allocating wf', 1)
    ALLOCATE(wqf(nqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_frequencies', 'Error allocating wqf', 1)
    ALLOCATE(xqf(3, nqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_frequencies', 'Error allocating xqf', 1)
    wf(:, :) = zero
    wqf(:) = 1.d0 / DBLE(nqtotf)
    xqf(:, :) = zero
    !
    IF (mpime == ionode_id) THEN
      DO iqq = 1, totq ! loop over q-points in fsthick
        iq = selecq(iqq)
        !
        !READ(iufilfreq, '(3f15.9)') xqf(:, iq)
        READ(iufilfreq) xqf(:, iq)
        DO imode = 1, nmodes
          !READ(iufilfreq, '(ES20.10)') wf(imode, iq)
          READ(iufilfreq) wf(imode, iq)
        ENDDO
      ENDDO
      CLOSE(iufilfreq)
      ! go from Ryd to eV
      wf(:, :) = wf(:, :) * ryd2ev ! in eV
      wsphmax = 1.1d0 * MAXVAL(wf(:, :)) ! increase by 10%
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    CALL mp_bcast(wf, ionode_id, inter_pool_comm)
    CALL mp_bcast(xqf, ionode_id, inter_pool_comm)
    CALL mp_bcast(wsphmax, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    ! create phonon grid
    !
    !dwsph = wsphmax / DBLE(nqstep - 1)
    dwsph = wsphmax / DBLE(nqstep)
    ALLOCATE(wsph(nqstep), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_frequencies', 'Error allocating wsph', 1)
    wsph(:) = 0.d0
    DO iwph = 1, nqstep
      !wsph(iwph) = DBLE(iwph - 1) * dwsph
      wsph(iwph) = DBLE(iwph) * dwsph
    ENDDO
    !
    ! HM: selecq will be used in read_kqmap and read_ephmat so we do not deallocate selecq here.
    !
    WRITE(stdout,'(/5x,a/)') 'Finish reading freq file'
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE read_frequencies
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_eigenvalues()
    !-----------------------------------------------------------------------
    !!
    !! read the eigenvalues obtained from a previous epw run
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_files,      ONLY : prefix, tmp_dir
    USE pwcom,         ONLY : ef
    USE epwcom,        ONLY : nkf1, nkf2, nkf3, degaussw, fsthick, mp_mesh_k
    USE eliashbergcom, ONLY : nkfs, nbndfs, dosef, ef0, ekfs, wkfs, xkfs, w0g, &
                              nkfs_all, ekfs_all, wkfs_all, xkfs_all, &
                              ixkf, ixkf_inv, &
                              ibnd_kfs_to_kfs_all, ibnd_kfs_all_to_kfs, nbndfs_all
    USE constants_epw, ONLY : ryd2ev, zero, eps8
    USE io_var,        ONLY : iufilegnv
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256) :: filegnv
    !! file name
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory where ikmap/egnv/freq/ephmat files are saved
    !
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER ::  nkftot
    !! Number of k-points
    INTEGER :: nkf1_, nkf2_, nkf3_
    !! Temporary variable for number of k-points along each direction
    INTEGER :: n, nbnd_
    !! Band indexes
    INTEGER :: ikdum
    !! Dummy
    INTEGER :: ikfs
    !! Counter on k-points for ekfs
    INTEGER :: ibndfs
    !! Counter on bands for ekfs
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: ebnd
    !! Local variable for energy
    REAL(KIND = DP) :: ebndmax1
    !! temporarily store the maximum value of ABS(ekfs-ekfs_all)
    REAL(KIND = DP) :: ebndmax2
    !! temporarily store the maximum value of ABS(ekfs-ekfs_all)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    REAL(KIND = DP), ALLOCATABLE :: ekf_(:, :)
    !! Temporary eigenvalues on the k point grid
    !
    IF (mpime == ionode_id) THEN
      !
      ! SP: Needs to be initialized
      nbnd_ = 0
      nkfs = 0
      !
      ! read eigenvalues on the irreducible fine k-mesh
      !
      dirname = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat'
      filegnv = TRIM(dirname) // '/' // 'egnv'
      !OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', &
           FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_eigenvalues', 'error opening file '//filegnv, iufilegnv)
      !
      !READ(iufilegnv, '(5i7)') nkftot, nkf1_, nkf2_, nkf3_, nkfs
      !READ(iufilegnv, '(i7, 5ES20.10)') nbnd_, ef, ef0, dosef, degaussw, fsthick
      READ(iufilegnv) nkftot, nkf1_, nkf2_, nkf3_, nkfs
      READ(iufilegnv) nbnd_, ef, ef0, dosef, degaussw, fsthick
      IF (nkf1 /= nkf1_ .OR. nkf2 /= nkf2_ .OR. nkf3 /= nkf3) &
        CALL errore('read_eigenvalues', 'e-ph mat elements were not calculated on the nkf1, nkf2, nkf3 mesh', 1)
      degaussw = degaussw * ryd2ev
      ef0 = ef0 * ryd2ev
      ef = ef * ryd2ev
      fsthick = fsthick * ryd2ev
      dosef = dosef / ryd2ev
      WRITE(stdout, '(5x, a32, ES20.10)') 'Fermi level (eV) = ', ef0
      WRITE(stdout, '(5x, a32, ES20.10)') 'DOS(states/spin/eV/Unit Cell) = ', dosef
      WRITE(stdout, '(5x, a32, ES20.10)') 'Electron smearing (eV) = ', degaussw
      WRITE(stdout,'(5x,a32,ES20.10)') 'Fermi window (eV) = ', fsthick
      IF (mp_mesh_k) THEN
        WRITE(stdout, '(5x, a, i9, a, i9/)') 'Nr irreducible k-points within the Fermi shell = ', nkfs, ' out of ', nkftot
      ELSE
        WRITE(stdout, '(5x, a, i9, a, i9/)') 'Nr k-points within the Fermi shell = ', nkfs, ' out of ', nkftot
      ENDIF
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    CALL mp_bcast(nkf1, ionode_id, inter_pool_comm)
    CALL mp_bcast(nkf2, ionode_id, inter_pool_comm)
    CALL mp_bcast(nkf3, ionode_id, inter_pool_comm)
    CALL mp_bcast(nkfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(degaussw, ionode_id, inter_pool_comm)
    CALL mp_bcast(ef0, ionode_id, inter_pool_comm)
    CALL mp_bcast(dosef, ionode_id, inter_pool_comm)
    CALL mp_bcast(fsthick, ionode_id, inter_pool_comm)
    CALL mp_bcast(ef, ionode_id, inter_pool_comm)
    !
    ALLOCATE(wkfs(nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error allocating wkfs', 1)
    ALLOCATE(xkfs(3, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error allocating xkfs', 1)
    wkfs(:) = zero
    xkfs(:, :) = zero
    !
    IF (mpime == ionode_id) THEN
      ! HM: nbndfs_all should be same with nbndfst.
      nbndfs_all = nbnd_
      nkfs_all = nkftot
    ENDIF
    CALL mp_bcast(nbndfs_all, ionode_id, inter_pool_comm)
    CALL mp_bcast(nkfs_all, ionode_id, inter_pool_comm)
    ALLOCATE(ekfs_all(nbndfs_all, nkfs_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error allocating ekfs_all', 1)
    ALLOCATE(wkfs_all(nkfs_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error allocating wkfs_all', 1)
    ALLOCATE(xkfs_all(3, nkfs_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error allocating xkfs_all', 1)
    xkfs_all(:, :) = zero
    wkfs_all(:) = zero
    ekfs_all(:, :) = zero
    ALLOCATE(ixkf(nkfs_all), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error allocating ixkf', 1)
    ALLOCATE(ixkf_inv(nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error allocating ixkf_inv', 1)
    ixkf(:) = 0
    ixkf_inv(:) = 0
    !
    IF (mpime == ionode_id) THEN
      !
      ! nbndfs - nr of bands within the Fermi shell
      !
      nbndfs = 0
      DO ik = 1, nkfs_all ! loop over irreducible k-points
        READ(iufilegnv) wkfs_all(ik), xkfs_all(:, ik), ikdum, ixkf(ik)
        IF (ik /= ikdum) CALL errore('read_eigenvalues', 'error reading file '//filegnv, ik)
        IF (ixkf(ik) /= 0) THEN
          ixkf_inv(ixkf(ik)) = ik
        ENDIF
        n = 0
        DO ibnd = 1, nbndfs_all
          READ(iufilegnv) ekfs_all(ibnd, ik)
          ! go from Ryd to eV
          ekfs_all(ibnd, ik) = ekfs_all(ibnd, ik) * ryd2ev
          IF (ABS(ekfs_all(ibnd, ik) - ef0) < fsthick) THEN
            n = n + 1
            IF (nbndfs < n) nbndfs = n
          ENDIF
        ENDDO
      ENDDO
      !
      WRITE(stdout, '(5x, i7, a/)') nbndfs, ' bands within the Fermi window'
      CLOSE(iufilegnv)
      !
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    CALL mp_bcast(nbndfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(ekfs_all, ionode_id, inter_pool_comm)
    CALL mp_bcast(wkfs_all, ionode_id, inter_pool_comm)
    CALL mp_bcast(xkfs_all, ionode_id, inter_pool_comm)
    CALL mp_bcast(ixkf_inv, ionode_id, inter_pool_comm)
    CALL mp_bcast(ixkf, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    ALLOCATE(ekfs(nbndfs, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error allocating ekfs', 1)
    ALLOCATE(w0g(nbndfs, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error allocating w0g', 1)
    ! sanity choice
    ekfs(:, :) = ef0 - 10.d0 * fsthick
    w0g(:, :) = zero
    !
    ALLOCATE(ibnd_kfs_all_to_kfs(nbndfs_all, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error allocating ibnd_kfs_all_to_kfs', 1)
    ALLOCATE(ibnd_kfs_to_kfs_all(nbndfs, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error allocating ibnd_kfs_to_kfs_all', 1)
    ibnd_kfs_all_to_kfs(:, :) = 0
    ibnd_kfs_to_kfs_all(:, :) = 0
    !
    IF (mpime == ionode_id) THEN
      DO ik = 1, nkfs_all
        IF (ixkf(ik) == 0) CYCLE
        ikfs = ixkf(ik)
        wkfs(ikfs) = wkfs_all(ik)
        xkfs(:, ikfs) = xkfs_all(:, ik)
        n = 0
        DO ibnd = 1, nbndfs_all
          IF (ABS(ekfs_all(ibnd, ik) - ef0) < fsthick) THEN
            n = n + 1
            ekfs(n, ikfs) = ekfs_all(ibnd, ik)
            w0g(n, ikfs)  = w0gauss((ekfs(n, ikfs) - ef0) / degaussw, 0) / degaussw
            ibnd_kfs_all_to_kfs(ibnd, ikfs) = n
            ibnd_kfs_to_kfs_all(n, ikfs) = ibnd
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    CALL mp_bcast(ekfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(wkfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(xkfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(w0g, ionode_id, inter_pool_comm)
    CALL mp_bcast(ibnd_kfs_all_to_kfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(ibnd_kfs_to_kfs_all, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout,'(/5x,a/)') 'Finish reading egnv file '
    !
    ebndmax1 = 0.0d0
    ebndmax2 = 0.0d0
    DO ik = 1, nkfs_all
      IF (ixkf(ik) == 0) CYCLE
      ikfs = ixkf(ik)
      DO ibnd = 1, nbndfs_all
        IF (ABS(ekfs_all(ibnd, ik) - ef0) < fsthick) THEN
          ibndfs = ibnd_kfs_all_to_kfs(ibnd, ikfs)
          ebnd = ekfs(ibndfs, ikfs) - ekfs_all(ibnd, ik)
          ebnd = ABS(ebnd)
          ebndmax1 = MAX(ebnd, ebndmax1)
          !!!!!!!! FOR DEBUG !!!!!!!!
          !WRITE(stdout, '(5x, f15.7, 1x, f15.7)') &
          !ekfs(ibndfs, ikfs), ekfs_all(ibnd, ik)
          !!!!!!!! FOR DEBUG !!!!!!!!
        ENDIF
      ENDDO
    ENDDO
    !
    DO ikfs = 1, nkfs
      ik = ixkf_inv(ikfs)
      DO ibndfs = 1, nbndfs
        IF (ABS(ekfs(ibndfs, ikfs) - ef0) < fsthick) THEN
          ibnd = ibnd_kfs_to_kfs_all(ibndfs, ikfs)
          ebnd = ekfs(ibndfs, ikfs) - ekfs_all(ibnd, ik)
          ebnd = ABS(ebnd)
          ebndmax2 = MAX(ebnd, ebndmax2)
          !!!!!!!! FOR DEBUG !!!!!!!!
          !WRITE(stdout, '(5x, f15.7, 1x, f15.7)') &
          !ekfs(ibndfs, ikfs), ekfs_all(ibnd, ik)
          !!!!!!!! FOR DEBUG !!!!!!!!
        ENDIF
      ENDDO
    ENDDO
    !
    IF ((ebndmax1 > eps8) .OR. (ebndmax2 > eps8)) THEN
      CALL errore('read_eigenvalues', 'Error: ekfs_all is not equal to ekfs.', 1)
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE read_eigenvalues
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_kqmap()
    !-----------------------------------------------------------------------
    !
    ! read the map index of k+(sign)q on the k-mesh
    ! this subroutine should be called after calling read_frequencies
    !
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout, ionode_id
    USE io_var,    ONLY : iufilikmap, iunselecq
    USE io_files,  ONLY : prefix, tmp_dir
    USE modes,     ONLY : nmodes
    USE epwcom,    ONLY : nkf1, nkf2, nkf3, nqstep
    USE elph2,     ONLY : nqtotf, xqf, totq, selecq, bztoibz
    USE grid,     ONLY : kpmq_map
    USE eliashbergcom, ONLY : ixkf, ixkff, xkff, xkfs, nkfs, ixkqf, ixqfs, &
                              nbndfs, nqfs, memlt_pool, ixkf
    USE constants_epw, ONLY : zero
    USE mp_global, ONLY : inter_pool_comm, npool
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE division,  ONLY : fkbounds
    USE low_lvl,   ONLY : mem_size_eliashberg
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256) :: filikmap
    !! Name of the file
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory where ikmap/egnv/freq/ephmat files are saved
    !
    INTEGER :: i, j, k, ik, nk, n
    !! Counter on k points
    INTEGER :: iq
    !! Counter on q points
    INTEGER :: nkq
    !! Index of k+sign*q on the fine k-mesh
    INTEGER :: nkftot
    !! Total number of k points
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k parallelization
    INTEGER :: nks
    !! Number of non-equivalent k points
    INTEGER :: ns
    !! Counter on rotation operations
    INTEGER :: iqq
    !! Q-point index from selecq.fmt window
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    INTEGER(8) :: imelt
    !! Memory allocated
    INTEGER, ALLOCATABLE :: index_(:, :)
    !! Index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell
    !
    REAL(KIND = DP) :: xk(3)
    !! coordinates of k points
    REAL(KIND = DP) :: xq(3)
    !! coordinates of q points
    !
    ALLOCATE(memlt_pool(npool), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_kqmap', 'Error allocating memlt_pool', 1)
    memlt_pool(:) = zero
    !
    ! get the size of arrays for frequency and eigenvalue variables allocated in
    ! read_frequencies and read_eigenvalues
    imelt = (nmodes + 4) * nqtotf + nqstep + (4 + 2 * nbndfs) * nkfs
    CALL mem_size_eliashberg(2, imelt)
    !
    nkftot = nkf1 * nkf2 * nkf3
    !
    ALLOCATE(ixkff(nkftot), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_kqmap', 'Error allocating ixkff', 1)
    ALLOCATE(bztoibz(nkftot), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_kqmap', 'Error allocating bztoibz', 1)
    ixkff(:) = 0
    bztoibz(:) = 0
    !
    IF (mpime == ionode_id) THEN
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! HM: no need to read selecq.fmt again since we already read selecq in read_frequencies
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! read 'selecq.fmt' file
      !OPEN(UNIT = iunselecq, FILE = 'selecq.fmt', STATUS = 'old', IOSTAT = ios)
      !READ(iunselecq, *) totq
      !ALLOCATE(selecq(totq), STAT = ierr)
      !IF (ierr /= 0) CALL errore('read_kqmap', 'Error allocating selecq', 1)
      !selecq(:) = 0
      !READ(iunselecq, *) nqtot
      !READ(iunselecq, *) selecq(:)
      !CLOSE(iunselecq)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      dirname = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat'
      filikmap = TRIM(dirname) // '/' // 'ikmap'
      !OPEN(UNIT = iufilikmap, FILE = filikmap, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufilikmap, FILE = filikmap, STATUS = 'unknown', &
           FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_kqmap', 'error opening file ' // filikmap, iufilikmap)
      !
      !READ(iufilikmap, *) ixkff(1:nkftot)
      READ(iufilikmap) ixkff(1:nkftot)
      READ(iufilikmap) bztoibz(1:nkftot)
      !
      CLOSE(iufilikmap)
    ENDIF
    !
    CALL mp_bcast(totq, ionode_id, inter_pool_comm)
    ! HM: no need to read selecq.fmt again since we already read selecq in read_frequencies
    !IF (mpime /= ionode_id) ALLOCATE(selecq(totq))
    CALL mp_bcast(ixkff, ionode_id, inter_pool_comm)
    CALL mp_bcast(bztoibz, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ! get the size of required memory for ixkqf, nqfs, index_
    imelt = nkfs + 2 * (upper_bnd - lower_bnd + 1) * nqtotf
    CALL mem_size_eliashberg(1, imelt)
    !
    ALLOCATE(ixkqf(lower_bnd:upper_bnd, nqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_kqmap', 'Error allocating ixkqf', 1)
    ALLOCATE(nqfs(nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_kqmap', 'Error allocating nqfs', 1)
    ALLOCATE(index_(lower_bnd:upper_bnd, nqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_kqmap', 'Error allocating index_', 1)
    ixkqf(:, :) = 0
    nqfs(:) = 0
    index_(:, :) = 0
    !
    ! nkfs - total nr. of k-points within the Fermi shell (fine mesh)
    !      - these are irreducible k-points if mp_mesh_k=.TRUE.
    ! nqtotf - total nr of q-points on the fine mesh
    !
    DO ik = lower_bnd, upper_bnd
      DO iqq = 1, totq
        iq = selecq(iqq)
        xk(:) = xkfs(:, ik)
        xq(:) = xqf(:, iq)
        !
        !  nkq - index of k+sign*q on the full fine k-mesh.
        !
        CALL kpmq_map(xk, xq, +1, nkq)
        !
        !  ixkqf(ik,iq) - index of k+sign*q on the fine k-mesh within the Fermi shell
        !
        ixkqf(ik, iq) = ixkff(nkq)
        !
        ! nqfs(ik) - nr of q-points at each k-point for which k+sign*q is within the Fermi shell
        ! index_   - index q-point on the full q-mesh for which k+sign*q is within the Fermi shell
        !
        IF (ixkqf(ik, iq) > 0) THEN
          nqfs(ik) = nqfs(ik) + 1
          index_(ik, nqfs(ik)) = iq
        ENDIF
      ENDDO
    ENDDO
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum(nqfs,  inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    ! get the size of required memory for ixqfs
    imelt = (upper_bnd - lower_bnd + 1) * MAXVAL(nqfs(:))
    CALL mem_size_eliashberg(1, imelt)
    !
    ALLOCATE(ixqfs(lower_bnd:upper_bnd, MAXVAL(nqfs(:))), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_kqmap', 'Error allocating ixqfs', 1)
    ixqfs(:, :) = 0
    !
    DO ik = lower_bnd, upper_bnd
      DO iq = 1, nqfs(ik)
        !
        ! ixqfs - index q-point on the full q-mesh for which k+sign*q is within the Fermi shell
        !
        ixqfs(ik, iq) = index_(ik, iq)
      ENDDO
    ENDDO
    !
    CALL mp_barrier(inter_pool_comm)
    !
    DEALLOCATE(index_, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_kqmap', 'Error deallocating index_', 1)
    !
    ! remove memory allocated for index_
    imelt = nqtotf * (upper_bnd - lower_bnd + 1)
    CALL mem_size_eliashberg(1, -imelt)
    !
    ! HM: selecq will be used in read_ephmat so we do not deallocate selecq here.
    !
    DO ik = 1, nkftot
      !!!!!!!! FOR DEBUG !!!!!!!!
      !WRITE(stdout, '(5x, I5, 1x, I5)') &
      !ixkff(ik), ixkf(bztoibz(ik))
      !!!!!!!! FOR DEBUG !!!!!!!!
      IF (ixkff(ik) /= ixkf(bztoibz(ik))) THEN
        CALL errore('read_kqmap', 'Error: ixkff(ik) is not equal to ixkf(bztoibz(ik)).', ik)
      ENDIF
    ENDDO
    !
    WRITE(stdout, '(/5x, a, i9/)') 'Max nr of q-points = ', MAXVAL(nqfs(:))
    WRITE(stdout, '(/5x, a/)') 'Finish reading ikmap files'
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE read_kqmap
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_ephmat()
    !-----------------------------------------------------------------------
    !!
    !! Read the electron-phonon matrix elements
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_var,        ONLY : iufilegnv, iufileph
    USE io_files,      ONLY : prefix, tmp_dir
    USE modes,         ONLY : nmodes
    USE elph2,         ONLY : nqtotf, wf, xqf, totq, selecq, bztoibz
    USE epwcom,        ONLY : eps_acustic, fsthick
    USE eliashbergcom, ONLY : nkfs, nbndfs, ef0, ekfs, g2, nqfs, &
                              xkfs, ixkf
    USE constants_epw, ONLY : ryd2ev, zero, eps8
    USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id
    USE division,      ONLY : fkbounds
    USE grid,          ONLY : kpmq_map
    USE low_lvl,       ONLY : set_ndnmbr, mem_size_eliashberg
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256) :: filephmat
    !! File name
    CHARACTER(LEN = 256) :: filegnv
    !! File name
    CHARACTER(LEN = 4) :: filelab
    !! File name
    !
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ikfs
    !! Counter on k-points within the fermi window
    INTEGER :: iq
    !! Counter on q-points
    INTEGER :: iqq
    !! Counter on q-points
    INTEGER :: nkq
    !! Index of k+sign*q on the fine k-mesh
    INTEGER :: ikqfs
    !! Index of k+sign*q on the irreducible k-mesh within the wiondow
    INTEGER :: ibnd, jbnd
    !! Counter on bands
    INTEGER :: imode
    !! Counter on modes
    INTEGER :: ipool
    !! Counter on pools
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k parallelization
    INTEGER :: tmp_pool_id
    !! Pool index read from file
    INTEGER(8) :: imelt
    !! Memory allocated
    INTEGER :: nmin, nmax
    !! Lower/upper bound index for .ephmat file read in current pool
    INTEGER :: nks
    !! Counter on k points within the Fermi shell
    INTEGER :: nnk
    !! Number of k-points within the Fermi shell
    INTEGER :: nnq(nkfs)
    !! Number of k+q points within the Fermi shell for a given k-point
    INTEGER :: npool_old
    !! npool which is read from the file output by write_ephmat
    INTEGER, ALLOCATABLE :: nkpool(:)
    !! nkpool(ipool) - sum of nr. of k points from pool 1 to pool ipool
    !
    REAL(KIND = DP) :: gmat
    !! Electron-phonon matrix element square
    REAL(KIND = DP) :: xk(3)
    !! coordinates of k points
    REAL(KIND = DP) :: xq(3)
    !! coordinates of q points
    !
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory where ephmat files are present
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ! get the size of the e-ph matrices that need to be stored in each pool
    imelt = (upper_bnd - lower_bnd + 1) * MAXVAL(nqfs(:)) * nbndfs**2 * nmodes
    CALL mem_size_eliashberg(2, imelt)
    !
    ALLOCATE(g2(lower_bnd:upper_bnd, MAXVAL(nqfs(:)), nbndfs, nbndfs, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_ephmat', 'Error allocating g2', 1)
    g2(:, :, :, :, :) = zero
    !
    ! go from Ryd to eV
    ! eps_acustic is given in units of cm-1 in the input file and converted to Ryd in epw_readin
    eps_acustic = eps_acustic * ryd2ev
    !
    WRITE(stdout, '(/5x, a/)') 'Start reading .ephmat files'
    !
    dirname = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat'
    !
    filegnv = TRIM(dirname) // '/' // 'npool'
    OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
    IF (ios /= 0) CALL errore('read_ephmat', 'error opening file ' // filegnv, iufilegnv)
    READ(iufilegnv) npool_old
    CLOSE(iufilegnv)
    !
    ALLOCATE(nkpool(npool_old), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_ephmat', 'Error allocating nkpool', 1)
    nkpool(:) = 0
    !
    DO ipool = 1, npool_old ! nr of pools
      CALL set_ndnmbr(0, ipool, 1, npool_old, filelab)
#if defined(__MPI)
      filephmat = TRIM(dirname) // '/' // 'ephmat' // filelab
#else
      filephmat = TRIM(dirname) // '/' // 'ephmat'
#endif
      !OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', &
           FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_ephmat', 'error opening file ' // filephmat, iufileph)
      !READ(iufileph, '(2i7)') tmp_pool_id, nkpool(ipool)
      READ(iufileph) tmp_pool_id, nkpool(ipool)
      IF (ipool /= tmp_pool_id)  CALL errore('read_ephmat', &
          'The header of .ephmat files is wrong.', 1)
      IF (ipool > 1) &
        nkpool(ipool) = nkpool(ipool) + nkpool(ipool - 1)
      !WRITE(stdout, '(2i7)') tmp_pool_id, nkpool(ipool)
      CLOSE(iufileph)
    ENDDO
    CALL mp_barrier(inter_pool_comm)
    !
    ! since the nkfs k-points within the Fermi shell are not evenly distrubed
    ! among the .ephmat files, we re-distribute them here among the npool-pools
    nmin = npool_old
    nmax = npool_old
    DO ipool = npool_old, 1, -1
      IF (lower_bnd <= nkpool(ipool)) THEN
        nmin = ipool
      ENDIF
      IF (upper_bnd <= nkpool(ipool)) THEN
        nmax = ipool
      ENDIF
    ENDDO
    !
    nnk = 0
    nnq(:) = 0
    DO ipool = 1, npool_old ! nr of pools
      CALL set_ndnmbr(0, ipool, 1, npool_old, filelab)
#if defined(__MPI)
      filephmat = TRIM(dirname) // '/' // 'ephmat' // filelab
#else
      filephmat = TRIM(dirname) // '/' // 'ephmat'
#endif
      !OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', &
           FORM ='unformatted', ACCESS = 'stream', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_ephmat', 'error opening file ' // filephmat, iufileph)
      !READ(iufileph, '(2i7)') tmp_pool_id, nks
      READ(iufileph) tmp_pool_id, nks
      IF (ipool >= nmin .AND. ipool <= nmax) THEN
        DO iqq = 1, totq ! loop over q-points
          iq = selecq(iqq)
          xq(:) = xqf(:, iq)
          DO ik = 1, nks ! loop over k-points in the pool
            ikfs = ik + nnk
            xk(:) = xkfs(:, ikfs)
            !
            !  nkq - index of k+sign*q on the full fine k-mesh.
            !
            CALL kpmq_map(xk, xq, +1, nkq)
            ikqfs = ixkf(bztoibz(nkq))
            IF (ikqfs > 0) THEN
              IF ((MINVAL(ABS(ekfs(:, ikfs) - ef0)) < fsthick) .AND. &
                  (MINVAL(ABS(ekfs(:, ikqfs) - ef0)) < fsthick)) THEN
                nnq(ikfs) = nnq(ikfs) + 1
                DO imode = 1, nmodes ! loop over phonon modes
                  DO ibnd = 1, nbndfs ! loop over iband's
                    IF (ABS(ekfs(ibnd, ikfs) - ef0) < fsthick) THEN
                      DO jbnd = 1, nbndfs ! loop over jband's
                        IF (ABS(ekfs(jbnd, ikqfs) - ef0) < fsthick) THEN
                          READ(iufileph) gmat
                          IF (ikfs >= lower_bnd .AND. ikfs <= upper_bnd) THEN
                            ! go from Ryd to eV
                            IF (wf(imode, iq) > eps_acustic) THEN
                              g2(ikfs, nnq(ikfs), ibnd, jbnd, imode) = gmat * ryd2ev * ryd2ev
                            ELSE
                              g2(ikfs, nnq(ikfs), ibnd, jbnd, imode) = zero
                            ENDIF
                          ENDIF
                        ENDIF ! ekq
                      ENDDO ! jbnd
                    ENDIF ! ekk
                  ENDDO ! ibnd
                ENDDO ! imode
              ENDIF ! ekk and ekq
            ENDIF ! ikqfs > 0
          ENDDO ! ik
        ENDDO ! iqq
        CLOSE(iufileph)
      ENDIF ! ipool
      nnk = nnk + nks
      IF (ipool == npool_old .AND. nnk /= nkfs)  CALL errore('read_ephmat', &
          'nnk should be equal to nkfs', 1)
    ENDDO ! ipool
    !
    CALL mp_barrier(inter_pool_comm)
    DEALLOCATE(selecq, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_ephmat', 'Error allocating selecq', 1)
    !
    WRITE(stdout, '(/5x, a/)') 'Finish reading .ephmat files'
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE read_ephmat
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE file_open_ephmat(lrepmatw2_restart)
    !----------------------------------------------------------------------------
    !!
    !! This routine opens all the files needed to save el-ph matrices for SC
    !! Adopted from io_transport.f90/iter_open
    !! 11/2022: Hari Paudyal
    !
    USE kinds,         ONLY : DP
    USE io_files,      ONLY : prefix, tmp_dir, delete_if_present
    USE io_var,        ONLY : iufileph, iufilegnv
    USE mp_global,     ONLY : my_pool_id, npool
    USE elph2,         ONLY : lrepmatw2_merge
    USE io_global,     ONLY : ionode_id
    USE low_lvl,       ONLY : set_ndnmbr
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: lrepmatw2_restart(npool)
    !! To restart opening files
    !
    ! Local variables
    !
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory to hold files
    CHARACTER(LEN = 256) :: filephmat
    !! Name e-ph mat file
    CHARACTER(LEN = 256) :: filegnv
    !! File name
    CHARACTER(LEN = 4) :: filelab
    !!
    LOGICAL :: exst
    !! Logical for existence of files
    LOGICAL :: exst2
    !! Logical for existence of files
    !!
    INTEGER :: dummy_int
    !! Dummy INTEGER for reading
    INTEGER :: ind
    !! Temp. index
    INTEGER :: npool_old
    !! npool which is read from the file output by write_ephmat
    INTEGER :: ios
    !! IO error message
    REAL(KIND = DP) :: dummy_real
    !! Dummy variable for reading
    !
    !
    dirname = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat'
    !
#if defined(__MPI)
    CALL set_ndnmbr(0, my_pool_id + 1, 1, npool, filelab)
    filephmat = TRIM(dirname) // '/' // 'ephmat' // filelab
#else
    filephmat = TRIM(dirname) // '/' // 'ephmat'
#endif
    !
    INQUIRE(FILE = 'restart.fmt', EXIST = exst)
    INQUIRE(FILE = filephmat, EXIST = exst2)
    !
    ! The restart.fmt exist - we try to restart
    IF (exst) THEN
      !
      IF (exst2) THEN
        !OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', &
        !            POSITION = 'append', FORM = 'formatted', IOSTAT = ios)
        OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', &
                    POSITION = 'rewind', FORM = 'unformatted', ACCESS = 'stream', ACTION = 'readwrite', IOSTAT = ios)
        IF (ios /= 0) CALL errore('file_open_ephmat', 'error opening file ' // filephmat, iufileph)
        ! This is done to move the pointer to the right position after a restart
        IF (lrepmatw2_restart(my_pool_id + 1) > 0) THEN
          READ(iufileph) dummy_int, dummy_int  ! my_pool_id, nks
          DO ind = 1, lrepmatw2_restart(my_pool_id + 1)
            READ(iufileph) dummy_real
          ENDDO
        ENDIF
        !
        ! HM: check the npool file
        filegnv = TRIM(dirname) // '/' // 'npool'
        OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
        IF (ios /= 0) CALL errore('read_ephmat', 'error opening file ' // filegnv, iufilegnv)
        READ(iufilegnv) npool_old
        CLOSE(iufilegnv)
        IF (npool /= npool_old) CALL errore('file_open_ephmat','Number of cores is different',1)
        !
      ELSE
        CALL errore('file_open_ephmat', 'A restart.fmt is present but not the prefix.ephmat folder', 1)
      ENDIF
      !
      lrepmatw2_merge = lrepmatw2_restart(my_pool_id + 1)
      !
    ELSE ! no restart file present
      !
      IF (exst2) THEN
        ! The file should not exist, we remove it
        CALL delete_if_present(filephmat, .TRUE.)
      ENDIF
      !
      !OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', &
           FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
      IF (ios /= 0) CALL errore('file_open_ephmat', 'error opening file ' // filephmat, iufileph)
      !
      lrepmatw2_merge = 0
      !
    ENDIF ! restart
    !
    lrepmatw2_restart(:) = 0
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE file_open_ephmat
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE write_ephmat(iqq, iq, lrepmatw2_restart, first_cycle)
    !-----------------------------------------------------------------------
    !!
    !!  This routine writes the elph matrix elements in a format required
    !!  by Eliashberg equations
    !!
    !!  Use matrix elements, electronic eigenvalues and phonon frequencies
    !!  from ep-wannier interpolation
    !!
    !-----------------------------------------------------------------------
    USE kinds,      ONLY : DP
    USE io_global,  ONLY : stdout, ionode_id
    USE io_var,     ONLY : iufilfreq, iufilegnv, iufileph, iunrestart
    USE io_files,   ONLY : prefix, tmp_dir
    USE modes,      ONLY : nmodes
    USE epwcom,     ONLY : nbndsub, fsthick, ngaussw, degaussw, shortrange, &
                           nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, efermi_read, &
                           fermi_energy, mp_mesh_k
    USE pwcom,      ONLY : ef
    USE elph2,      ONLY : etf, ibndmin, ibndmax, nkqf, epf17, wkf, nkf, &
                           nqtotf, wf, xqf, nkqtotf, efnew, nbndfst, nktotf, &
                           lrepmatw2_merge, totq, bztoibz
    USE eliashbergcom, ONLY : nkfs, ekfs, wkfs, xkfs, dosef, ixkf, nbndfs, &
                              ekfs_all, wkfs_all, xkfs_all
    USE constants_epw, ONLY : ryd2ev, ryd2mev, two, eps8, eps6, zero
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id, npool
    USE division,      ONLY : fkbounds
    USE grid,          ONLY : kpmq_map
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Q-point index from selecq.fmt window
    INTEGER, INTENT(in) :: iq
    !! Q-point index from full grid
    INTEGER, INTENT(inout) :: lrepmatw2_restart(npool)
    !! Current position inside the file during writing
    LOGICAL, INTENT(in) :: first_cycle
    !! Check wheter this is the first cycle after a restart.
    !
    ! Local variables
    !
    CHARACTER(LEN = 256) :: filfreq
    !! Name freq file
    CHARACTER(LEN = 256) :: filegnv
    !! Name eigenvalue file
    CHARACTER(LEN = 256) :: filephmat
    !! Name e-ph mat file
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory to save ikmap/egnv/freq/ephmat files
    !
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! k+q-point index
    INTEGER :: ikfs
    !! Counter on the k-point index for ekfs
    INTEGER :: nkq
    !! Index of k+sign*q on the fine k-mesh
    INTEGER :: ikqfs
    !! Index of k+sign*q on the irreducible k-mesh within the wiondow
    INTEGER :: ibnd, jbnd, pbnd
    !! Counter on bands
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: ipool
    !! Counter on npool
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: nkftot
    !! Total number of k+q points
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k or q parallelization
    INTEGER :: nks
    !! Number of k-point on the current pool
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    INTEGER :: ifil
    !! Temporary running index
    INTEGER :: n
    !! Integer for the degenerate average over eigenstates
    INTEGER :: ind(npool)
    !! Temporary index
    INTEGER :: dummy
    !! Dummy variable for writing
    REAL(KIND = DP) :: xk(3)
    !! coordinates of k points
    REAL(KIND = DP) :: xq(3)
    !! coordinates of q points
    REAL(KIND = DP) :: epf2_deg(nbndfst, nbndfst, nmodes)
    !! Epc in degeneracies
    REAL(KIND = DP) :: w_1, w_2
    !! Temporary electronic energy
    REAL(KIND = DP) :: tmp_g2(nbndfst * nbndfst * nmodes * nkf)
    !! Temporary index
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: wq
    !! phonon freq on the fine grid
    REAL(KIND = DP) :: inv_wq
    !! inverse phonon greq
    REAL(KIND = DP):: g2
    !! Electron-phonon matrix element square
    REAL(KIND = DP) :: dummy_real
    !! Dummy variable for reading
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! Function to compute the density of states at the Fermi level
    REAL(KIND = DP), EXTERNAL :: efermig
    !! Return the fermi energy
    REAL(KIND = DP) :: dummy_real3(3)
    !! Dummy variable for reading
    !
    ind(:)    = 0
    tmp_g2(:) = 0
    dummy     = 0
    !
    dirname = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat'
    !
    ! write phonon frequencies to file
    IF (my_pool_id == 0) THEN
      filfreq = TRIM(dirname) // '/' // 'freq'
      IF (iqq == 1) THEN
        !OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
        OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', &
             FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
        IF (ios /= 0) CALL errore('write_ephmat', 'error opening file ' // filfreq, iufilfreq)
      ELSEIF ((iqq .NE. 1) .AND. first_cycle) THEN
        !OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', POSITION = 'append', FORM = 'formatted', IOSTAT = ios)
        !OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', POSITION = 'append', FORM = 'unformatted', IOSTAT = ios)
        OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', &
             POSITION = 'rewind', FORM = 'unformatted', ACCESS = 'stream', ACTION = 'readwrite', IOSTAT = ios)
        IF (ios /= 0) CALL errore('write_ephmat', 'error opening file ' // filfreq, iufilfreq)
        !
        ! HM: iqq - 1 should be the same with the index read from restart.fmt
        !WRITE(stdout, '(5x,a,i8)')'iqq - 1 =', (iqq - 1)
        ! This is done to move the pointer to the right position after a restart
        READ(iufilfreq) dummy, dummy, dummy, dummy, dummy
        DO ifil = 1, iqq - 1
          READ(iufilfreq) dummy_real3(:)
          DO imode = 1, nmodes
            !WRITE(iufilfreq, '(ES20.10)') wf(imode, iq)
            READ(iufilfreq) dummy_real
          ENDDO
        ENDDO
      ENDIF
      !IF (iqq == 1) WRITE(iufilfreq, '(5i7)') nqtotf, nqf1, nqf2, nqf3, nmodes
      IF (iqq == 1) WRITE(iufilfreq) nqtotf, nqf1, nqf2, nqf3, nmodes
      !WRITE(iufilfreq, '(3f15.9)') xqf(:, iq)
      WRITE(iufilfreq) xqf(:, iq)
      DO imode = 1, nmodes
        !WRITE(iufilfreq, '(ES20.10)') wf(imode, iq)
        WRITE(iufilfreq) wf(imode, iq)
      ENDDO
      !
      FLUSH(iufilfreq)
      !
    ENDIF
    !
    ! Fermi level and corresponding DOS
    !
    ! since wkf(:,ikq) = 0 these bands do not bring any contribution to ef0 or dosef
    !
    IF (efermi_read) THEN
      ef0 = fermi_energy
    ELSE
      ef0 = efnew
      !ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk)
      ! if some bands are skipped (nbndskip /= 0), nelec has already been recalculated
      ! in ephwann_shuffle
    ENDIF
    !
    dosef = dos_ef(ngaussw, degaussw, ef0, etf, wkf, nkqf, nbndsub)
    ! N(Ef) in the equation for lambda is the DOS per spin
    dosef = dosef / two
    !
    ! find the bounds of k-dependent arrays in the parallel case
    nkftot = nktotf
    CALL fkbounds(nkftot, lower_bnd, upper_bnd)
    !
    IF (iqq == 1) THEN
      !
      ! find fermicount - nr of k-points within the Fermi shell per pool
      ! for mp_mesh_k=true. femicount is the nr of irreducible k-points within the Fermi shell per pool
      !
      fermicount = 0
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        !
        IF (MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) THEN
          fermicount = fermicount + 1
        ENDIF
        !
      ENDDO
      !
      ! nks = irr nr of k-points within the Fermi shell (fine mesh)
      nks = fermicount
      !
      ! collect contributions from all pools (sum over k-points)
      CALL mp_sum(nks, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
      IF (my_pool_id == 0) THEN
        !
        ! write eigenvalues to file
        filegnv = TRIM(dirname) // '/' // 'egnv'
        !OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
        OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', &
             FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
        IF (ios /= 0) CALL errore('write_ephmat', 'error opening file ' // filegnv, iufilegnv)
        IF (nks /= nkfs) CALL errore('write_ephmat', &
          'nks should be equal to nr. of irreducible k-points within the Fermi shell on the fine mesh', 1)
        !WRITE(iufilegnv, '(5i7)') nkftot, nkf1, nkf2, nkf3, nks
        !WRITE(iufilegnv, '(i7,5ES20.10)') nbndfst, ef, ef0, dosef, degaussw,fsthick
        WRITE(iufilegnv) nkftot, nkf1, nkf2, nkf3, nks
        WRITE(iufilegnv) nbndfst, ef, ef0, dosef, degaussw, fsthick
        ikfs = 0
        DO ik = 1, nkftot
          !WRITE(iufilegnv, '(4f15.9)') wkfs(ik), xkfs(:, ik)
          IF (MINVAL(ABS(ekfs_all(:, ik) - ef)) < fsthick) THEN
            ikfs = ikfs + 1
            WRITE(iufilegnv) wkfs_all(ik), xkfs_all(:, ik), ik, ikfs
          ELSE
            WRITE(iufilegnv) wkfs_all(ik), xkfs_all(:, ik), ik, 0
          ENDIF
          DO ibnd = 1, nbndfst
            !WRITE(iufilegnv, '(ES20.10)') ekfs(ibnd, ik)
            WRITE(iufilegnv) ekfs_all(ibnd, ik)
          ENDDO
        ENDDO
        CLOSE(iufilegnv)
        !
        ! write npool to file
        filegnv = TRIM(dirname) // '/' // 'npool'
        OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', &
             FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
        IF (ios /= 0) CALL errore('write_ephmat', 'error opening file ' // filegnv, iufilegnv)
        WRITE(iufilegnv) npool
        CLOSE(iufilegnv)
      ENDIF
      !
    ENDIF ! iqq
    !
    ! write the e-ph matrix elements in the Bloch representation on the fine mesh
    ! in .ephmat files (one for each pool)
    !
    !IF (iqq == 1) WRITE(iufileph, '(2i7)') my_pool_id+1, fermicount
    IF (iqq == 1) WRITE(iufileph) my_pool_id + 1, fermicount
    !
    xq(:) = xqf(:, iq)
    !
    ! nkf - nr of k-points in the pool (fine mesh)
    ! for mp_mesh_k = true nkf is nr of irreducible k-points in the pool
    !
    DO ik = 1, nkf
      !
      ikk = 2 * ik - 1
      ikq = ikk + 1
      !
      ! Average the el-ph matrix elements on degenerate bands and phonon modes.
      ! HP: Adapted from: io_transport.f90
      !
      ! Average over the k electrons
      DO imode = 1, nmodes
        DO jbnd = 1, nbndfst
          DO ibnd = 1, nbndfst
            w_1 = etf(ibndmin - 1 + ibnd, ikk)
            g2  = zero
            n   = 0
            DO pbnd = 1, nbndfst
              w_2 = etf(ibndmin - 1 + pbnd, ikk)
              IF (ABS(w_2-w_1) < eps6) THEN
                n = n + 1
                g2 = g2 + ABS(epf17(jbnd, pbnd, imode, ik))**two
              ENDIF
            ENDDO
            epf2_deg(jbnd, ibnd, imode) = DSQRT(g2 / FLOAT(n))
          ENDDO
        ENDDO
      ENDDO
      epf17(:, :, :, ik) = epf2_deg(:, :, :)
      !
      ! Average over the k+q electrons
      DO imode = 1, nmodes
        DO jbnd = 1, nbndfst
          DO ibnd = 1, nbndfst
            w_1 = etf(ibndmin - 1 + jbnd, ikq)
            g2 = 0.d0
            n  = 0
            DO pbnd = 1, nbndfst
              w_2 = etf(ibndmin - 1 + pbnd, ikq)
              IF (ABS(w_2 - w_1) < eps6) THEN
                n = n + 1
                g2 = g2 + ABS(epf17(pbnd, ibnd, imode, ik))**two
              ENDIF
            ENDDO
            epf2_deg(jbnd, ibnd, imode) = g2 / FLOAT(n)
          ENDDO
        ENDDO
      ENDDO
      !
      ! Note that we already took the square above
      epf17(:, :, :, ik) = epf2_deg(:, :, :)
      !
      ! go only over irreducible k-points
      !
      ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
      !
      ikfs = ixkf(lower_bnd + ik - 1)
      IF (ikfs > 0) THEN
        xk(:) = xkfs(:, ikfs)
        !
        !  nkq - index of k+sign*q on the full fine k-mesh.
        !
        CALL kpmq_map(xk, xq, +1, nkq)
        !
        IF (mp_mesh_k) THEN
          !
          ikqfs = ixkf(bztoibz(nkq))
          !
        ELSE
          ! full k-point grid
          ikqfs = ixkf(nkq)
          !
        ENDIF
        !
        IF (ikqfs > 0) THEN
          DO imode = 1, nmodes ! phonon modes
            wq = wf(imode, iq)
            inv_wq =  1.0d0 / (two * wq)
            !
            DO ibnd = 1, nbndfst
              IF (ABS(ekfs(ibnd, ikfs) - ef0) < fsthick) THEN
                DO jbnd = 1, nbndfst
                  IF (ABS(ekfs(jbnd, ikqfs) - ef0) < fsthick) THEN
                    !
                    ! here we take into account the zero-point DSQRT(hbar/2M\omega)
                    ! with hbar = 1 and M already contained in the eigenmodes
                    ! g2 is Ry^2, wkf must already account for the spin factor
                    !
                    IF (shortrange .AND. (ABS(xqf(1, iq)) > eps8 .OR. ABS(xqf(2, iq)) > eps8 &
                         .OR. ABS(xqf(3, iq)) > eps8 )) THEN
                      ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary
                      !     number, in which case its square will be a negative number.
                      g2 = REAL(epf17(jbnd, ibnd, imode, ik) * inv_wq)
                    ELSE
                      g2 = ABS(epf17(jbnd, ibnd, imode, ik)) * inv_wq
                    ENDIF
                    ind(my_pool_id + 1) = ind(my_pool_id + 1) + 1
                    tmp_g2(ind(my_pool_id + 1)) = g2
                  ENDIF
                ENDDO ! jbnd
              ENDIF
            ENDDO ! ibnd
          ENDDO ! imode
        ENDIF ! ikqfs
      ENDIF ! ixkf
      !
    ENDDO ! ik's
    !
    lrepmatw2_merge = lrepmatw2_merge + ind(my_pool_id + 1)
    !
    ! Save to file restart information in formatted way for possible restart
    lrepmatw2_restart(:) = 0
    lrepmatw2_restart(my_pool_id + 1) = lrepmatw2_merge
    CALL mp_sum(lrepmatw2_restart, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (ind(my_pool_id + 1) > 0) THEN
      !WRITE(*, *) my_pool_id + 1, ind(my_pool_id + 1), lrepmatw2_merge
      DO ifil = 1, ind(my_pool_id + 1)
        !WRITE(iufileph, '(ES20.10)') tmp_g2(ifil)
        WRITE(iufileph) tmp_g2(ifil)
      ENDDO
      !
      FLUSH(iufileph)
      !
    ENDIF
    !
    IF (my_pool_id == 0) THEN
      dummy = 0
      ! format is compatible with IBTE
      OPEN(UNIT = iunrestart, FILE = 'restart.fmt')
      WRITE(iunrestart, *) iqq
      WRITE(iunrestart, *) dummy
      WRITE(iunrestart, *) dummy
      WRITE(iunrestart, *) npool
      DO ipool = 1, npool
        WRITE(iunrestart, *) lrepmatw2_restart(ipool)
      ENDDO
      DO ipool = 1, npool
       WRITE(iunrestart, *) dummy
      ENDDO
      CLOSE(iunrestart)
    ENDIF
    !
    IF (iqq == totq) THEN
      IF (my_pool_id == 0) CLOSE(iufilfreq)
      CLOSE(iufileph)
      DEALLOCATE(ekfs, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_ephmat', 'Error deallocating ekfs', 1)
      DEALLOCATE(wkfs, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_ephmat', 'Error deallocating wkfs', 1)
      DEALLOCATE(xkfs, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_ephmat', 'Error deallocating xkfs', 1)
      DEALLOCATE(ekfs_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_ephmat', 'Error deallocating ekfs_all', 1)
      DEALLOCATE(wkfs_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_ephmat', 'Error deallocating wkfs_all', 1)
      DEALLOCATE(xkfs_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_ephmat', 'Error deallocating xkfs_all', 1)
      DEALLOCATE(ixkf, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_ephmat', 'Error deallocating ixkf', 1)
      !
      WRITE(stdout, '(5x, a32, d24.15)') 'Fermi level (eV) = ', ef0 * ryd2ev
      WRITE(stdout, '(5x, a32, d24.15)') 'DOS(states/spin/eV/Unit Cell) = ', dosef / ryd2ev
      WRITE(stdout, '(5x, a32, d24.15)') 'Electron smearing (eV) = ', degaussw * ryd2ev
      WRITE(stdout, '(5x, a32, d24.15)') 'Fermi window (eV) = ', fsthick * ryd2ev
      WRITE(stdout, '(5x, a)')           ' '
      WRITE(stdout, '(5x, a/)')           'Finish writing .ephmat files'
      !
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE write_ephmat
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE count_kpoints
    !-----------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE epwcom,    ONLY : nbndsub, fsthick, ngaussw, degaussw, &
                          efermi_read, fermi_energy, mp_mesh_k
    USE pwcom,     ONLY : nelec, ef
    USE klist_epw, ONLY : isk_dummy
    USE elph2,     ONLY : etf, nkqf, wkf, nkf, nkqtotf, nktotf
    USE constants_epw, ONLY : two
    USE mp,        ONLY : mp_barrier, mp_sum
    USE mp_global, ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    ! Local variables
    !
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: nks
    !! Number of k-point on the current pool
    !
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: dosef
    !! density of states at the Fermi level
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! Function to compute the density of states at the Fermi level
    REAL(KIND = DP), EXTERNAL :: efermig
    !! Return the fermi energy
    !
    ! Fermi level and corresponding DOS
    !
    ! since wkf(:,ikq) = 0 these bands do not bring any contribution to ef0 or dosef
    !
    IF (efermi_read) THEN
      ef0 = fermi_energy
    ELSE
      ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk_dummy)
    ENDIF
    !
    dosef = dos_ef(ngaussw, degaussw, ef0, etf, wkf, nkqf, nbndsub)
    ! N(Ef) in the equation for lambda is the DOS per spin
    dosef = dosef / two
    !
    ! fermicount = nr of k-points within the Fermi shell per pool
    !
    fermicount = 0
    DO ik = 1, nkf
      !
      ikk = 2 * ik - 1
      !
      IF (MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) &
        fermicount = fermicount + 1
      !
    ENDDO
    !
    ! nks =  nr of k-points within the Fermi shell (fine mesh)
    nks = fermicount
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum(nks, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (mp_mesh_k) THEN
      WRITE(stdout, '(/5x, a, i9, a, i9/)') 'Nr irreducible k-points within the Fermi shell = ', nks, ' out of ', nktotf
    ELSE
      WRITE(stdout, '(/5x, a, i9, a, i9/)') 'Nr k-points within the Fermi shell = ', nks, ' out of ', nktotf
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE count_kpoints
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kmesh_fine()
    !-----------------------------------------------------------------------
    !!
    !!   This routine defines the nr. of k-points on the fine k-mesh
    !!   within the Fermi shell
    !!
    USE kinds,     ONLY : DP
    USE io_global, ONLY : ionode_id, stdout
    USE io_files,  ONLY : prefix, tmp_dir
    USE epwcom,    ONLY : nkf1, nkf2, nkf3, fsthick, mp_mesh_k
    USE pwcom,     ONLY : ef
    USE elph2,     ONLY : xkf, wkf, etf, nkf, nkqtotf, ibndmin, ibndmax, nktotf, nbndfst
    USE eliashbergcom, ONLY : nkfs, ixkf, xkfs, wkfs, ekfs, nbndfs, xkfs_all,  &
                              wkfs_all, ekfs_all
    USE constants_epw, ONLY : zero
    USE mp_global, ONLY : inter_pool_comm, npool
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    USE division,  ONLY : fkbounds
    USE kfold,     ONLY : backtoBZ
    !
    IMPLICIT NONE
    !
    INTEGER :: nk
    !! Counter on k points
    INTEGER :: nks
    !! Counter on k points within the Fermi shell
    INTEGER :: ikk
    !! k-point index
    INTEGER :: nkf_mesh
    !! Total number of k points
    !! These are irreducible k-points if mp_mesh_k = .TRUE.
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k parallelization
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: xx, yy, zz
    !! Temporary variables
    REAL(KIND = DP), ALLOCATABLE :: xkf_(:, :)
    !! Temporary k point grid
    REAL(KIND = DP), ALLOCATABLE :: wkf_(:)
    !! Temporary weights on the k point grid
    REAL(KIND = DP), ALLOCATABLE :: ekf_(:, :)
    !! Temporary eigenvalues on the k point grid
    !
    nkf_mesh = nktotf
    nbndfs = nbndfst
    !
    ALLOCATE(ekf_(nbndfs, nkf_mesh), STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error allocating ekf_', 1)
    ALLOCATE(wkf_(nkf_mesh), STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error allocating wkf_', 1)
    ALLOCATE(xkf_(3, nkf_mesh), STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error allocating xkf_', 1)
    ALLOCATE(ixkf(nkf_mesh), STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error allocating ixkf', 1)
    xkf_(:, :) = zero
    ekf_(:, :) = zero
    wkf_(:) = zero
    ixkf(:) = 0
    !
    CALL fkbounds(nkf_mesh, lower_bnd, upper_bnd)
    !
    ! nkf - nr of k-blocks in the pool (fine grid)
    !
    DO nk = 1, nkf
      ikk = 2 * nk - 1
      xkf_(:, lower_bnd + nk - 1) = xkf(:, ikk)
      wkf_(lower_bnd + nk - 1) = wkf(ikk)
      ekf_(:, lower_bnd + nk - 1) = etf(ibndmin:ibndmax, ikk)
    ENDDO
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum(ekf_, inter_pool_comm)
    CALL mp_sum(xkf_, inter_pool_comm)
    CALL mp_sum(wkf_, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (mpime == ionode_id) THEN
      !
      IF (mp_mesh_k) THEN
        WRITE(stdout, '(/5x, a, i9/)') 'Nr. of irreducible k-points on the uniform grid: ', nkf_mesh
      ELSE
        WRITE(stdout, '(/5x, a, i9/)') 'Nr. of k-points on the uniform grid: ', nkf_mesh
      ENDIF
      !
      nkfs = 0
      DO nk = 1, nkf_mesh
        IF (MINVAL(ABS(ekf_(:, nk) - ef)) < fsthick) THEN
          nkfs = nkfs + 1
          ixkf(nk) = nkfs
        ELSE
          ixkf(nk) = 0
        ENDIF
        !  bring back into to the first BZ
        xx = xkf_(1, nk) * nkf1
        yy = xkf_(2, nk) * nkf2
        zz = xkf_(3, nk) * nkf3
        CALL backtoBZ(xx, yy, zz, nkf1, nkf2, nkf3)
        xkf_(1, nk) = xx / DBLE(nkf1)
        xkf_(2, nk) = yy / DBLE(nkf2)
        xkf_(3, nk) = zz / DBLE(nkf3)
      ENDDO
      !
    ENDIF
    CALL mp_bcast(nkfs, ionode_id, inter_pool_comm)
    !
    ALLOCATE(ekfs(nbndfs, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error allocating ekfs', 1)
    ALLOCATE(wkfs(nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error allocating wkfs', 1)
    ALLOCATE(xkfs(3, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error allocating xkfs', 1)
    xkfs(:, :) = zero
    wkfs(:) = zero
    ekfs(:, :) = zero
    !
    ALLOCATE(ekfs_all(nbndfs, nkf_mesh), STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error allocating ekfs_all', 1)
    ALLOCATE(wkfs_all(nkf_mesh), STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error allocating wkfs_all', 1)
    ALLOCATE(xkfs_all(3, nkf_mesh), STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error allocating xkfs_all', 1)
    xkfs_all(:, :) = zero
    wkfs_all(:) = zero
    ekfs_all(:, :) = zero
    !
    IF (mpime == ionode_id) THEN
      nks = 0
      DO nk = 1, nkf_mesh
        IF (MINVAL(ABS(ekf_(:, nk) - ef)) < fsthick) THEN
          nks = nks + 1
          IF (nks > nkf_mesh) CALL errore('kmesh_fine', 'too many k-points', 1)
          wkfs(nks)    = wkf_(nk)
          xkfs(:, nks) = xkf_(:, nk)
          ekfs(:, nks) = ekf_(:, nk)
        ENDIF
      ENDDO
      wkfs_all(:)    = wkf_(:)
      xkfs_all(:, :) = xkf_(:, :)
      ekfs_all(:, :) = ekf_(:, :)
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    CALL mp_bcast(ixkf, ionode_id, inter_pool_comm)
    CALL mp_bcast(xkfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(wkfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(ekfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(xkfs_all, ionode_id, inter_pool_comm)
    CALL mp_bcast(wkfs_all, ionode_id, inter_pool_comm)
    CALL mp_bcast(ekfs_all, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    DEALLOCATE(ekf_, STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error deallocating ekf_', 1)
    DEALLOCATE(wkf_, STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error deallocating wkf_', 1)
    DEALLOCATE(xkf_, STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error deallocating xkf_', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE kmesh_fine
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kqmap_fine()
    !-----------------------------------------------------------------------
    !!
    !! this routine finds the index of k+sign*q on the fine k-mesh
    !!
    USE kinds,     ONLY : DP
    USE symm_base, ONLY : s, t_rev, time_reversal, set_sym_bl
    USE epwcom,    ONLY : nkf1, nkf2, nkf3, mp_mesh_k
    USE elph2,     ONLY : nqtotf, nktotf, xqf, map_rebal, bztoibz
    USE eliashbergcom, ONLY : ixkff, ixkf, xkfs, nkfs, nqfs
    USE constants_epw, ONLY : eps5, zero
    USE io_global, ONLY : stdout, ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    USE division,  ONLY : fkbounds
    USE grid,      ONLY : kpmq_map, kpoint_grid_epw
    USE io_files,  ONLY : prefix, tmp_dir, create_directory
    USE io_var,    ONLY : iufilikmap
    USE mp_world,  ONLY : mpime
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256) :: filikmap
    !! Name of the file
    !
    INTEGER :: i, j, k, ik, nk, n, ikbz
    !! Counter on k points
    INTEGER :: iq
    !! Counter on q points
    INTEGER :: nkq
    !! Index of k+sign*q on the fine k-mesh
    INTEGER :: nkftot
    !! Total number of k points
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k parallelization
    INTEGER :: ns
    !! Counter on rotation operations
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios
    !! IO error message
    INTEGER, ALLOCATABLE :: index_(:, :)
    !! Index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell
    INTEGER, ALLOCATABLE :: bztoibz_tmp(:)
    !! Dummy for bztoibz
    !
    REAL(KIND = DP) :: xk(3)
    !! coordinates of k points
    REAL(KIND = DP) :: xq(3)
    !! coordinates of q points
    !
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory to save ikmap/egnv/freq/ephmat files
    !
    nkftot = nkf1 * nkf2 * nkf3
    !
    ALLOCATE(ixkff(nkftot), STAT = ierr)
    IF (ierr /= 0) CALL errore('kqmap_fine', 'Error allocating ixkff', 1)
    ixkff(:) = 0
    !
    ! to map k+q onto k we need to define the index of k on the full mesh (ixkff)
    ! using index of the k-point within the Fermi shell (ixkf)
    !
    IF (mp_mesh_k) THEN
      ! SP - July 2020
      ! We should not recompute bztoibz
      DO ikbz = 1, nkftot
        ixkff(ikbz) = ixkf(bztoibz(ikbz))
      ENDDO
      !
    ELSE
      !
      ALLOCATE(bztoibz_tmp(nkftot), STAT = ierr)
      IF (ierr /= 0) CALL errore('kqmap_fine', 'Error allocating bztoibz_tmp', 1)
      ! full k-point grid
      DO nk = 1, nkftot
        ixkff(nk) = ixkf(nk)
        bztoibz_tmp(nk) = nk
      ENDDO
      !
    ENDIF
    !
    dirname = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat'
    CALL create_directory(TRIM(dirname))
    !
    IF (mpime == ionode_id) THEN
      filikmap = TRIM(dirname) // '/' // 'ikmap'
      !OPEN(UNIT = iufilikmap, FILE = filikmap, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufilikmap, FILE = filikmap, STATUS = 'unknown', &
           FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
      IF (ios /= 0) CALL errore('kqmap_fine', 'error opening file ' // filikmap, iufilikmap)
      !
      !WRITE(iufilikmap, *) ixkff(1:nkftot)
      WRITE(iufilikmap) ixkff(1:nkftot)
      IF (mp_mesh_k) THEN
        WRITE(iufilikmap) bztoibz(1:nkftot)
      ELSE
        WRITE(iufilikmap) bztoibz_tmp(1:nkftot)
      ENDIF
      !
      CLOSE(iufilikmap)
    ENDIF
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ALLOCATE(nqfs(nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('kqmap_fine', 'Error allocating nqfs', 1)
    ALLOCATE(index_(lower_bnd:upper_bnd, nqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('kqmap_fine', 'Error allocating index_', 1)
    nqfs(:) = 0
    index_(:, :) = 0
    !
    ! find the index of k+sign*q on the fine k-mesh
    ! nkfs - total nr. of k-points within the Fermi shell
    !      - these are irreducible k-points if mp_mesh_k = true
    ! nqtotf - total nr of q-points on the fine mesh
    !
    DO ik = lower_bnd, upper_bnd
      DO iq = 1, nqtotf
        xk(:) = xkfs(:, ik)
        xq(:) = xqf(:, iq)
        !
        ! find nkq - index of k+sign*q on the full fine k-mesh.
        !
        CALL kpmq_map(xk, xq, +1, nkq)
        !
        ! ixkff(nkq) - index of k+sign*q on the fine k-mesh within the Fermi shell
        ! nqfs(ik) - nr of q-points at each k-point for which k+sign*q is within the Fermi shell
        ! index_   - index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell
        !
        IF (ixkff(nkq) > 0) THEN
          nqfs(ik) = nqfs(ik) + 1
          index_(ik, nqfs(ik)) = iq
        ENDIF
      ENDDO ! loop over full set of q-points (fine mesh)
    ENDDO ! loop over k-points within the Fermi shell in each pool (fine mesh)
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum(nqfs, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    DEALLOCATE(ixkff, STAT = ierr)
    IF (ierr /= 0) CALL errore('kqmap_fine', 'Error deallocating ixkff', 1)
    DEALLOCATE(index_, STAT = ierr)
    IF (ierr /= 0) CALL errore('kqmap_fine', 'Error deallocating index_', 1)
    DEALLOCATE(nqfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('kqmap_fine', 'Error deallocating nqfs', 1)
    IF (.NOT. mp_mesh_k) THEN
      DEALLOCATE(bztoibz_tmp, STAT = ierr)
      IF (ierr /= 0) CALL errore('kqmap_fine', 'Error deallocating bztoibz_tmp', 1)
    ENDIF
    !
    IF (mp_mesh_k) THEN
      WRITE(stdout, '(/5x, a/)') 'Finish mapping k+sign*q onto the fine irreducibe &
                                  k-mesh and writing .ikmap file'
    ELSE
      WRITE(stdout, '(/5x, a/)') 'Finish mapping k+sign*q onto the fine k-mesh &
                                  and and writing .ikamp file'
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE kqmap_fine
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE check_restart_ephwrite(iq_restart)
    !-----------------------------------------------------------------------
    !!
    !!   This routine checks the variables in restart while writing ephmat
    !!   06/2020 : Hari Paudyal
    !!   11/2022 : Hari Paudyal (updated)
    !!
    USE kinds,         ONLY : DP
    USE io_files,      ONLY : prefix, tmp_dir
    USE epwcom,        ONLY : nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, fsthick, mp_mesh_k
    USE io_var,        ONLY : iufilfreq, iufilegnv
    USE modes,         ONLY : nmodes
    USE elph2,         ONLY : nqtotf
    USE constants_epw, ONLY : zero
    USE mp_global,     ONLY : my_pool_id, npool
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq_restart
    !! restart step of q-point
    !
    ! local variables
    CHARACTER(LEN = 256) :: filfreq
    !! file name
    CHARACTER(LEN = 256) :: filegnv
    !! file name
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory to save ikmap/egnv/freq/ephmat files
    !
    LOGICAL :: exst, exst2
    !! Logical for existence of files
    !
    INTEGER :: iq
    !! Counter on q points
    INTEGER :: imode
    !! Counter on modes
    INTEGER :: iqq
    !! Q-point index from selecq.fmt window
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    INTEGER ::  nkftot_
    !! Temporary variable for number of k-points
    INTEGER ::  nkfs_
    !! Temporary variable for number of irr k-points
    INTEGER :: nkf1_, nkf2_, nkf3_
    !! Temporary variable for number of k-points along each direction
    INTEGER ::  nqtotf_
    !! Temporary variable for number of q-points
    INTEGER ::  nmodes_
    !! Temporary variable for number of modes
    INTEGER :: nqf1_, nqf2_, nqf3_
    !! Temporary variable for number of q-points along each direction
    !
    REAL(KIND = DP), ALLOCATABLE :: wf_tmp(:, :)
    ! temp. interpolated eigenfrequencie
    REAL(KIND = DP), ALLOCATABLE :: xqf_tmp(:, :)
    ! temp. fine q point grid
    !
    dirname = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat'
    !
    INQUIRE(FILE = TRIM(dirname) // '/' // 'ikmap', EXIST = exst)
    !
    IF (.NOT. exst) THEN
      CALL errore('check_restart_ephwrite', 'A restart.fmt is present but the directory ' // TRIM(prefix) // '.ephmat' // &
                  ' is not found. Remove the restart.fmt file and restart.', 1)
    ENDIF
    !
    ! read header of egnv file
    filegnv = TRIM(dirname) // '/' // 'egnv'
    !OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
    OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', &
         FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
    IF (ios /= 0) CALL errore('check_restart_ephwrite', 'error opening file '//filegnv, iufilegnv)
    !
    !READ(iufilegnv, '(5i7)') nkftot_, nkf1_, nkf2_, nkf3_, nkfs_
    READ(iufilegnv) nkftot_, nkf1_, nkf2_, nkf3_, nkfs_
    IF (nkf1 /= nkf1_ .OR. nkf2 /= nkf2_ .OR. nkf3 /= nkf3_) &
      CALL errore('check_restart_ephwrite', 'nkf1, nkf2, nkf3 is not consistent with restart.fmt', 1)
    CLOSE(iufilegnv)
    !
    ! read freq file
    filfreq = TRIM(dirname) // '/' // 'freq'
    !OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
    OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', &
         FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
    IF (ios /= 0) CALL errore('check_restart_ephwrite', 'error opening file ' // filfreq, iufilfreq)
    !READ(iufilfreq, '(5i7)') nqtotf_, nqf1_, nqf2_, nqf3_, nmodes_
    READ(iufilfreq) nqtotf_, nqf1_, nqf2_, nqf3_, nmodes_
    IF (nqf1 /= nqf1_ .OR. nqf2 /= nqf2_ .OR. nqf3 /= nqf3_) &
      CALL errore('check_restart_ephwrite', 'nqf1, nqf2, nqf3 is not consistent with restart.fmt', 1)
    !
    ! read up to iq_restart and write fresh freq file to avoid overwritting
    ALLOCATE(wf_tmp(nmodes, iq_restart), STAT = ierr)
    IF (ierr /= 0) CALL errore('check_restart_ephwrite', 'Error allocating wf_tmp', 1)
    ALLOCATE(xqf_tmp(3, iq_restart), STAT = ierr)
    IF (ierr /= 0) CALL errore('check_restart_ephwrite', 'Error allocating xqf_tmp', 1)
    wf_tmp(:, :) = zero
    xqf_tmp(:, :) = zero
    !
    DO iqq = 1, iq_restart ! read only up to iq_restart
      !READ(iufilfreq, '(3f15.9)') xqf_tmp(:, iqq)
      READ(iufilfreq) xqf_tmp(:, iqq)
      DO imode = 1, nmodes
        !READ(iufilfreq, '(ES20.10)') wf_tmp(imode, iqq)
        READ(iufilfreq) wf_tmp(imode, iqq)
      ENDDO
    ENDDO
    CLOSE(iufilfreq)
    !
    ! Write fresh freq file
    IF (my_pool_id == 0) THEN
      !
      !OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', &
           FORM = 'unformatted', ACCESS = 'stream', IOSTAT = ios)
      IF (ios /= 0) CALL errore('check_restart_ephwrite', 'error opening file ' // filfreq, iufilfreq)
      !
      !WRITE(iufilfreq, '(5i7)') nqtotf, nqf1, nqf2, nqf3, nmodes
      WRITE(iufilfreq) nqtotf, nqf1, nqf2, nqf3, nmodes
      DO iqq = 1, iq_restart ! write only up to iq_restart
        !WRITE(iufilfreq, '(3f15.9)') xqf_tmp(:, iqq)
        WRITE(iufilfreq) xqf_tmp(:, iqq)
        DO imode = 1, nmodes
          !WRITE(iufilfreq, '(ES20.10)') wf_tmp(imode, iqq)
          WRITE(iufilfreq) wf_tmp(imode, iqq)
        ENDDO
      ENDDO
      CLOSE(iufilfreq)
    ENDIF
    !
    DEALLOCATE(wf_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('check_restart_ephwrite', 'Error deallocating wf_tmp', 1)
    DEALLOCATE(xqf_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('check_restart_ephwrite', 'Error deallocating xqf_tmp', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE check_restart_ephwrite
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gap_distribution_FS(itemp, cname)
    !-----------------------------------------------------------------------
    !
    ! This routine writes to files the distribution of the superconducting
    ! gap on the Fermi surface
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : fsthick
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : agap, nkfs, nbndfs, ef0, ekfs, w0g
    USE constants_epw, ONLY : kelvin2eV, zero, eps4, eps5
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    CHARACTER(LEN = 256), INTENT(in) :: cname
    !! character in output file name
    !
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! output file name
    !
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ibnd
    !! Counter on band
    INTEGER :: ibin
    !! Counter on bins
    INTEGER :: nbin
    !! Number of bins
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: temp
    !! Temperature in K
    REAL(KIND = DP) :: dbin
    !! Step size in nbin
    REAL(KIND = DP) :: delta_min, delta_max
    !! Min/Max value of superconducting gap
    REAL(KIND = DP) :: weight
    !! Variable for weight
    REAL(KIND = DP), ALLOCATABLE :: delta_k_bin(:)
    !! Histogram superconducting gap
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    temp = gtemp(itemp) / kelvin2eV
    !
    delta_min = MINVAL(agap(:, :))
    delta_max = MAXVAL(agap(:, :))
    !
    WRITE(stdout, '(5x, a, 2f12.6, a)') 'Min. / Max. values of superconducting gap = ', &
      delta_min * 1000.d0, delta_max * 1000.d0, ' meV'
    !
    IF (delta_min > zero) THEN
      delta_min = 0.9d0 * delta_min
    ELSE
      delta_min = 1.1d0 * delta_min
    ENDIF
    !
    IF (delta_max > zero) THEN
      delta_max = 1.1d0 * delta_max
    ELSE
      delta_max = 0.9d0 * delta_max
    ENDIF
    !
    !nbin = NINT((delta_max - delta_min) / eps4) + 1
    !dbin = (delta_max - delta_min) / DBLE(nbin)
    dbin = 3.0d-5 !eV
    nbin = NINT((delta_max - delta_min) / dbin) + 1
    !
    ALLOCATE(delta_k_bin(nbin), STAT = ierr)
    IF (ierr /= 0) CALL errore('gap_distribution_FS', 'Error allocating delta_k_bin', 1)
    delta_k_bin(:) = zero
    !
    DO ik = 1, nkfs
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          ibin = NINT((agap(ibnd, ik) - delta_min) / dbin) + 1
          !ibin = NINT(agap(ibnd, ik) / dbin) + 1
          weight = w0g(ibnd, ik)
          delta_k_bin(ibin) = delta_k_bin(ibin) + weight
        ENDIF
      ENDDO
    ENDDO
    !
    IF (temp < 10.d0) THEN
      WRITE(name1, '(a, a1, a4, a14, f4.2)') TRIM(prefix), '.', cname, '_aniso_gap0_00', temp
    ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
      WRITE(name1, '(a, a1, a4, a13, f5.2)') TRIM(prefix), '.', cname, '_aniso_gap0_0', temp
    ELSEIF (temp >= 100.d0) THEN
      WRITE(name1, '(a, a1, a4, a12, f6.2)') TRIM(prefix), '.', cname, '_aniso_gap0_', temp
    ENDIF
    !
    OPEN(UNIT = iufilgap, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
    IF (ios /= 0) CALL errore('gap_distribution_FS', 'error opening file ' // name1, iufilgap)
    WRITE(iufilgap, '(2a20)') '#     T [K]    ', '\rho(delta_nk) [meV]'
    DO ibin = 1, nbin
      IF (delta_k_bin(ibin) > eps5) &
        WRITE(iufilgap,'(2ES20.10)') temp + delta_k_bin(ibin) / MAXVAL(delta_k_bin(:)), &
                                     (dbin * DBLE(ibin) + delta_min) * 1000.d0
    ENDDO
    CLOSE(iufilgap)
    !
    DEALLOCATE(delta_k_bin, STAT = ierr)
    IF (ierr /= 0) CALL errore('gap_distribution_FS', 'Error deallocating delta_k_bin', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE gap_distribution_FS
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gap_FS(itemp)
    !-----------------------------------------------------------------------
    !
    ! This routine writes to files the superconducting gap on the Fermi surface
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgapFS
    USE io_files,      ONLY : prefix
    USE cell_base,     ONLY : bg
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : fsthick, nkf1, nkf2, nkf3
    USE elph2,         ONLY : gtemp, bztoibz
    USE eliashbergcom, ONLY : agap, nkfs, nbndfs, ef0, ekfs, ixkff, ekfs_all, &
                              nbndfs_all, ibnd_kfs_all_to_kfs, ixkf, nkfs_all
    USE constants_epw, ONLY : kelvin2eV, zero
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    !
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! output file name
    CHARACTER(LEN = 256) :: cname
    !! character in output file name
    !
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ikfs
    !! Counter on k-points: ikfs = ixkf(bztoibz(ik))
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: i, j, k
    !! Counter on grid points nkf1, nkf2, nkf3
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: temp
    !! Temperature in K
    REAL(KIND = DP) :: rdum
    !! Dummy for real numbers
    REAL(KIND = DP) :: x1, x2, x3
    !! Cartesian coordinates of grid points nkf1, nkf2, nkf3
    REAL(KIND = DP), ALLOCATABLE :: agap_tmp(:, :)
    !! Temporary array for superconducting gap at ik, ibnd
    !
    temp = gtemp(itemp) / kelvin2eV
    !
    cname = 'imag'
    !
    ! RM - If the k-point is outside the Fermi shell,
    ! ixkff(ik)=0 and agap_tmp(:,0) = 0.0
    !
    !ALLOCATE(agap_tmp(nbndfs, 0:nkfs), STAT = ierr)
    ALLOCATE(agap_tmp(0:nbndfs, 0:nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('gap_FS', 'Error allocating agap_tmp', 1)
    agap_tmp(1:nbndfs, 1:nkfs) = agap(1:nbndfs, 1:nkfs)
    agap_tmp(:, 0) = zero
    agap_tmp(0, :) = zero
    !
    ! SP & RM: .cube file for VESTA plotting (only if iverbosity = 2)
    !
    IF (iverbosity == 2) THEN
      !
      DO ibnd = 1, nbndfs
        !
        IF (ibnd < 10) THEN
          ! We make the assumption that there are no superconductor with Tc
          ! higher than 999 K.
          IF (temp < 10.d0) THEN
            WRITE(name1, '(a, a1, a4, a14, f4.2, a1, i1, a5)') TRIM(prefix), '.', cname, '_aniso_gap0_00', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 100.d0) THEN
            WRITE(name1, '(a, a1, a4, a13, f5.2, a1, i1, a5)') TRIM(prefix), '.', cname, '_aniso_gap0_0', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 1000.d0) THEN
            WRITE(name1, '(a, a1, a4, a12, f6.2, a1, i1, a5)') TRIM(prefix), '.', cname, '_aniso_gap0_', temp, '_', ibnd, '.cube'
          ENDIF
        ELSEIF (ibnd < 100) THEN
          IF (temp < 10.d0) THEN
            WRITE(name1, '(a, a1, a4, a14, f4.2, a1, i2, a5)') TRIM(prefix), '.', cname, '_aniso_gap0_00', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 100.d0 .AND. temp > 9.9999d0) THEN
            WRITE(name1, '(a, a1, a4, a13, f5.2, a1, i2, a5)') TRIM(prefix), '.', cname, '_aniso_gap0_0', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 1000.d0 .AND. temp > 99.9999d0) THEN
            WRITE(name1, '(a, a1, a4, a12, f6.2, a1, i2, a5)') TRIM(prefix), '.', cname, '_aniso_gap0_', temp, '_', ibnd, '.cube'
          ENDIF
        ELSEIF (ibnd < 1000) THEN
          IF (temp < 10.d0) THEN
            WRITE(name1, '(a, a1, a4, a14, f4.2, a1, i3, a5)') TRIM(prefix), '.', cname, '_aniso_gap0_00', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 100.d0 .AND. temp > 9.9999d0) THEN
            WRITE(name1, '(a, a1, a4, a13, f5.2, a1, i3, a5)') TRIM(prefix), '.', cname, '_aniso_gap0_0', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 1000.d0 .AND. temp > 99.9999d0) THEN
            WRITE(name1, '(a, a1, a4, a12, f6.2, a1, i3, a5)') TRIM(prefix), '.', cname, '_aniso_gap0_', temp, '_', ibnd, '.cube'
          ENDIF
        ELSE
          CALL errore('gap_FS', 'Too many bands', 1)
        ENDIF
        !
        OPEN(UNIT = iufilgapFS, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('gap_FS', 'error opening file ' // name1, iufilgapFS)
        WRITE(iufilgapFS, *) 'Cubfile created from EPW calculation'
        WRITE(iufilgapFS, *) 'gap'
        WRITE(iufilgapFS, '(i5, 3f12.6)') 1, 0.0d0, 0.0d0, 0.0d0
        WRITE(iufilgapFS, '(i5, 3f12.6)') nkf1, (bg(i, 1) / DBLE(nkf1), i = 1, 3)
        WRITE(iufilgapFS, '(i5, 3f12.6)') nkf2, (bg(i, 2) / DBLE(nkf2), i = 1, 3)
        WRITE(iufilgapFS, '(i5, 3f12.6)') nkf3, (bg(i, 3) / DBLE(nkf3), i = 1, 3)
        WRITE(iufilgapFS, '(i5, 4f12.6)') 1, 1.0d0, 0.0d0, 0.0d0, 0.0d0
        ! agap_tmp is written to file in meV
        WRITE(iufilgapFS, '(6f12.6)') (agap_tmp(ibnd, ixkff(ik)) * 1000.d0, ik = 1, nkf1 * nkf2 * nkf3)
        CLOSE(iufilgapFS)
      ENDDO
      !
      ! HP: Write in .frmsf format compatible with fermisurfer program
      IF (temp < 10.d0) THEN
        WRITE(name1, '(a, a1, a4, a14, f4.2, a6)') TRIM(prefix), '.', cname, '_aniso_gap0_00', temp, '.frmsf'
      ELSEIF (temp < 100.d0 .AND. temp > 9.9999d0) THEN
        WRITE(name1, '(a, a1, a4, a13, f5.2, a6)') TRIM(prefix), '.', cname, '_aniso_gap0_0', temp, '.frmsf'
      ELSEIF (temp < 1000.d0 .AND. temp > 99.9999d0) THEN
        WRITE(name1, '(a, a1, a4, a12, f6.2, a6)') TRIM(prefix), '.', cname, '_aniso_gap0_', temp, '.frmsf'
      ENDIF
      OPEN(UNIT = iufilgapFS, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('gap_FS', 'error opening file ' // name1, iufilgapFS)
      !
      WRITE(iufilgapFS, '(3i5)') nkf1, nkf2, nkf3
      WRITE(iufilgapFS, '(i5)') 1
      WRITE(iufilgapFS, '(i5)') nbndfs
      WRITE(iufilgapFS, '(3f12.6)') (bg(i, 1), i = 1, 3)
      WRITE(iufilgapFS, '(3f12.6)') (bg(i, 2), i = 1, 3)
      WRITE(iufilgapFS, '(3f12.6)') (bg(i, 3), i = 1, 3)
      ! HM: Outputting a dummy value in the .frmsf file may form fake Fermi surfaces.
      !     To avoid using a dummy value for the states outside of fsthick window, 
      !     use ekfs_all instead of ekfs.
      !WRITE(iufilgapFS, '(6f12.6)') ((ekfs(ibnd, ixkff(ik)) - ef0, ik = 1, nkf1 * nkf2 * nkf3), ibnd = 1, nbndfs)
      !WRITE(iufilgapFS, '(6f12.6)') ((agap_tmp(ibnd, ixkff(ik)) * 1000.d0, ik = 1, nkf1 * nkf2 * nkf3), ibnd = 1, nbndfs)
      i = 1
      DO ibnd = 1, nbndfs_all
        DO ik = 1, nkf1 * nkf2 * nkf3
          IF (i == 6) THEN
            WRITE(iufilgapFS, '(f12.6)') ekfs_all(ibnd, bztoibz(ik)) - ef0
          ELSE
            WRITE(iufilgapFS, '(f12.6)', advance='no') ekfs_all(ibnd, bztoibz(ik)) - ef0
          ENDIF
          i = i + 1
          IF (i == 7) i = 1
        ENDDO
      ENDDO
      i = 1
      DO ibnd = 1, nbndfs_all
        DO ik = 1, nkf1 * nkf2 * nkf3
          ikfs = ixkf(bztoibz(ik))
          IF (ikfs == 0) THEN
            rdum = 0.0d0
          ELSE
            rdum = agap_tmp(ibnd_kfs_all_to_kfs(ibnd, ikfs), ikfs) * 1000.d0
          ENDIF
          IF (i == 6) THEN
            WRITE(iufilgapFS, '(f12.6)') rdum
          ELSE
            WRITE(iufilgapFS, '(f12.6)', advance='no') rdum
          ENDIF
          i = i + 1
          IF (i == 7) i = 1
        ENDDO
      ENDDO
      CLOSE(iufilgapFS)
      !
    ENDIF
    !
    ! SP & RM : Write on file the superconducting gap close to the Fermi surface
    ! along with
    !     Cartesian coordinate, band index, energy distance from Fermi level and
    !     gap value.
    !
    IF (temp < 10.d0) THEN
      WRITE(name1, '(a, a1, a4, a16, f4.2)') TRIM(prefix), '.', cname, '_aniso_gap_FS_00', temp
    ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
      WRITE(name1, '(a, a1, a4, a15, f5.2)') TRIM(prefix), '.', cname, '_aniso_gap_FS_0', temp
    ELSEIF (temp >= 100.d0) THEN
      WRITE(name1, '(a, a1, a4, a14, f6.2)') TRIM(prefix), '.', cname, '_aniso_gap_FS_', temp
    ENDIF
    OPEN(UNIT = iufilgapFS, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
    IF (ios /= 0) CALL errore('gap_FS', 'error opening file ' // name1, iufilgapFS)
    WRITE(iufilgapFS, '(a78)') '#               k-point                  Band Enk-Ef [eV]        delta(0) [meV]'
    DO i = 1, nkf1
      DO j = 1, nkf2
        DO k = 1, nkf3
          ik = k + (j - 1) * nkf3 + (i - 1) * nkf2 * nkf3
          IF (ixkff(ik) > 0) THEN
            DO ibnd = 1, nbndfs
              ! RM: Everything is in eV here.
              ! SP: Here take a 0.2 eV interval around the FS.
              IF (ABS(ekfs(ibnd, ixkff(ik)) - ef0) < fsthick) THEN
              !IF (ABS(ekfs(ibnd, ixkff(ik)) - ef0 ) < 0.2) THEN
                 x1 = bg(1, 1) * (i - 1) /nkf1 + bg(1, 2) * (j - 1) / nkf2 + bg(1, 3) *(k - 1) / nkf3
                 x2 = bg(2, 1) * (i - 1) /nkf1 + bg(2, 2) * (j - 1) / nkf2 + bg(2, 3) *(k - 1) / nkf3
                 x3 = bg(3, 1) * (i - 1) /nkf1 + bg(3, 2) * (j - 1) / nkf2 + bg(3, 3) *(k - 1) / nkf3
                 WRITE(iufilgapFS,'(3f12.6, i8, f12.6, f24.15)') x1, x2, x3, ibnd, &
                       ekfs(ibnd, ixkff(ik)) - ef0, agap_tmp(ibnd, ixkff(ik)) * 1000.d0
              ENDIF
            ENDDO ! ibnd
          ENDIF
        ENDDO  ! k
      ENDDO ! j
    ENDDO ! i
    CLOSE(iufilgapFS)
    !
    DEALLOCATE(agap_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('gap_FS', 'Error deallocating agap_tmp', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE gap_FS
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    !
  END MODULE io_eliashberg

  !
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
    !----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nstemp, fsthick
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nsiw, gap0, gap, agap, wsi, nznormi, znormi, deltai, &
                              aznormi, naznormi, adeltai, adeltaip, nkfs, nbndfs, ef0, ekfs, &
                              dosef, wkfs, w0g
    USE constants_epw, ONLY : kelvin2eV, eps6, zero
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE supercond, ONLY : free_energy
    USE low_lvl,   ONLY : mem_size_eliashberg
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
    INTEGER :: imelt
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
    !
    ! get the size of required allocated memory
    imelt = (1 + nbndfs * nkfs) * nstemp + (3 + 4 * nbndfs * nkfs) * nsiw(itemp)
    CALL mem_size_eliashberg(2, imelt)
    !
    ALLOCATE(gap(nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating gap', 1)
    ALLOCATE(agap(nbndfs, nkfs, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating agap', 1)
    ALLOCATE(deltai(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating deltai', 1)
    ALLOCATE(znormi(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating znormi', 1)
    ALLOCATE(nznormi(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating nznormi', 1)
    ALLOCATE(adeltai(nbndfs, nkfs, nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating adeltai', 1)
    ALLOCATE(adeltaip(nbndfs, nkfs, nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating adeltaip', 1)
    ALLOCATE(aznormi(nbndfs, nkfs, nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating aznormi', 1)
    ALLOCATE(naznormi(nbndfs, nkfs, nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'Error allocating naznormi', 1)
    gap(:) = zero
    agap(:, :, :) = zero
    deltai(:) = zero
    znormi(:) = zero
    nznormi(:) = zero
    adeltai(:, :, :) = zero
    adeltaip(:, :, :) = zero
    aznormi(:, :, :) = zero
    naznormi(:, :, :) = zero
    !
    IF (mpime == ionode_id) THEN
      !
      temp = gtemp(itemp) / kelvin2eV
      ! anisotropic case
      IF (temp < 10.d0) THEN
        WRITE(name1, 101) TRIM(prefix), '.imag_aniso_00', temp
      ELSEIF (temp >= 10.d0) THEN
        WRITE(name1, 102) TRIM(prefix), '.imag_aniso_0', temp
      ELSEIF (temp >= 100.d0) THEN
        WRITE(name1, 103) TRIM(prefix), '.imag_aniso_', temp
      ENDIF
      !
      OPEN(UNIT = iufilgap, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('eliashberg_read_aniso_iaxis', 'error opening file ' // name1, iufilgap)
      READ(iufilgap, '(a)') word
      DO iw = 1, nsiw(itemp) ! loop over omega
        DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              READ(iufilgap, '(5ES20.10)') omega, eband, &
                   aznormi(ibnd, ik, iw), adeltai(ibnd, ik, iw), naznormi(ibnd, ik, iw)
              IF (iw == 1) &
                agap(ibnd, ik, itemp) = adeltai(ibnd, ik, 1)
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
              znormi(iw) = znormi(iw) + weight * aznormi(ibnd, ik, iw)
              deltai(iw) = deltai(iw) + weight * adeltai(ibnd, ik, iw)
              nznormi(iw) = nznormi(iw) + weight * naznormi(ibnd, ik, iw)
            ENDIF
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! iw
      gap(itemp) = deltai(1)
      gap0 = gap(itemp)
      !
      CALL gap_FS(itemp)
      !
      IF (iverbosity == 2) &
        CALL free_energy(itemp)
      !
    ENDIF
    CALL mp_bcast(gap0, ionode_id, inter_pool_comm)
    CALL mp_bcast(gap, ionode_id, inter_pool_comm)
    CALL mp_bcast(agap, ionode_id, inter_pool_comm)
    CALL mp_bcast(deltai, ionode_id, inter_pool_comm)
    CALL mp_bcast(znormi, ionode_id, inter_pool_comm)
    CALL mp_bcast(nznormi, ionode_id, inter_pool_comm)
    CALL mp_bcast(adeltai, ionode_id, inter_pool_comm)
    CALL mp_bcast(adeltaip, ionode_id, inter_pool_comm)
    CALL mp_bcast(aznormi, ionode_id, inter_pool_comm)
    CALL mp_bcast(naznormi, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    101 FORMAT(a, a14, f4.2)
    102 FORMAT(a, a13, f5.2)
    103 FORMAT(a, a12, f6.2)
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
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : fsthick, laniso, liso
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nsiw, agap, wsi, &
                              naznormi, aznormi, adeltai, nznormi, znormi, &
                              deltai, nkfs, nbndfs, ef0, ekfs
    USE constants_epw, ONLY : kelvin2eV
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
    INTEGER :: ios
    !! IO error message
    !
    REAL(KIND = DP) :: temp
    !! Temperature in K
    !
    temp = gtemp(itemp) / kelvin2eV
    !
    cname = 'imag'
    !
    IF (laniso) THEN
      !
      IF (temp < 10.d0) THEN
        WRITE(name1, 101) TRIM(prefix), '.', cname, '_aniso_00', temp
      ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
        WRITE(name1, 102) TRIM(prefix), '.', cname, '_aniso_0', temp
      ELSEIF (temp >= 100.d0) THEN
        WRITE(name1, 103) TRIM(prefix), '.', cname, '_aniso_', temp
      ENDIF
      OPEN(UNIT = iufilgap, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('eliashberg_write_iaxis', 'error opening file ' // name1, iufilgap)
      !
      WRITE(iufilgap,'(5a20)') '#        w [eV]', 'Enk-Ef [eV]', 'znorm(w)', 'delta(w) [eV]', 'nznorm(w)'
      DO iw = 1, nsiw(itemp) ! loop over omega
        DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              WRITE(iufilgap, '(5ES20.10)') wsi(iw), ekfs(ibnd, ik)- ef0, &
                    aznormi(ibnd, ik, iw), adeltai(ibnd, ik, iw), naznormi(ibnd, ik, iw)
              IF (iw == 1) &
                agap(ibnd, ik, itemp) = adeltai(ibnd, ik, iw)
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
    ENDIF
    !
    ! isotropic case
    ! SP: Only write isotropic for laniso if user really wants that
    IF ((laniso .AND. iverbosity == 2) .OR. liso) THEN
      IF (temp < 10.d0) THEN
        WRITE(name1, 104) TRIM(prefix), '.', cname, '_iso_00', temp
      ELSEIF (temp >= 10.d0 .AND. temp < 100.d0 ) THEN
        WRITE(name1, 105) TRIM(prefix), '.', cname, '_iso_0', temp
      ELSEIF (temp >= 100.d0) THEN
        WRITE(name1, 106) TRIM(prefix), '.', cname, '_iso_', temp
      ENDIF
      OPEN(UNIT = iufilgap, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('eliashberg_write_iaxis', 'error opening file ' // name1, iufilgap)
      !
      WRITE(iufilgap, '(4a20)') 'w [eV]', 'znorm(w)', 'delta(w) [eV]', 'nznorm(w)'
      DO iw = 1, nsiw(itemp) ! loop over omega
        WRITE(iufilgap, '(4ES20.10)') wsi(iw), znormi(iw), deltai(iw), nznormi(iw)
      ENDDO
      CLOSE(iufilgap)
    ENDIF
    !
    101 FORMAT(a, a1, a4, a9, f4.2)
    102 FORMAT(a, a1, a4, a8, f5.2)
    103 FORMAT(a, a1, a4, a7, f6.2)
    104 FORMAT(a, a1, a4, a7, f4.2)
    105 FORMAT(a, a1, a4, a6, f5.2)
    106 FORMAT(a, a1, a4, a5, f6.2)
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
    !
    !
    ! This routine writes to files results from the solutions of the Eliashberg
    ! equations on the real-axis
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nqstep, fsthick, laniso, liso
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nsw, ws, gap, agap, delta, znorm, adelta, aznorm, &
                              nkfs, nbndfs, ef0, ekfs
    USE constants_epw, ONLY : kelvin2eV
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter for temperature
    CHARACTER(len=256), INTENT(in) :: cname
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
    INTEGER :: ios
    !! IO error message
    !
    REAL(KIND = DP) :: temp
    !! Temperature in K
    REAL(KIND = DP) :: var1, var2, var3, var4
    !! Temporary working variables
    !
    temp = gtemp(itemp) / kelvin2eV
    !
    IF (laniso) THEN
      IF (iverbosity == 2) THEN
        IF (temp < 10.d0) THEN
          WRITE(name1, 101) TRIM(prefix), '.', cname, '_aniso_00', temp
        ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
          WRITE(name1, 102) TRIM(prefix), '.', cname, '_aniso_0', temp
        ELSEIF (temp >= 100.d0) THEN
          WRITE(name1, 103) TRIM(prefix), '.', cname, '_aniso_', temp
        ENDIF
        OPEN(UNIT = iufilgap, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('eliashberg_write_raxis', 'error opening file ' // name1, iufilgap)
        WRITE(iufilgap, '(6a20)') '#        w [eV]', 'Enk-Ef [eV]', 'Re[znorm(w)]', 'Im[znorm(w)]', &
                                  'Re[delta(w)] [eV]', 'Im[delta(w)] [eV]'
      ENDIF
      !
      DO ik = 1, nkfs
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
            lgap = .TRUE.
            ! DO iw = 1, nsw
            DO iw = 1, nsw - 1   ! FG: this change is to prevent segfault in ws(iw+1) and adelta(*,*,iw+1)
              var1 = REAL(adelta(ibnd, ik, iw))
              var2 = REAL(adelta(ibnd, ik, iw + 1))
              var3 = var1 - ws(iw)
              var4 = var2 - ws(iw + 1)
              IF (lgap .AND. iw < nqstep .AND. var1 > 0.d0 .AND. var2 > 0.d0 .AND. var3 * var4 < 0.d0) THEN
                agap(ibnd, ik, itemp) = (var3 * ws(iw + 1) - var4 * ws(iw)) / (var3 - var4)
                lgap = .FALSE.
              ENDIF
              IF (iverbosity == 2) THEN
                WRITE(iufilgap, '(6ES20.10)') ws(iw), ekfs(ibnd, ik) - ef0, &
                      REAL(aznorm(ibnd, ik, iw)), AIMAG(aznorm(ibnd, ik, iw)), &
                      REAL(adelta(ibnd, ik, iw)), AIMAG(adelta(ibnd, ik, iw))
              ENDIF
            ENDDO ! iw
            IF (lgap) &
              agap(ibnd,ik,itemp) = REAL(adelta(ibnd,ik,1))
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
      IF (iverbosity == 2) &
        CLOSE(iufilgap)
      !
      CALL gap_distribution_FS(itemp, cname)
      !
    ENDIF
    !
    ! isotropic case
    ! SP: Only write isotropic for laniso if user really wants that
    IF ((laniso .AND. iverbosity == 2) .OR. liso) THEN
      IF (temp < 10.d0) THEN
        WRITE(name1, 104) TRIM(prefix), '.', cname, '_iso_00', temp
      ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
        WRITE(name1, 105) TRIM(prefix), '.', cname, '_iso_0', temp
      ELSEIF (temp >= 100.d0) THEN
        WRITE(name1, 106) TRIM(prefix), '.', cname, '_iso_', temp
      ENDIF
      OPEN(UNIT = iufilgap, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('eliashberg_write_raxis', 'error opening file ' // name1, iufilgap)
      WRITE(iufilgap, '(5a20)') 'w [eV]', 'Re[znorm(w)]', 'Im[znorm(w)]', 'Re[delta(w)] [eV]', 'Im[delta(w)] [eV]'
      !
      lgap = .TRUE.
      ! DO iw = 1, nsw
      DO iw = 1, nsw - 1   ! this change is to prevent segfault in delta(iw+1) and ws(iw+1)
        var1 = REAL(delta(iw))
        var2 = REAL(delta(iw + 1))
        var3 = var1 - ws(iw)
        var4 = var2 - ws(iw + 1)
        IF (lgap .AND. iw < nqstep .AND. var1 > 0.d0 .AND. var2 > 0.d0 .AND. var3 * var4 < 0.d0) THEN
            gap(itemp) = (var3 * ws(iw + 1) - var4 * ws(iw)) / (var3 - var4)
          lgap = .FALSE.
        ENDIF
        WRITE(iufilgap, '(5ES20.10)') ws(iw), REAL(znorm(iw)), AIMAG(znorm(iw)), &
                                      REAL(delta(iw)), AIMAG(delta(iw))
      ENDDO ! iw
      CLOSE(iufilgap)
      IF (lgap) &
        gap(itemp) = REAL(delta(1))
    ENDIF
    !
    101 FORMAT(a, a1, a4, a9, f4.2)
    102 FORMAT(a, a1, a4, a8, f5.2)
    103 FORMAT(a, a1, a4, a7, f6.2)
    104 FORMAT(a, a1, a4, a7, f4.2)
    105 FORMAT(a, a1, a4, a6, f5.2)
    106 FORMAT(a, a1, a4, a5, f6.2)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE eliashberg_write_raxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_write_cont_raxis(itemp, cname)
    !-----------------------------------------------------------------------
    !
    !
    ! This routine writes to files results from the solutions of the Eliashberg
    ! equations on the real-axis
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nqstep, fsthick, laniso, liso
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nsw, ws, gap, agap, delta, znorm, adelta, aznorm, &
                              nkfs, nbndfs, ef0, ekfs
    USE constants_epw, ONLY : kelvin2eV
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter for temperature
    CHARACTER(len=256), INTENT(in) :: cname
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
    INTEGER :: ios
    !! IO error message
    !
    REAL(KIND = DP) :: temp
    !! Temperature in K
    REAL(KIND = DP) :: var1, var2, var3, var4
    !! Temporary working variables
    !
    temp = gtemp(itemp) / kelvin2eV
    !
    IF (laniso) THEN
      IF (iverbosity == 2) THEN
        IF (temp < 10.d0) THEN
          WRITE(name1, 101) TRIM(prefix), '.', cname, '_aniso_00', temp
        ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
          WRITE(name1, 102) TRIM(prefix), '.', cname, '_aniso_0', temp
        ELSEIF (temp >= 100.d0) THEN
          WRITE(name1, 103) TRIM(prefix), '.', cname, '_aniso_', temp
        ENDIF
      OPEN(UNIT = iufilgap, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('eliashberg_write_cont_raxis', 'error opening file ' // name1, iufilgap)
        WRITE(iufilgap, '(6a20)') '#        w [eV]', 'Enk-Ef [eV]', 'Re[znorm(w)]', 'Im[znorm(w)]',&
                                 'Re[delta(w)] [eV]', 'Im[delta(w)] [eV]'
      ENDIF
      !
      DO ik = 1, nkfs
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd,ik) - ef0) < fsthick) THEN
            lgap = .TRUE.
            ! DO iw = 1, nsw
            DO iw = 1, nsw - 1   ! FG: this change is to prevent segfault in ws(iw+1) and adelta(*,*,iw+1)
              var1 = REAL(adelta(ibnd, ik, iw))
              var2 = REAL(adelta(ibnd, ik, iw + 1))
              var3 = var1 - ws(iw)
              var4 = var2 - ws(iw + 1)
              IF (lgap .AND. iw < nqstep .AND. var1 > 0.d0 .AND. var2 > 0.d0 .AND. var3 * var4 < 0.d0) THEN
                agap(ibnd, ik, itemp) = (var3 * ws(iw + 1) - var4 * ws(iw)) / (var3 - var4)
                lgap = .FALSE.
              ENDIF
              IF (iverbosity == 2) THEN
                WRITE(iufilgap, '(6ES20.10)') ws(iw), ekfs(ibnd, ik) - ef0, &
                      REAL(aznorm(ibnd, ik, iw)), AIMAG(aznorm(ibnd, ik, iw)), &
                      REAL(adelta(ibnd, ik, iw)), AIMAG(adelta(ibnd, ik, iw))
              ENDIF
            ENDDO ! iw
            IF (lgap) &
              agap(ibnd,ik,itemp) = REAL(adelta(ibnd,ik,1))
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
      IF (iverbosity == 2) &
        CLOSE(iufilgap)
      !
      CALL gap_distribution_FS(itemp, cname)
      !
    ENDIF
    !
    ! isotropic case
    ! SP: Only write isotropic for laniso if user really wants that
    IF ((laniso .AND. iverbosity == 2) .OR. liso) THEN
      IF (temp < 10.d0) THEN
        WRITE(name1, 104) TRIM(prefix), '.', cname, '_iso_00', temp
      ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
        WRITE(name1, 105) TRIM(prefix), '.', cname, '_iso_0', temp
      ELSEIF (temp >= 100.d0) THEN
        WRITE(name1, 106) TRIM(prefix), '.', cname, '_iso_', temp
      ENDIF
      OPEN(UNIT = iufilgap, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('eliashberg_write_cont_raxis', 'error opening file ' // name1, iufilgap)
      !
      WRITE(iufilgap,'(5a20)') 'w [eV]', 'Re[znorm(w)]', 'Im[znorm(w)]', 'Re[delta(w)] [eV]', 'Im[delta(w)] [eV]'
      lgap = .TRUE.
      ! DO iw = 1, nsw
      DO iw = 1, nsw-1   ! this change is to prevent segfault in delta(iw+1) and ws(iw+1)
        var1 = REAL(delta(iw))
        var2 = REAL(delta(iw + 1))
        var3 = var1 - ws(iw)
        var4 = var2 - ws(iw + 1)
        IF (lgap .AND. iw < nqstep .AND. var1 > 0.d0 .AND. var2 > 0.d0 .AND. var3 * var4 < 0.d0) THEN
            gap(itemp) = (var3 * ws(iw + 1) - var4 * ws(iw)) / (var3 - var4)
          lgap = .FALSE.
        ENDIF
        WRITE(iufilgap, '(5ES20.10)') ws(iw), REAL(znorm(iw)), AIMAG(znorm(iw)), &
                                      REAL(delta(iw)), AIMAG(delta(iw))
      ENDDO ! iw
      CLOSE(iufilgap)
      IF (lgap ) &
        gap(itemp) = REAL(delta(1))
    ENDIF
    !
    101 FORMAT(a, a1, a4, a9, f4.2)
    102 FORMAT(a, a1, a4, a8, f5.2)
    103 FORMAT(a, a1, a4, a7, f6.2)
    104 FORMAT(a, a1, a4, a7, f4.2)
    105 FORMAT(a, a1, a4, a6, f5.2)
    106 FORMAT(a, a1, a4, a5, f6.2)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE eliashberg_write_cont_raxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_a2f()
    !-----------------------------------------------------------------------
    !!
    !! Read the eliashberg spectral function from fila2f
    !!
    USE epwcom,        ONLY : nqstep, fila2f
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
    IF (fila2f == ' ') WRITE(fila2f, '(a, a8)') TRIM(prefix), '.a2f_iso'
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
    USE elph2,     ONLY : nqtotf, wf, wqf, xqf
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
    INTEGER :: totq
    !! Total number of q-points inside fsthick
    INTEGER :: nqf1_, nqf2_, nqf3_
    !! Temporary variable for number of q-points along each direction
    INTEGER, ALLOCATABLE :: selecq(:)
    !! List of selected q-points
    !
    IF (mpime == ionode_id) THEN
      ! read 'selecq.fmt' file
      OPEN(UNIT = iunselecq, FILE = 'selecq.fmt', STATUS = 'old', IOSTAT = ios)
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
      OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', FORM = 'unformatted', IOSTAT = ios)
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
    DEALLOCATE(selecq, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_frequencies', 'Error deallocating selecq', 1)
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
    USE eliashbergcom, ONLY : nkfs, nbndfs, dosef, ef0, ekfs, wkfs, xkfs, w0g
    USE constants_epw, ONLY : ryd2ev, zero
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
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    !
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
      OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', FORM = 'unformatted', IOSTAT = ios)
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
      !
      ! at each k-point keep only the bands within the Fermi shell
      !
      ALLOCATE(ekf_(nbnd_, nkfs), STAT = ierr)
      IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error allocating ekf_', 1)
      ekf_(:, :) = zero
      !
      ! nbndfs - nr of bands within the Fermi shell
      !
      nbndfs = 0
      DO ik = 1, nkfs ! loop over irreducible k-points
        !READ(iufilegnv,'(4f15.9)') wkfs(ik), xkfs(:, ik)
        READ(iufilegnv) wkfs(ik), xkfs(:, ik)
        DO ibnd = 1, nbnd_
          !READ(iufilegnv, '(ES20.10)') ekf_(ibnd, ik)
          READ(iufilegnv) ekf_(ibnd, ik)
        ENDDO
        n = 0
        DO ibnd = 1, nbnd_
          ! go from Ryd to eV
          ekf_(ibnd, ik) = ekf_(ibnd, ik) * ryd2ev
          IF (ABS(ekf_(ibnd, ik) - ef0) < fsthick) THEN
            n = n + 1
            IF (nbndfs < n) nbndfs = n
          ENDIF
        ENDDO
      ENDDO
      WRITE(stdout, '(5x, i7, a/)') nbndfs, ' bands within the Fermi window'
      CLOSE(iufilegnv)
      !
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    CALL mp_bcast(nbndfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(wkfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(xkfs, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    ALLOCATE(ekfs(nbndfs, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error allocating ekfs', 1)
    ALLOCATE(w0g(nbndfs, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error allocating w0g', 1)
    ! sanity choice
    ekfs(:, :) = ef0 - 10.d0 * fsthick
    w0g(:, :) = zero
    IF (mpime == ionode_id) THEN
      DO ik = 1, nkfs ! loop over k-points
        n = 0
        DO ibnd = 1, nbnd_
          IF (ABS(ekf_(ibnd, ik) - ef0 ) < fsthick) THEN
            n = n + 1
            ekfs(n, ik) = ekf_(ibnd, ik)
            w0g(n, ik) = w0gauss((ekfs(n, ik) - ef0) / degaussw, 0) / degaussw
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE(ekf_, STAT = ierr)
      IF (ierr /= 0) CALL errore('read_eigenvalues', 'Error deallocating ekf_', 1)
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    CALL mp_bcast(ekfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(w0g, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout,'(/5x,a/)') 'Finish reading egnv file '
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
    !
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout, ionode_id
    USE io_var,    ONLY : iufilikmap, iunselecq
    USE io_files,  ONLY : prefix, tmp_dir
    USE symm_base, ONLY : t_rev, time_reversal, s, set_sym_bl
    USE modes,     ONLY : nmodes
    USE epwcom,    ONLY : nkf1, nkf2, nkf3, nqstep
    USE elph2,     ONLY : nqtotf, xqf
    USE grid,     ONLY : kpmq_map
    USE eliashbergcom, ONLY : ixkf, ixkff, xkff, xkfs, nkfs, ixkqf, ixqfs, nbndfs, nqfs, memlt_pool
    USE constants_epw, ONLY : zero
    USE symm_base, ONLY : nrot
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
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    INTEGER :: imelt
    !! Memory allocated
    INTEGER, ALLOCATABLE :: index_(:, :)
    !! Index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell
    !
    REAL(KIND = DP) :: xk(3)
    !! coordinates of k points
    REAL(KIND = DP) :: xq(3)
    !! coordinates of q points
    !
    INTEGER :: iqq
    !! Q-point index from selecq.fmt window
    INTEGER :: totq
    !! Total number of q-points inside fsthick
    INTEGER :: nqtot
    !! Total number of q-point for verifications
    INTEGER, ALLOCATABLE :: selecq(:)
    !! List of selected q-points
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
    ixkff(:) = 0
    !
    IF (mpime == ionode_id) THEN
      ! read 'selecq.fmt' file
      OPEN(UNIT = iunselecq, FILE = 'selecq.fmt', STATUS = 'old', IOSTAT = ios)
      READ(iunselecq, *) totq
      ALLOCATE(selecq(totq), STAT = ierr)
      IF (ierr /= 0) CALL errore('read_kqmap', 'Error allocating selecq', 1)
      selecq(:) = 0
      READ(iunselecq, *) nqtot
      READ(iunselecq, *) selecq(:)
      CLOSE(iunselecq)
      !
      dirname = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat'
      filikmap = TRIM(dirname) // '/' // 'ikmap'
      !OPEN(UNIT = iufilikmap, FILE = filikmap, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufilikmap, FILE = filikmap, STATUS = 'unknown', FORM = 'unformatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_kqmap', 'error opening file ' // filikmap, iufilikmap)
      !
      !READ(iufilikmap, *) ixkff(1:nkftot)
      READ(iufilikmap) ixkff(1:nkftot)
      !
      CLOSE(iufilikmap)
    ENDIF
    !
    CALL mp_bcast(totq, ionode_id, inter_pool_comm)
    IF (mpime /= ionode_id) ALLOCATE(selecq(totq))
    CALL mp_bcast(selecq, ionode_id, inter_pool_comm)
    CALL mp_bcast(ixkff, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ! get the size of required memory for ixkqf, nqfs, index_
    imelt = (nqtotf + 1) * nkfs + (upper_bnd - lower_bnd + 1) * nqtotf
    CALL mem_size_eliashberg(1, imelt)
    !
    ALLOCATE(ixkqf(nkfs, nqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_kqmap', 'Error allocating ixkqf', 1)
    ALLOCATE(nqfs(nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_kqmap', 'Error allocating nqfs', 1)
    ALLOCATE(index_(lower_bnd:upper_bnd, nqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_kqmap', 'Error allocating index_', 1)
    ixkqf(:, :) = 0
    nqfs(:) = 0
    index_(:, :) = 0
    !
    ! find the index of k+sign*q on the fine k-mesh
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
    CALL mp_sum(ixkqf, inter_pool_comm)
    CALL mp_sum(nqfs,  inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    ! get the size of required memory for ixqfs
    imelt = nkfs * MAXVAL(nqfs(:))
    CALL mem_size_eliashberg(1, imelt)
    !
    ALLOCATE(ixqfs(nkfs, MAXVAL(nqfs(:))), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_kqmap', 'Error allocating ixqfs', 1)
    ixqfs(:, :) = 0
    !
    DO ik = lower_bnd, upper_bnd
      DO iq = 1, nqfs(ik)
        !
        ! ixqfs - index q-point on the full q-mesh for which k+sign*q is within the Fermi shell
        !
        ixqfs(ik,iq) = index_(ik,iq)
      ENDDO
    ENDDO
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum(ixqfs, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    DEALLOCATE(index_, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_kqmap', 'Error deallocating index_', 1)
    !
    ! remove memory allocated for index_
    imelt = nqtotf * (upper_bnd - lower_bnd + 1)
    CALL mem_size_eliashberg(1, -imelt)
    !
    DEALLOCATE(selecq, STAT = ierr)
    IF (ierr /= 0) CALL errore('read_kqmap', 'Error allocating selecq', 1)
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
    USE io_var,        ONLY : iufileph
    USE io_files,      ONLY : prefix, tmp_dir
    USE modes,         ONLY : nmodes
    USE elph2,         ONLY : nqtotf, wf
    USE epwcom,        ONLY : eps_acustic, fsthick
    USE eliashbergcom, ONLY : nkfs, nbndfs, ef0, ekfs, g2, ixkqf, nqfs
    USE constants_epw, ONLY : ryd2ev, zero
    USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, npool
    USE division,      ONLY : fkbounds
    USE low_lvl,       ONLY : set_ndnmbr, mem_size_eliashberg
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256) :: filephmat
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
    INTEGER :: iq
    !! Counter on q-points
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
    INTEGER :: imelt
    !! Memory allocated
    INTEGER :: nmin, nmax
    !! Lower/upper bound index for .ephmat file read in current pool
    INTEGER :: nks
    !! Counter on k points within the Fermi shell
    INTEGER :: nnk
    !! Number of k-points within the Fermi shell
    INTEGER :: nnq(nkfs)
    !! Number of k+q points within the Fermi shell for a given k-point
    INTEGER :: nkpool(npool)
    !! nkpool(ipool) - sum of nr. of k points from pool 1 to pool ipool
    !
    REAL(KIND = DP) :: gmat
    !! Electron-phonon matrix element square
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
    DO ipool = 1, npool ! nr of pools
      CALL set_ndnmbr(0, ipool, 1, npool, filelab)
#if defined(__MPI)
      filephmat = TRIM(dirname) // '/' // 'ephmat' // filelab
#else
      filephmat = TRIM(dirname) // '/' // 'ephmat'
#endif
      !OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', FORM = 'unformatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_ephmat', 'error opening file ' // filephmat, iufileph)
      !READ(iufileph, '(2i7)') tmp_pool_id, nkpool(ipool)
      READ(iufileph) tmp_pool_id, nkpool(ipool)
      IF (ipool /= tmp_pool_id)  CALL errore('read_ephmat', &
          'npool should be equal to the number of .ephmat files', 1)
      IF (ipool > 1) &
        nkpool(ipool) = nkpool(ipool) + nkpool(ipool - 1)
      !WRITE(stdout, '(2i7)') tmp_pool_id, nkpool(ipool)
      CLOSE(iufileph)
    ENDDO
    CALL mp_barrier(inter_pool_comm)
    !
    ! since the nkfs k-points within the Fermi shell are not evenly distrubed
    ! among the .ephmat files, we re-distribute them here among the npool-pools
    nmin = npool
    nmax = npool
    DO ipool = npool, 1, -1
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
    DO ipool = 1, npool ! nr of pools
      CALL set_ndnmbr(0, ipool, 1, npool, filelab)
#if defined(__MPI)
      filephmat = TRIM(dirname) // '/' // 'ephmat' // filelab
#else
      filephmat = TRIM(dirname) // '/' // 'ephmat'
#endif
      !OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', FORM ='unformatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_ephmat', 'error opening file ' // filephmat, iufileph)
      !READ(iufileph, '(2i7)') tmp_pool_id, nks
      READ(iufileph) tmp_pool_id, nks
      IF (ipool >= nmin .AND. ipool <= nmax) THEN
        DO iq = 1, nqtotf ! loop over q-points
          DO ik = 1, nks ! loop over k-points in the pool
            IF (ixkqf(ik+nnk, iq) > 0) THEN
              nnq(ik + nnk) = nnq(ik + nnk) + 1
              DO imode = 1, nmodes ! loop over phonon modes
                DO ibnd = 1, nbndfs ! loop over iband's
                  IF (ABS(ekfs(ibnd, ik+nnk) - ef0) < fsthick) THEN
                    DO jbnd = 1, nbndfs ! loop over jband's
                      IF (ABS(ekfs(jbnd, ixkqf(ik + nnk, iq)) - ef0) < fsthick) THEN
                        !READ(iufileph, '(ES20.10)') gmat
                        READ(iufileph) gmat
                        IF (ik+nnk >= lower_bnd .AND. ik+nnk <= upper_bnd) THEN
                          ! go from Ryd to eV
                          IF (wf(imode, iq) > eps_acustic) THEN
                            g2(ik + nnk, nnq(ik + nnk), ibnd, jbnd, imode) = gmat * ryd2ev * ryd2ev
                          ELSE
                            g2(ik + nnk, nnq(ik + nnk), ibnd, jbnd, imode) = zero
                          ENDIF
                        ENDIF
                      ENDIF ! ekq
                    ENDDO ! jbnd
                  ENDIF ! ekk
                ENDDO ! ibnd
              ENDDO ! imode
            ENDIF ! ekk and ekq
          ENDDO ! ik
        ENDDO ! iq
        CLOSE(iufileph)
      ENDIF ! ipool
      nnk = nnk + nks
      IF (ipool == npool .AND. nnk /= nkfs)  CALL errore('read_ephmat', &
          'nnk should be equal to nkfs', 1)
    ENDDO ! ipool
    !
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout, '(/5x, a/)') 'Finish reading .ephmat files'
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE read_ephmat
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE write_ephmat(iqq, iq, totq)
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
                           fermi_energy
    USE pwcom,      ONLY : ef
    USE elph2,      ONLY : etf, ibndmin, ibndmax, nkqf, epf17, wkf, nkf, &
                           nqtotf, wf, xqf, nkqtotf, efnew, nbndfst, nktotf
    USE eliashbergcom, ONLY : nkfs, ekfs, wkfs, xkfs, dosef, ixkf, ixkqf, nbndfs
    USE constants_epw, ONLY : ryd2ev, ryd2mev, two, eps8
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id, npool
    USE division,      ONLY : fkbounds
    USE low_lvl,       ONLY : set_ndnmbr
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Q-point index from selecq.fmt window
    INTEGER, INTENT(in) :: iq
    !! Q-point index from full grid
    INTEGER :: totq
    !! Total number of q-points inside fsthick
    !
    ! Local variables
    !
    CHARACTER(LEN = 256) :: filfreq
    !! Name freq file
    CHARACTER(LEN = 256) :: filegnv
    !! Name eigenvalue file
    CHARACTER(LEN = 256) :: filephmat
    !! Name e-ph mat file
    CHARACTER(LEN = 4) :: filelab
    !!
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ibnd, jbnd
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
    INTEGER :: ind(npool)
    !! Temporary index
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
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! Function to compute the density of states at the Fermi level
    REAL(KIND = DP), EXTERNAL :: efermig
    !! Return the fermi energy
    INTEGER :: dummy
    !! Dummy variable for writing
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory to save ikmap/egnv/freq/ephmat files
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
      IF (iq == 1) THEN
        !OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
        OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', FORM = 'unformatted', IOSTAT = ios)

      ELSE
        !OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', POSITION = 'append', FORM = 'formatted', IOSTAT = ios)
        OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', POSITION = 'append', FORM = 'unformatted', IOSTAT = ios)
      ENDIF
      IF (ios /= 0) CALL errore('write_ephmat', 'error opening file ' // filfreq, iufilfreq)
      !IF (iq == 1) WRITE(iufilfreq, '(5i7)') nqtotf, nqf1, nqf2, nqf3, nmodes
      IF (iq == 1) WRITE(iufilfreq) nqtotf, nqf1, nqf2, nqf3, nmodes
      !WRITE(iufilfreq, '(3f15.9)') xqf(:, iq)
      WRITE(iufilfreq) xqf(:, iq)
      DO imode = 1, nmodes
        !WRITE(iufilfreq, '(ES20.10)') wf(imode, iq)
        WRITE(iufilfreq) wf(imode, iq)
      ENDDO
      CLOSE(iufilfreq)
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
    IF (iq == 1) THEN
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
      ! write eigenvalues to file
      IF (my_pool_id == 0) THEN
        filegnv = TRIM(dirname) // '/' // 'egnv'
        !OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
        OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', FORM = 'unformatted', IOSTAT = ios)
        IF (ios /= 0) CALL errore('write_ephmat', 'error opening file ' // filegnv, iufilegnv)
        IF (nks /= nkfs) CALL errore('write_ephmat', &
          'nks should be equal to nr. of irreducible k-points within the Fermi shell on the fine mesh', 1)
        !WRITE(iufilegnv, '(5i7)') nkftot, nkf1, nkf2, nkf3, nks
        !WRITE(iufilegnv, '(i7,5ES20.10)') nbndfst, ef, ef0, dosef, degaussw,fsthick
        WRITE(iufilegnv) nkftot, nkf1, nkf2, nkf3, nks
        WRITE(iufilegnv) nbndfst, ef, ef0, dosef, degaussw, fsthick
        DO ik = 1, nks
          !WRITE(iufilegnv, '(4f15.9)') wkfs(ik), xkfs(:, ik)
          WRITE(iufilegnv) wkfs(ik), xkfs(:, ik)
          DO ibnd = 1, nbndfst
            !WRITE(iufilegnv, '(ES20.10)') ekfs(ibnd, ik)
            WRITE(iufilegnv) ekfs(ibnd, ik)
          ENDDO
        ENDDO
        CLOSE(iufilegnv)
      ENDIF
      !
    ENDIF ! iq
    !
    ! write the e-ph matrix elements in the Bloch representation on the fine mesh
    ! in .ephmat files (one for each pool)
    !
#if defined(__MPI)
    CALL set_ndnmbr(0, my_pool_id + 1, 1, npool, filelab)
    filephmat = TRIM(dirname) // '/' // 'ephmat' // filelab
#else
    filephmat = TRIM(dirname) // '/' // 'ephmat'
#endif
    !
    IF (iq == 1) THEN
      !OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', FORM = 'unformatted', IOSTAT = ios)
    ELSE
      !OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', POSITION = 'append', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufileph, FILE = filephmat, STATUS = 'unknown', POSITION = 'append', FORM = 'unformatted', IOSTAT = ios)
    ENDIF
    IF (ios /= 0) CALL errore('write_ephmat', 'error opening file ' // filephmat, iufileph)
    !
    !IF (iq == 1) WRITE(iufileph, '(2i7)') my_pool_id+1, fermicount
    IF (iq == 1) WRITE(iufileph) my_pool_id+1, fermicount
    !
    ! nkf - nr of k-points in the pool (fine mesh)
    ! for mp_mesh_k = true nkf is nr of irreducible k-points in the pool
    !
    DO ik = 1, nkf
      !
      ikk = 2 * ik - 1
      !
      ! go only over irreducible k-points
      !
      ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
      !
      DO imode = 1, nmodes ! phonon modes
        wq = wf(imode, iq)
        inv_wq =  1.0 / (two * wq)
        !
        DO ibnd = 1, nbndfst
          IF (ABS(ekfs(ibnd, ixkf(lower_bnd + ik - 1)) - ef0) < fsthick) THEN
            DO jbnd = 1, nbndfst
              IF (ABS(ekfs(jbnd, ixkqf(ixkf(lower_bnd + ik - 1), iq)) - ef0) < fsthick) THEN
                !
                ! here we take into account the zero-point DSQRT(hbar/2M\omega)
                ! with hbar = 1 and M already contained in the eigenmodes
                ! g2 is Ry^2, wkf must already account for the spin factor
                !
                IF (shortrange .AND. (ABS(xqf(1, iq)) > eps8 .OR. ABS(xqf(2, iq)) > eps8 &
                     .OR. ABS(xqf(3, iq)) > eps8 )) THEN
                  ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary
                  !     number, in which case its square will be a negative number.
                  g2 = REAL((epf17(jbnd, ibnd, imode, ik)**two) * inv_wq)
                ELSE
                  g2 = ABS(epf17(jbnd, ibnd, imode, ik))**two * inv_wq
                ENDIF
                ind(my_pool_id + 1) = ind(my_pool_id + 1) + 1
                tmp_g2(ind(my_pool_id + 1)) = g2
              ENDIF
            ENDDO ! jbnd
          ENDIF
        ENDDO ! ibnd
      ENDDO ! imode
      !
    ENDDO ! ik's
    !
    IF (ind(my_pool_id + 1) > 0) THEN
      DO ifil = 1, ind(my_pool_id + 1)
        !WRITE(iufileph, '(ES20.10)') tmp_g2(ifil)
        WRITE(iufileph) tmp_g2(ifil)
      ENDDO
    ENDIF
    CLOSE(iufileph)
    !
    IF (my_pool_id == 0) THEN
      ! format is compatible with IBTE
      OPEN(UNIT = iunrestart, FILE = 'restart.fmt')
      WRITE(iunrestart, *) iqq
      WRITE(iunrestart, *) dummy
      WRITE(iunrestart, *) dummy
      WRITE(iunrestart, *) npool
      DO ipool = 1, npool
        WRITE(iunrestart, *) dummy
      ENDDO
      DO ipool = 1, npool
       WRITE(iunrestart, *) dummy
      ENDDO
      CLOSE(iunrestart)
    ENDIF
    !
    IF (iqq == totq) THEN
      DEALLOCATE(ekfs, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_ephmat', 'Error deallocating ekfs', 1)
      DEALLOCATE(wkfs, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_ephmat', 'Error deallocating wkfs', 1)
      DEALLOCATE(xkfs, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_ephmat', 'Error deallocating xkfs', 1)
      DEALLOCATE(ixkqf, STAT = ierr)
      IF (ierr /= 0) CALL errore('write_ephmat', 'Error deallocating ixkqf', 1)
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
    USE eliashbergcom, ONLY : nkfs, ixkf, xkfs, wkfs, ekfs, nbndfs
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
      wkf_(lower_bnd + nk - 1)   = wkf(ikk)
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
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error allocating wkf_', 1)
    ALLOCATE(xkfs(3, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('kmesh_fine', 'Error allocating xkfs', 1)
    xkfs(:, :) = zero
    wkfs(:) = zero
    ekfs(:, :) = zero
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
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    CALL mp_bcast(ixkf, ionode_id, inter_pool_comm)
    CALL mp_bcast(xkfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(wkfs, ionode_id, inter_pool_comm)
    CALL mp_bcast(ekfs, ionode_id, inter_pool_comm)
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
    USE epwcom,    ONLY : nkf1, nkf2, nkf3, mp_mesh_k, fixsym
    USE elph2,     ONLY : nqtotf, nktotf, xqf, map_rebal
    USE eliashbergcom, ONLY : ixkff, ixkf, xkfs, nkfs, ixkqf, ixqfs, nqfs
    USE constants_epw, ONLY : eps5, zero
    USE symm_base, ONLY : nrot
    USE io_global, ONLY : stdout, ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    USE division,  ONLY : fkbounds
    USE kinds_epw, ONLY : SIK2
    USE grid,      ONLY : kpmq_map, kpoint_grid_epw
    USE io_files,  ONLY : prefix, tmp_dir, create_directory
    USE io_var,    ONLY : iufilikmap
    USE low_lvl,   ONLY : fix_sym
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
    INTEGER :: nks
    !! Number of non-equivalent k points
    INTEGER :: ns
    !! Counter on rotation operations
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios
    !! IO error message
    INTEGER, ALLOCATABLE :: index_(:, :)
    !! Index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell
    !
    REAL(KIND = DP) :: xk(3)
    !! coordinates of k points
    REAL(KIND = DP) :: xq(3)
    !! coordinates of q points
    !
    INTEGER :: bztoibz(nkf1 * nkf2 * nkf3)
    !! BZ to IBZ mapping
    INTEGER :: bztoibz_tmp(nkf1 * nkf2 * nkf3)
    !! Temporary mapping
    INTEGER(SIK2) :: s_bztoibz(nkf1 * nkf2 * nkf3)
    !! symmetry matrix for each k-point from the full BZ
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
      bztoibz(:)   = 0
      s_bztoibz(:) = 0
      !
      CALL set_sym_bl()
      IF (fixsym) CALL fix_sym(.TRUE.)
      CALL kpoint_grid_epw(nrot, time_reversal, .FALSE., s, t_rev, nkf1, nkf2, nkf3, bztoibz, s_bztoibz)
      ! 
      bztoibz_tmp(:) = 0
      DO ikbz = 1, nkftot
        bztoibz_tmp(ikbz) = map_rebal(bztoibz(ikbz))
        ixkff(ikbz) = ixkf(map_rebal(bztoibz(ikbz)))
      ENDDO
      bztoibz(:) = bztoibz_tmp(:)
      !
    ELSE
      ! full k-point grid
      DO nk = 1, nkftot
        ixkff(nk) = ixkf(nk)
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
      OPEN(UNIT = iufilikmap, FILE = filikmap, STATUS = 'unknown', FORM = 'unformatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('kqmap_fine', 'error opening file ' // filikmap, iufilikmap)
      !
      !WRITE(iufilikmap, *) ixkff(1:nkftot)
      WRITE(iufilikmap) ixkff(1:nkftot)
      !
      CLOSE(iufilikmap)
    ENDIF
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ALLOCATE(ixkqf(nkfs, nqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('kqmap_fine', 'Error allocating ixkqf', 1)
    ALLOCATE(nqfs(nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('kqmap_fine', 'Error allocating nqfs', 1)
    ALLOCATE(index_(lower_bnd:upper_bnd, nqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('kqmap_fine', 'Error allocating index_', 1)
    ixkqf(:, :) = 0
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
        ! find ixkqf(ik,iq) - index of k+sign*q on the fine k-mesh within the Fermi shell
        !
        ixkqf(ik, iq) = ixkff(nkq)
        !
        ! nqfs(ik) - nr of q-points at each k-point for which k+sign*q is within the Fermi shell
        ! index_   - index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell
        !
        IF (ixkqf(ik, iq) > 0) THEN
          nqfs(ik) = nqfs(ik) + 1
          index_(ik, nqfs(ik)) = iq
        ENDIF
      ENDDO ! loop over full set of q-points (fine mesh)
    ENDDO ! loop over k-points within the Fermi shell in each pool (fine mesh)
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum(ixkqf, inter_pool_comm)
    CALL mp_sum(nqfs, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    ALLOCATE(ixqfs(nkfs, MAXVAL(nqfs(:))), STAT = ierr)
    IF (ierr /= 0) CALL errore('kqmap_fine', 'Error allocating ixqfs', 1)
    ixqfs(:, :) = 0
    !
    DO ik = lower_bnd, upper_bnd
      DO iq = 1, nqfs(ik)
        !
        ! ixqfs - index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell
        !
        ixqfs(ik, iq) = index_(ik, iq)
      ENDDO
    ENDDO
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum(ixqfs, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    DEALLOCATE(ixkff, STAT = ierr)
    IF (ierr /= 0) CALL errore('kqmap_fine', 'Error deallocating ixkff', 1)
    DEALLOCATE(ixqfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('kqmap_fine', 'Error deallocating ixqfs', 1)
    DEALLOCATE(index_, STAT = ierr)
    IF (ierr /= 0) CALL errore('kqmap_fine', 'Error deallocating index_', 1)
    DEALLOCATE(nqfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('kqmap_fine', 'Error deallocating nqfs', 1)
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
    SUBROUTINE check_restart_ephwrite()
    !-----------------------------------------------------------------------
    !!
    !!   This routine checks the variables in restart while writing ephmat 
    !!   6/28/2020 Hari Paudyal
    !!
    USE io_files,  ONLY : prefix, tmp_dir
    USE epwcom,    ONLY : nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, fsthick, mp_mesh_k
    USE io_var,     ONLY : iufilfreq, iufilegnv
    !
    IMPLICIT NONE
    !
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
    LOGICAL :: exst
    !! Logical for existence of files
    LOGICAL :: exst2
    !! Logical for existence of files
    INTEGER :: ios
    !! IO error message
    CHARACTER(LEN = 256) :: filfreq
    !! file name
    CHARACTER(LEN = 256) :: filegnv
    !! file name
    CHARACTER(LEN = 256) :: dirname
    !! Name of the directory to save ikmap/egnv/freq/ephmat files
    !
    dirname = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat'
    !
    INQUIRE(FILE = 'restart.fmt', EXIST = exst)
    !
    IF (exst) THEN
      INQUIRE(FILE = TRIM(dirname) // '/' // 'ikmap', EXIST = exst2)
      IF (.NOT. exst2) THEN
        CALL errore('check_restart_ephwrite', 'A restart.fmt is present but the directory ' // TRIM(prefix) // '.ephmat' // &
                    ' is not found. Remove the restart.fmt file and restart.', 1)
      ENDIF
      !
      ! read header of egnv file
      filegnv = TRIM(dirname) // '/' // 'egnv'
      !OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufilegnv, FILE = filegnv, STATUS = 'unknown', FORM = 'unformatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('check_restart_ephwrite', 'error opening file '//filegnv, iufilegnv)
      !
      !READ(iufilegnv, '(5i7)') nkftot_, nkf1_, nkf2_, nkf3_, nkfs_
      READ(iufilegnv) nkftot_, nkf1_, nkf2_, nkf3_, nkfs_
      IF (nkf1 /= nkf1_ .OR. nkf2 /= nkf2_ .OR. nkf3 /= nkf3_) &
        CALL errore('check_restart_ephwrite', 'nkf1, nkf2, nkf3 is not consistent with restart.fmt', 1)
      CLOSE(iufilegnv)
      !
      ! read header of freq file
      filfreq = TRIM(dirname) // '/' // 'freq'
      !OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
      OPEN(UNIT = iufilfreq, FILE = filfreq, STATUS = 'unknown', FORM = 'unformatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('kmap_fine', 'error opening file ' // filfreq, iufilfreq)
      !READ(iufilfreq, '(5i7)') nqtotf_, nqf1_, nqf2_, nqf3_, nmodes_
      READ(iufilfreq) nqtotf_, nqf1_, nqf2_, nqf3_, nmodes_
      IF (nqf1 /= nqf1_ .OR. nqf2 /= nqf2_ .OR. nqf3 /= nqf3_) &
        CALL errore('check_restart_ephwrite', 'nqf1, nqf2, nqf3 is not consistent with restart.fmt', 1)
      CLOSE(iufilfreq)
      !
    ENDIF
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
    USE epwcom,        ONLY : fsthick
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : agap, nkfs, nbndfs, ef0, ekfs, w0g
    USE constants_epw, ONLY : kelvin2eV, zero, eps5
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
    REAL(KIND = DP) :: delta_max
    !! Max value of superconducting gap
    REAL(KIND = DP) :: weight
    !! Variable for weight
    REAL(KIND = DP), ALLOCATABLE :: delta_k_bin(:)
    !! Histogram superconducting gap
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    temp = gtemp(itemp) / kelvin2eV
    !
    delta_max = 1.1d0 * MAXVAL(agap(:,:,itemp))
    nbin = NINT(delta_max / eps5) + 1
    dbin = delta_max / DBLE(nbin)
    ALLOCATE(delta_k_bin(nbin), STAT = ierr)
    IF (ierr /= 0) CALL errore('gap_distribution_FS', 'Error allocating delta_k_bin', 1)
    delta_k_bin(:) = zero
    !
    DO ik = 1, nkfs
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          ibin = NINT(agap(ibnd, ik, itemp) / dbin) + 1
          weight = w0g(ibnd, ik)
          delta_k_bin(ibin) = delta_k_bin(ibin) + weight
        ENDIF
      ENDDO
    ENDDO
    !
    IF (temp < 10.d0) THEN
       WRITE(name1, 101) TRIM(prefix), '.', cname, '_aniso_gap0_00', temp
    ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
       WRITE(name1, 102) TRIM(prefix), '.', cname, '_aniso_gap0_0', temp
    ELSEIF (temp >= 100.d0) THEN
       WRITE(name1, 103) TRIM(prefix), '.', cname, '_aniso_gap0_', temp
    ENDIF
    !
    OPEN(UNIT = iufilgap, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
    IF (ios /= 0) CALL errore('gap_distribution_FS', 'error opening file ' // name1, iufilgap)
    DO ibin = 1, nbin
      WRITE(iufilgap,'(2ES20.10)') temp + delta_k_bin(ibin) / MAXVAL(delta_k_bin(:)), dbin * DBLE(ibin)
    ENDDO
    CLOSE(iufilgap)
    !
    DEALLOCATE(delta_k_bin, STAT = ierr)
    IF (ierr /= 0) CALL errore('gap_distribution_FS', 'Error deallocating delta_k_bin', 1)
    !
    101 FORMAT(a, a1, a4, a14, f4.2)
    102 FORMAT(a, a1, a4, a13, f5.2)
    103 FORMAT(a, a1, a4, a12, f6.2)
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
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : agap, nkfs, nbndfs, ef0, ekfs, ixkff
    USE constants_epw, ONLY : kelvin2eV, zero
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
    ALLOCATE(agap_tmp(nbndfs, 0:nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('gap_FS', 'Error allocating agap_tmp', 1)
    agap_tmp(:, 1:nkfs) = agap(:, 1:nkfs, itemp)
    agap_tmp(:, 0) = zero
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
             WRITE(name1, 101) TRIM(prefix), '.', cname, '_aniso_gap0_00', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 100.d0) THEN
             WRITE(name1, 102) TRIM(prefix), '.', cname, '_aniso_gap0_0', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 1000.d0) THEN
             WRITE(name1, 103) TRIM(prefix), '.', cname, '_aniso_gap0_', temp, '_', ibnd, '.cube'
          ENDIF
        ELSEIF (ibnd < 100) THEN
          IF (temp < 10.d0) THEN
             WRITE(name1, 104) TRIM(prefix), '.', cname, '_aniso_gap0_00', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 100.d0 .AND. temp > 9.9999d0) THEN
             WRITE(name1, 105) TRIM(prefix), '.', cname, '_aniso_gap0_0', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 1000.d0 .AND. temp > 99.9999d0) THEN
             WRITE(name1, 106) TRIM(prefix), '.', cname, '_aniso_gap0_', temp, '_', ibnd, '.cube'
          ENDIF
        ELSEIF (ibnd < 1000) THEN
          IF (temp < 10.d0) THEN
             WRITE(name1, 107) TRIM(prefix), '.', cname, '_aniso_gap0_00', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 100.d0 .AND. temp > 9.9999d0) THEN
             WRITE(name1, 108) TRIM(prefix), '.', cname, '_aniso_gap0_0', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 1000.d0 .AND. temp > 99.9999d0) THEN
             WRITE(name1, 109) TRIM(prefix), '.', cname, '_aniso_gap0_', temp, '_', ibnd, '.cube'
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
        WRITE(iufilgapFS, '(6f12.6)') (agap_tmp(ibnd, ixkff(ik)), ik = 1, nkf1 * nkf2 * nkf3)
        CLOSE(iufilgapFS)
      ENDDO
      !
    ENDIF
    !
    ! SP & RM : Write on file the superconducting gap close to the Fermi surface
    ! along with
    !     Cartesian coordinate, band index, energy distance from Fermi level and
    !     gap value.
    !
    IF (temp < 10.d0) THEN
      WRITE(name1, 110) TRIM(prefix), '.', cname, '_aniso_gap_FS_00', temp
    ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
      WRITE(name1, 111) TRIM(prefix), '.', cname, '_aniso_gap_FS_0', temp
    ELSEIF (temp >= 100.d0) THEN
      WRITE(name1, 112) TRIM(prefix), '.', cname, '_aniso_gap_FS_', temp
    ENDIF
    OPEN(UNIT = iufilgapFS, FILE = name1, STATUS = 'unknown', FORM = 'formatted', IOSTAT = ios)
    IF (ios /= 0) CALL errore('gap_FS', 'error opening file ' // name1, iufilgapFS)
    WRITE(iufilgapFS, '(a78)') '#               k-point                  Band Enk-Ef [eV]        delta(0) [eV]'
    DO i = 1, nkf1
      DO j = 1, nkf2
        DO k = 1, nkf3
          ik = k + (j - 1) * nkf3 + (i - 1) * nkf2 * nkf3
          !IF (ixkff(ik) > 0) THEN
            DO ibnd = 1, nbndfs
              ! RM: Everything is in eV here.
              ! SP: Here take a 0.2 eV interval around the FS.
              IF (ABS(ekfs(ibnd, ixkff(ik)) - ef0) < fsthick) THEN
              !IF (ABS(ekfs(ibnd, ixkff(ik)) - ef0 ) < 0.2) THEN
                 x1 = bg(1, 1) * (i - 1) /nkf1 + bg(1, 2) * (j - 1) / nkf2 + bg(1, 3) *(k - 1) / nkf3
                 x2 = bg(2, 1) * (i - 1) /nkf1 + bg(2, 2) * (j - 1) / nkf2 + bg(2, 3) *(k - 1) / nkf3
                 x3 = bg(3, 1) * (i - 1) /nkf1 + bg(3, 2) * (j - 1) / nkf2 + bg(3, 3) *(k - 1) / nkf3
                 WRITE(iufilgapFS,'(3f12.6, i8, f12.6, f24.15)') x1, x2, x3, ibnd, &
                       ekfs(ibnd, ixkff(ik)) - ef0, agap_tmp(ibnd, ixkff(ik))
              ENDIF
            ENDDO ! ibnd
          !ENDIF
        ENDDO  ! k
      ENDDO ! j
    ENDDO ! i
    CLOSE(iufilgapFS)
    !
    DEALLOCATE(agap_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('gap_FS', 'Error deallocating agap_tmp', 1)
    !
    101 FORMAT(a, a1, a4, a14, f4.2, a1, i1, a5)
    102 FORMAT(a, a1, a4, a13, f5.2, a1, i1, a5)
    103 FORMAT(a, a1, a4, a12, f6.2, a1, i1, a5)
    104 FORMAT(a, a1, a4, a14, f4.2, a1, i2, a5)
    105 FORMAT(a, a1, a4, a13, f5.2, a1, i2, a5)
    106 FORMAT(a, a1, a4, a12, f6.2, a1, i2, a5)
    107 FORMAT(a, a1, a4, a14, f4.2, a1, i3, a5)
    108 FORMAT(a, a1, a4, a13, f5.2, a1, i3, a5)
    109 FORMAT(a, a1, a4, a12, f6.2, a1, i3, a5)
    110 FORMAT(a, a1, a4, a16, f4.2)
    111 FORMAT(a, a1, a4, a15, f5.2)
    112 FORMAT(a, a1, a4, a14, f6.2)
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

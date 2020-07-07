  !
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !
  !----------------------------------------------------------------------
  MODULE io_selfen
  !----------------------------------------------------------------------
  !!
  !! This module contains various writing or reading routines related to self-energies.
  !! Most of them are for restart purposes.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !----------------------------------------------------------------------------
    SUBROUTINE selfen_el_write(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all)
    !----------------------------------------------------------------------------
    !!
    !! Write self-energy
    !!
    USE kinds,     ONLY : DP
    USE elph2,     ONLY : lower_bnd, upper_bnd, nbndfst
    USE io_var,    ONLY : iufilsigma_all
    USE io_files,  ONLY : diropn
    USE constants_epw, ONLY : zero
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : mpime
    USE io_global, ONLY : ionode_id
    USE epwcom,    ONLY : nstemp
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(inout) :: sigmar_all(nbndfst, nktotf, nstemp)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(inout) :: sigmai_all(nbndfst, nktotf, nstemp)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(inout) :: zi_all(nbndfst, nktotf, nstemp)
    !! Z parameter of electron-phonon self-energy accross all pools
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: lsigma_all
    !! Length of the vector
    INTEGER :: itemp
    !! Counter on temperatures
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf * nstemp + 2)
    !! Vector to store the array
    !
    IF (mpime == ionode_id) THEN
      !
      lsigma_all = 3 * nbndfst * nktotf * nstemp + 2
      ! First element is the current q-point
      aux(1) = REAL(iqq - 1, KIND = DP) ! we need to start at the next q
      ! Second element is the total number of q-points
      aux(2) = REAL(totq, KIND = DP)
      !
      i = 2
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            aux(i) = sigmar_all(ibnd, ik, itemp)
          ENDDO
        ENDDO
      ENDDO
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            aux(i) = sigmai_all(ibnd, ik, itemp)
          ENDDO
        ENDDO
      ENDDO
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            aux(i) = zi_all(ibnd, ik, itemp)
          ENDDO
        ENDDO
      ENDDO
      CALL diropn(iufilsigma_all, 'sigma_restart', lsigma_all, exst)
      CALL davcio(aux, lsigma_all, iufilsigma_all, 1, +1)
      CLOSE(iufilsigma_all)
    ENDIF
    !
    ! Make everythin 0 except the range of k-points we are working on
    IF (lower_bnd > 1) THEN
      sigmar_all(:, 1:lower_bnd - 1, :) = zero
      sigmai_all(:, 1:lower_bnd - 1, :) = zero
      zi_all(:, 1:lower_bnd - 1, :) = zero
    ENDIF
    IF (upper_bnd < nktotf) THEN
      sigmar_all(:, upper_bnd + 1:nktotf, :) = zero
      sigmai_all(:, upper_bnd + 1:nktotf, :) = zero
      zi_all(:, upper_bnd + 1:nktotf, :) = zero
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE selfen_el_write
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE selfen_el_read(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all)
    !----------------------------------------------------------------------------
    !!
    !! Self-energy reading
    !!
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE elph2,     ONLY : lower_bnd, upper_bnd, nbndfst
    USE io_var,    ONLY : iufilsigma_all
    USE io_files,  ONLY : prefix, tmp_dir, diropn
    USE constants_epw, ONLY :  zero
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    USE epwcom,    ONLY : nstemp
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: iqq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(out) :: sigmar_all(nbndfst, nktotf, nstemp)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: sigmai_all(nbndfst, nktotf, nstemp)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: zi_all(nbndfst, nktotf, nstemp)
    !! Z parameter of electron-phonon self-energy accross all pools
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: lsigma_all
    !! Length of the vector
    INTEGER :: nqtotf_read
    !! Total number of q-point read
    INTEGER :: itemp
    !! Counter on temperatures
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf * nstemp + 2)
    !! Vector to store the array
    !
    CHARACTER(LEN = 256) :: name1
    !
    IF (mpime == ionode_id) THEN
      !
      ! First inquire if the file exists
#if defined(__MPI)
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.sigma_restart1'
#else
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.sigma_restart'
#endif
      INQUIRE(FILE = name1, EXIST = exst)
      !
      IF (exst) THEN ! read the file
        !
        lsigma_all = 3 * nbndfst * nktotf * nstemp + 2
        CALL diropn(iufilsigma_all, 'sigma_restart', lsigma_all, exst)
        CALL davcio(aux, lsigma_all, iufilsigma_all, 1, -1)
        !
        ! First element is the iteration number
        iqq = INT(aux(1))
        iqq = iqq + 1 ! we need to start at the next q
        nqtotf_read = INT(aux(2))
        IF (nqtotf_read /= totq) CALL errore('selfen_el_read', &
          &'Error: The current total number of q-point is not the same as the read one. ', 1)
        !
        i = 2
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              sigmar_all(ibnd, ik, itemp) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              sigmai_all(ibnd, ik, itemp) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              zi_all(ibnd, ik, itemp) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iufilsigma_all)
      ENDIF
    ENDIF
    !
    CALL mp_bcast(exst, ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iqq, ionode_id, world_comm)
      CALL mp_bcast(sigmar_all, ionode_id, world_comm)
      CALL mp_bcast(sigmai_all, ionode_id, world_comm)
      CALL mp_bcast(zi_all, ionode_id, world_comm)
      !
      ! Make everythin 0 except the range of k-points we are working on
      IF (lower_bnd > 1) THEN
        sigmar_all(:, 1:lower_bnd - 1, :) = zero
        sigmai_all(:, 1:lower_bnd - 1, :) = zero
        zi_all(:, 1:lower_bnd - 1, :) = zero
      ENDIF
      IF (upper_bnd < nktotf) THEN
        sigmar_all(:, upper_bnd + 1:nktotf, :) = zero
        sigmai_all(:, upper_bnd + 1:nktotf, :) = zero
        zi_all(:, upper_bnd + 1:nktotf, :) = zero
      ENDIF
      !
      WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from: ', iqq,'/', totq
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE selfen_el_read
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE spectral_write(iqq, totq, nktotf, esigmar_all, esigmai_all)
    !----------------------------------------------------------------------------
    !!
    !! Write self-energy
    !!
    USE kinds,     ONLY : DP
    USE elph2,     ONLY : lower_bnd, upper_bnd, nbndfst
    USE epwcom,    ONLY : nstemp, wmin_specfun, wmax_specfun, nw_specfun
    USE io_var,    ONLY : iufilesigma_all
    USE io_files,  ONLY : diropn
    USE constants_epw, ONLY : zero
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : mpime
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(inout) :: esigmar_all(nbndfst, nktotf, nw_specfun, nstemp)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(inout) :: esigmai_all(nbndfst, nktotf, nw_specfun, nstemp)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: iw
    !! Counter on the frequency
    INTEGER :: lesigma_all
    !! Length of the vector
    INTEGER :: itemp
    !! Counter on temperature
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP) :: ww(nw_specfun)
    !! Current frequency
    REAL(KIND = DP) :: aux(2 * nbndfst * nktotf * nw_specfun * nstemp + 2)
    !! Vector to store the array
    !
    IF (mpime == ionode_id) THEN
      !
      ! energy range and spacing for spectral function
      !
      dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1.d0)
      DO iw = 1, nw_specfun
        ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
      ENDDO
      !
      lesigma_all = 2 * nbndfst * nktotf * nw_specfun * nstemp + 2
      ! First element is the current q-point
      aux(1) = REAL(iqq - 1, KIND = DP) ! we need to start at the next q
      ! Second element is the total number of q-points
      aux(2) = REAL(totq, KIND = DP)
      !
      i = 2
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            DO iw = 1, nw_specfun
              i = i + 1
              aux(i) = esigmar_all(ibnd, ik, iw, itemp)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            DO iw = 1, nw_specfun
              i = i + 1
              aux(i) = esigmai_all(ibnd, ik, iw, itemp)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL diropn(iufilesigma_all, 'esigma_restart', lesigma_all, exst)
      CALL davcio(aux, lesigma_all, iufilesigma_all, 1, +1)
      CLOSE(iufilesigma_all)
    ENDIF
    !
    ! Make everythin 0 except the range of k-points we are working on
    IF (lower_bnd > 1) THEN
      esigmar_all(:, 1:lower_bnd - 1, :, :) = zero
      esigmai_all(:, 1:lower_bnd - 1, :, :) = zero
    ENDIF
    IF (upper_bnd < nktotf) THEN
      esigmar_all(:, upper_bnd + 1:nktotf, :, :) = zero
      esigmai_all(:, upper_bnd + 1:nktotf, :, :) = zero
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE spectral_write
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE spectral_read(iqq, totq, nktotf, esigmar_all, esigmai_all)
    !----------------------------------------------------------------------------
    !!
    !! Self-energy reading
    !!
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE elph2,     ONLY : lower_bnd, upper_bnd, nbndfst
    USE epwcom,    ONLY : nstemp, wmin_specfun, wmax_specfun, nw_specfun
    USE io_var,    ONLY : iufilesigma_all
    USE io_files,  ONLY : prefix, tmp_dir, diropn
    USE constants_epw, ONLY : zero
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: iqq
    !! Current q-point
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    REAL(KIND = DP), INTENT(out) :: esigmar_all(nbndfst, nktotf, nw_specfun, nstemp)
    !! Real part of the electron-phonon self-energy accross all pools
    REAL(KIND = DP), INTENT(out) :: esigmai_all(nbndfst, nktotf, nw_specfun, nstemp)
    !! Imaginary part of the electron-phonon self-energy accross all pools
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Iterative index
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: iw
    !! Counter on the frequency
    INTEGER :: lesigma_all
    !! Length of the vector
    INTEGER :: nqtotf_read
    !! Total number of q-point read
    INTEGER :: itemp
    !! Counter on temperatures
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP) :: ww(nw_specfun)
    !! Current frequency
    REAL(KIND = DP) :: aux(2 * nbndfst * nktotf * nw_specfun * nstemp + 2)
    !! Vector to store the array
    !
    CHARACTER(LEN = 256) :: name1
    !
    !
    IF (mpime == ionode_id) THEN
      !
      ! First inquire if the file exists
#if defined(__MPI)
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.esigma_restart1'
#else
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.esigma_restart'
#endif
      INQUIRE(FILE = name1, EXIST = exst)
      !
      IF (exst) THEN ! read the file
        !
        lesigma_all = 2 * nbndfst * nktotf * nw_specfun * nstemp + 2
        CALL diropn(iufilesigma_all, 'esigma_restart', lesigma_all, exst)
        CALL davcio(aux, lesigma_all, iufilesigma_all, 1, -1)
        !
        ! First element is the iteration number
        iqq = INT(aux(1))
        iqq = iqq + 1 ! we need to start at the next q
        nqtotf_read = INT(aux(2))
        IF (nqtotf_read /= totq) CALL errore('electron_read',&
          &'Error: The current total number of q-point is not the same as the read one. ', 1)
        !
        i = 2
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              DO iw = 1, nw_specfun
                i = i + 1
                esigmar_all(ibnd, ik, iw, itemp) = aux(i)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              DO iw = 1, nw_specfun
                i = i + 1
                esigmai_all(ibnd, ik, iw, itemp) = aux(i)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iufilesigma_all)
      ENDIF
    ENDIF
    !
    CALL mp_bcast(exst, ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iqq, ionode_id, world_comm)
      CALL mp_bcast(esigmar_all, ionode_id, world_comm)
      CALL mp_bcast(esigmai_all, ionode_id, world_comm)
      !
      ! Make everythin 0 except the range of k-points we are working on
      IF (lower_bnd > 1) THEN
        esigmar_all(:, 1:lower_bnd - 1, :, :) = zero
        esigmai_all(:, 1:lower_bnd - 1, :, :) = zero
      ENDIF
      IF (upper_bnd < nktotf) THEN
        esigmar_all(:, upper_bnd + 1:nktotf, :, :) = zero
        esigmai_all(:, upper_bnd + 1:nktotf, :, :) = zero
      ENDIF
      !
      WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from: ', iqq,'/', totq
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE spectral_read
    !----------------------------------------------------------------------------
    !
  !------------------------------------------------------------------------------
  END MODULE io_selfen
  !------------------------------------------------------------------------------

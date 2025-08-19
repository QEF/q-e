  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
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
    SUBROUTINE selfen_ph_write()
    !----------------------------------------------------------------------------
    !!
    !! SP: Added lambda and phonon lifetime writing to file.
    !!
    USE kinds,         ONLY : DP
    USE global_var,    ONLY : gtemp, nqtotf, lambda_all, wf, gamma_all
    USE ep_constants,  ONLY : ryd2ev, kelvin2eV, ryd2mev
    USE input,         ONLY : nstemp
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE io_var,        ONLY : lambda_phself, linewidth_phself
    USE modes,         ONLY : nmodes
    !
    IMPLICIT NONE
    !
    ! Local variable
    CHARACTER(LEN = 20) :: tp
    !! string for temperature
    CHARACTER(LEN = 256) :: filephselfen
    !! file name of phonon selfenergy
    CHARACTER(LEN = 30)  :: myfmt
    !! Variable used for formatting output
    CHARACTER(LEN = 256) :: filephlinewid
    !! file name of phonon linewidth

    INTEGER :: itempphen
    !! Temperature counter for writing phonon selfen
    INTEGER :: iqq
    !! Counter on coarse q-point grid
    INTEGER :: imode
    !! Counter on mode
    !
    IF (mpime == ionode_id) THEN
      !
      DO itempphen = 1, nstemp
        WRITE(tp, "(f8.3)") gtemp(itempphen) * ryd2ev / kelvin2eV
        filephselfen = 'lambda.phself.' // trim(adjustl(tp)) // 'K'
        OPEN(UNIT = lambda_phself, FILE = filephselfen)
        WRITE(lambda_phself, '(/2x,a/)') '#Lambda phonon self-energy'
        WRITE(lambda_phself, *) '#Modes     ',(imode, imode = 1, nmodes)
        DO iqq = 1, nqtotf
          !
          myfmt = "(1000(3x,E15.5))"
          WRITE(lambda_phself,'(i9,4x)', ADVANCE = 'no') iqq
          WRITE(lambda_phself, FMT = myfmt) (REAL(lambda_all(imode, iqq, 1, itempphen)), imode = 1, nmodes)
          !
        ENDDO
        CLOSE(lambda_phself)
        !
        ! SP - 03/2019
        ! \Gamma = 1/\tau = phonon lifetime
        ! \Gamma = - 2 * Im \Pi^R where \Pi^R is the retarted phonon self-energy.
        ! Im \Pi^R = pi*k-point weight*[f(E_k+q) - f(E_k)]*delta[E_k+q - E_k - w_q]
        ! Since gamma_all = pi*k-point weight*[f(E_k) - f(E_k+q)]*delta[E_k+q - E_k - w_q] we have
        ! \Gamma = 2 * gamma_all
        filephlinewid = 'linewidth.phself.' // trim(adjustl(tp)) // 'K'
        OPEN(UNIT = linewidth_phself, FILE = filephlinewid)
        WRITE(linewidth_phself, '(a)') '# Phonon frequency and phonon lifetime in meV '
        WRITE(linewidth_phself, '(a)') '# Q-point  Mode   Phonon freq (meV)   Phonon linewidth (meV)'
        DO iqq = 1, nqtotf
          DO imode = 1, nmodes
            WRITE(linewidth_phself, '(i9,i6,E20.8,E22.10)') iqq, imode, &
                                   ryd2mev * wf(imode, iqq), 2.0d0 * ryd2mev * REAL(gamma_all(imode, iqq, 1, itempphen))
          ENDDO
        ENDDO
        CLOSE(linewidth_phself)
      ENDDO ! itempphen
    ENDIF ! mpime
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE selfen_ph_write
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE selfen_el_write(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all)
    !----------------------------------------------------------------------------
    !!
    !! Write self-energy
    !!
    USE kinds,         ONLY : DP
    USE global_var,    ONLY : lower_bnd, upper_bnd, nbndfst
    USE io_var,        ONLY : iufilsigma_all
    USE io_files,      ONLY : diropn
    USE ep_constants,  ONLY : zero
    USE mp,            ONLY : mp_barrier
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : meta_ionode, meta_ionode_id
    USE input,         ONLY : nstemp
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
    IF (meta_ionode) THEN
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
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE global_var,    ONLY : lower_bnd, upper_bnd, nbndfst
    USE io_var,        ONLY : iufilsigma_all
    USE io_files,      ONLY : prefix, tmp_dir, diropn
    USE ep_constants,  ONLY :  zero
    USE mp,            ONLY : mp_barrier, mp_bcast
    USE mp_world,      ONLY : mpime, world_comm
    USE io_global,     ONLY : meta_ionode, meta_ionode_id
    USE input,         ONLY : nstemp
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
    IF (meta_ionode) THEN
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
    CALL mp_bcast(exst, meta_ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iqq, meta_ionode_id, world_comm)
      CALL mp_bcast(sigmar_all, meta_ionode_id, world_comm)
      CALL mp_bcast(sigmai_all, meta_ionode_id, world_comm)
      CALL mp_bcast(zi_all, meta_ionode_id, world_comm)
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
    SUBROUTINE selfen_el_write_wfpt(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all, sigmar_dw_all)
    !----------------------------------------------------------------------------
    !!
    !! Write self-energy for WFPT. The only difference is that we write sigmar_dw_all.
    !! TODO (JML) : merge with selfen_el_write
    !!
    USE kinds,         ONLY : DP
    USE global_var,    ONLY : lower_bnd, upper_bnd, nbndfst
    USE io_var,        ONLY : iufilsigma_all
    USE io_files,      ONLY : diropn
    USE ep_constants,  ONLY : zero
    USE mp,            ONLY : mp_barrier
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE input,         ONLY : nstemp
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
    REAL(KIND = DP), INTENT(inout) :: sigmar_dw_all(nbndfst, nktotf, nstemp)
    !! Debye-Waller electron-phonon self-energy accross all pools
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
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            aux(i) = sigmar_dw_all(ibnd, ik, itemp)
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
      sigmar_dw_all(:, 1:lower_bnd - 1, :) = zero
    ENDIF
    IF (upper_bnd < nktotf) THEN
      sigmar_all(:, upper_bnd + 1:nktotf, :) = zero
      sigmai_all(:, upper_bnd + 1:nktotf, :) = zero
      zi_all(:, upper_bnd + 1:nktotf, :) = zero
      sigmar_dw_all(:, upper_bnd + 1:nktotf, :) = zero
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE selfen_el_write_wfpt
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE selfen_el_read_wfpt(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all, sigmar_dw_all)
    !----------------------------------------------------------------------------
    !!
    !! Self-energy reading for WFPT.  The only difference is that we read sigmar_dw_all.
    !! TODO (JML) : merge with selfen_el_read
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE global_var,    ONLY : lower_bnd, upper_bnd, nbndfst
    USE io_var,        ONLY : iufilsigma_all
    USE io_files,      ONLY : prefix, tmp_dir, diropn
    USE ep_constants,  ONLY :  zero
    USE mp,            ONLY : mp_barrier, mp_bcast
    USE mp_world,      ONLY : mpime, world_comm
    USE io_global,     ONLY : ionode_id
    USE input,         ONLY : nstemp
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
    REAL(KIND = DP), INTENT(out) :: sigmar_dw_all(nbndfst, nktotf, nstemp)
    !! Debyw-Waller electron-phonon self-energy accross all pools
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
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              sigmar_dw_all(ibnd, ik, itemp) = aux(i)
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
      CALL mp_bcast(sigmar_dw_all, ionode_id, world_comm)
      !
      ! Make everythin 0 except the range of k-points we are working on
      IF (lower_bnd > 1) THEN
        sigmar_all(:, 1:lower_bnd - 1, :) = zero
        sigmai_all(:, 1:lower_bnd - 1, :) = zero
        zi_all(:, 1:lower_bnd - 1, :) = zero
        sigmar_dw_all(:, 1:lower_bnd - 1, :) = zero
      ENDIF
      IF (upper_bnd < nktotf) THEN
        sigmar_all(:, upper_bnd + 1:nktotf, :) = zero
        sigmai_all(:, upper_bnd + 1:nktotf, :) = zero
        zi_all(:, upper_bnd + 1:nktotf, :) = zero
        sigmar_dw_all(:, upper_bnd + 1:nktotf, :) = zero
      ENDIF
      !
      WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from: ', iqq,'/', totq
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE selfen_el_read_wfpt
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE spectral_write(iqq, totq, nktotf, esigmar_all, esigmai_all)
    !----------------------------------------------------------------------------
    !!
    !! Write self-energy
    !!
    USE kinds,     ONLY : DP
    USE global_var,ONLY : lower_bnd, upper_bnd, nbndfst
    USE input,     ONLY : nstemp, wmin_specfun, wmax_specfun, nw_specfun
    USE io_var,    ONLY : iufilesigma_all
    USE io_files,  ONLY : diropn
    USE ep_constants,      ONLY : zero
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : mpime
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : my_pool_id
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
    IF (my_pool_id == ionode_id) THEN
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
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE global_var,    ONLY : lower_bnd, upper_bnd, nbndfst
    USE input,         ONLY : nstemp, nw_specfun
    USE io_var,        ONLY : iufilesigma_all
    USE io_files,      ONLY : prefix, tmp_dir, diropn
    USE ep_constants,  ONLY : zero
    USE mp,            ONLY : mp_barrier, mp_bcast
    USE mp_world,      ONLY : world_comm
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : my_pool_id
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
    REAL(KIND = DP) :: aux(2 * nbndfst * nktotf * nw_specfun * nstemp + 2)
    !! Vector to store the array
    !
    CHARACTER(LEN = 256) :: name1
    !
    !
    IF (my_pool_id == ionode_id) THEN
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

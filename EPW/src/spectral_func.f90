  !
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE spectral_func
  !----------------------------------------------------------------------
  !!
  !! This module contains the various spectral function routines
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE spectral_func_el_q(iqq, iq, totq, first_cycle)
    !-----------------------------------------------------------------------
    !!
    !!  Compute the electron spectral function including the  electron-
    !!  phonon interaction in the Migdal approximation.
    !!
    !!  We take the trace of the spectral function to simulate the photoemission
    !!  intensity. I do not consider the c-axis average for the time being.
    !!  The main approximation is constant dipole matrix element and diagonal
    !!  selfenergy. The diagonality can be checked numerically.
    !!
    !!  Use matrix elements, electronic eigenvalues and phonon frequencies
    !!  from ep-wannier interpolation
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iospectral_sup, iospectral
    USE modes,         ONLY : nmodes
    USE epwcom,        ONLY : nbndsub, eps_acustic, fsthick, ngaussw, &
                              degaussw, wmin_specfun, wmax_specfun, nw_specfun, &
                              shortrange, efermi_read, fermi_energy, restart, &
                              restart_step, nstemp
    USE pwcom,         ONLY : ef
    USE elph2,         ONLY : gtemp, etf, ibndmin, ibndmax, nkqf, xqf, nktotf, efnew, &
                              epf17, wkf, nkf, wf, wqf, xkf, nkqtotf, adapt_smearing, &
                              esigmar_all, esigmai_all, a_all, nbndfst, lower_bnd
    USE constants_epw, ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, ci, eps8
    USE constants,     ONLY : pi
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE io_selfen,     ONLY : spectral_write
    USE poolgathering, ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: first_cycle
    !! Use to determine weather this is the first cycle after restart
    INTEGER, INTENT(in) :: iqq
    !! Current q-point index in selecq
    INTEGER, INTENT(in) :: iq
    !! Current q-point index
    INTEGER, INTENT(in) :: totq
    !! Total number of q-point in window
    !
    ! Local variables
    CHARACTER(LEN = 20) :: tp
    !! String for temperatures
    CHARACTER(LEN = 256) :: filespec
    !! File name for spectral function
    CHARACTER(LEN = 256) :: filespecsup
    !! File name for supporting information
    INTEGER :: iw
    !! Counter on the frequency
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared in Ry^2
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
    REAL(KIND = DP) :: fact1
    !! Temporary variable to store $f_{mk+q}(T) + n_{q\nu}(T)$
    REAL(KIND = DP) :: fact2
    !! Temporary variable to store $1 - f_{mk+q}(T) + n_{q\nu}(T)$
    REAL(KIND = DP) :: weight
    !! Self-energy factor
    !!$$ N_q \Re( \frac{f_{mk+q}(T) + n_{q\nu}(T)}{ \varepsilon_{nk} - \varepsilon_{mk+q} + \omega_{q\nu} - i\delta }) $$
    !!$$ + N_q \Re( \frac{1- f_{mk+q}(T) + n_{q\nu}(T)}{ \varepsilon_{nk} - \varepsilon_{mk+q} - \omega_{q\nu} - i\delta }) $$
    REAL(KIND = DP) :: inv_eptemp
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse of the smearing for efficiency reasons
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP) :: specfun_sum
    !! Sum of spectral function
    REAL(KIND = DP) :: esigmar0
    !! static SE
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP) :: wq(nmodes)
    !! Phonon frequency on the fine grid
    REAL(KIND = DP) :: inv_wq(nmodes)
    !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
    REAL(KIND = DP) :: wgq(nmodes, nstemp)
    !! Bose occupation factor $n_{q\nu}(T)$
    REAL(KIND = DP) :: g2_tmp(nmodes)
    !! If the phonon frequency is too small discart g
    REAL(KIND = DP) :: fermi(nw_specfun)
    !! Spectral function
    REAL(KIND = DP) :: ww(nw_specfun)
    !! Current frequency
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Collect eigenenergies from all pools in parallel case
    !
    COMPLEX(KIND = DP) :: etmp1
    !! Temporary variable to store etmp1 = ekq - wq + ci * degaussw
    COMPLEX(KIND = DP) :: etmp2
    !! Temporary variable to strore etmp2 = ekq + wq + ci * degaussw
    COMPLEX(KIND = DP) :: etmpw1
    !! Temporary variable to store etmpw1 = ww - etmp1
    COMPLEX(KIND = DP) :: etmpw2
    !! Temporary variable to store etmpw1 = ww - etmp2
    COMPLEX(KIND = DP) :: fact
    !! SE factor
    !
    IF (adapt_smearing) CALL errore('spectral_func_el_q', 'adapt_smearing cannot be used with spectral functions ', 1)
    !
    DO itemp = 1, nstemp !loop over temperatures
      ! SP: Define the inverse so that we can efficiently multiply instead of dividing
      inv_eptemp   = one / gtemp(itemp)
      inv_degaussw = one / degaussw
      !
      ! energy range and spacing for spectral function
      !
      dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
      DO iw = 1, nw_specfun
        ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
      ENDDO
      !
      ! Now pre-treat phonon modes for efficiency
      ! Treat phonon frequency and Bose occupation
      wq(:) = zero
      DO imode = 1, nmodes
        wq(imode) = wf(imode, iq)
        IF (wq(imode) > eps_acustic) THEN
          g2_tmp(imode) = one
          wgq(imode, itemp)    = wgauss(-wq(imode) * inv_eptemp, -99)
          wgq(imode, itemp)    = wgq(imode, itemp) / (one - two * wgq(imode, itemp))
          inv_wq(imode) = one / (two * wq(imode))
        ELSE
          g2_tmp(imode) = zero
          wgq(imode, itemp)    = zero
          inv_wq(imode) = zero
        ENDIF
      ENDDO
      !
      IF (iqq == 1) THEN
        !
        WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
        WRITE(stdout, '(5x, "Electron Spectral Function in the Migdal Approximation")')
        WRITE(stdout, '(5x, a/)') REPEAT('=', 67)
        !
        IF (fsthick < 1.d3) WRITE(stdout, '(/5x, a ,f10.6, a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
        WRITE(stdout, '(/5x, a, f10.6, a)' ) 'Golden Rule strictly enforced with T = ', gtemp(itemp) * ryd2ev, ' eV'
        !
      ENDIF
      !
      ! Fermi level and corresponding DOS
      !
      IF (efermi_read) THEN
        ef0 = fermi_energy
      ELSE
        ef0 = efnew
      ENDIF
      !
      IF (iqq == 1) THEN
        WRITE(stdout, 100) degaussw * ryd2ev, ngaussw
        WRITE(stdout, '(a)') ' '
      ENDIF
      !
      ! SP: Sum rule added to conserve the number of electron.
      IF (iqq == 1) THEN
        WRITE(stdout, '(5x, a)') 'The sum rule to conserve the number of electron is enforced.'
        WRITE(stdout, '(5x, a)') 'The self energy is rescaled so that its real part is zero at the Fermi level.'
        WRITE(stdout, '(5x, a)') 'The sum rule replace the explicit calculation of the Debye-Waller term.'
        WRITE(stdout, '(a)') ' '
      ENDIF
      !
      IF (restart) THEN
        ! Make everythin 0 except the range of k-points we are working on
        esigmar_all(:, 1:lower_bnd - 1, :, :) = zero
        esigmar_all(:, lower_bnd + nkf:nktotf, :, :) = zero
        esigmai_all(:, 1:lower_bnd - 1, :, :) = zero
        esigmai_all(:, lower_bnd + nkf:nktotf, :, :) = zero
        !
      ENDIF
      !
      ! In the case of a restart do not add the first step
      IF (first_cycle .and. itemp == nstemp) THEN
        first_cycle = .FALSE.
      ELSE
        !
        ! loop over all k points of the fine mesh
        !
        fermicount = 0
        DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
          ! (but in this case they are the same)
          !
          IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
              (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
            !
            fermicount = fermicount + 1
            DO imode = 1, nmodes
              DO ibnd = 1, nbndfst
                !
                ! the energy of the electron at k (relative to Ef)
                ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
                !
                DO jbnd = 1, nbndfst
                  !
                  ! the fermi occupation for k+q
                  ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
                  wgkq = wgauss(-ekq * inv_eptemp, -99)
                  !
                  ! here we take into account the zero-point DSQRT(hbar/2M\omega)
                  ! with hbar = 1 and M already contained in the eigenmodes
                  ! g2 is Ry^2, wkf must already account for the spin factor
                  !
                  IF (shortrange .AND. (ABS(xqf(1, iq))> eps8 .OR. ABS(xqf(2, iq))> eps8 &
                     .OR. ABS(xqf(3, iq))> eps8)) THEN
                  !  ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary
                  !  !     number, in which case its square will be a negative number.
                    g2 = REAL((epf17(jbnd, ibnd, imode, ik)**two) * inv_wq(imode) * g2_tmp(imode))
                  ELSE
                    g2 = (ABS(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq(imode) * g2_tmp(imode)
                  ENDIF
                  !
                  fact1 =       wgkq + wgq(imode, itemp)
                  fact2 = one - wgkq + wgq(imode, itemp)
                  etmp1 = ekq - wq(imode) + ci * degaussw
                  etmp2 = ekq + wq(imode) + ci * degaussw
                  !
                  DO iw = 1, nw_specfun
                    !
                    etmpw1 = ww(iw) - etmp1
                    etmpw2 = ww(iw) - etmp2
                    !
                    fact = (fact1 / etmpw1) + (fact2 / etmpw2)
                    !
                    weight = wqf(iq) * REAL(fact)
                    !
                    ! \Re\Sigma [Eq. 3 in Comput. Phys. Commun. 209, 116 (2016)]
                    esigmar_all(ibnd, ik + lower_bnd - 1, iw, itemp) = esigmar_all(ibnd, ik + lower_bnd - 1, iw, itemp) + &
                                                                       g2 * weight
                    !
                    ! SP : Application of the sum rule
                    esigmar0 = - g2 *  wqf(iq) * REAL((fact1 / etmp1) + (fact2 / etmp2))
                    esigmar_all(ibnd, ik + lower_bnd - 1, iw, itemp) = esigmar_all(ibnd, ik + lower_bnd - 1, iw, itemp) - esigmar0
                    !
                    weight = wqf(iq) * AIMAG(fact)
                    !
                    ! \Im\Sigma [Eq. 3 in Comput. Phys. Commun. 209, 116 (2016)]
                    esigmai_all(ibnd, ik + lower_bnd - 1, iw, itemp) = esigmai_all(ibnd, ik + lower_bnd - 1, iw, itemp) + &
                                                                       g2 * weight
                    !
                  ENDDO
                ENDDO !jbnd
              ENDDO !ibnd
            ENDDO !imode
          ENDIF ! endif  fsthick
        ENDDO ! end loop on k
        !
        ! Creation of a restart point
        IF (restart) THEN
          IF (MOD(iqq, restart_step) == 0 .and. itemp == nstemp) THEN
            WRITE(stdout, '(5x, a, i10)' ) 'Creation of a restart point at ', iqq
            CALL mp_sum(esigmar_all, inter_pool_comm)
            CALL mp_sum(esigmai_all, inter_pool_comm)
            CALL mp_sum(fermicount, inter_pool_comm)
            CALL mp_barrier(inter_pool_comm)
            CALL spectral_write(iqq, totq, nktotf, esigmar_all, esigmai_all)
          ENDIF
        ENDIF
      ENDIF ! in case of restart, do not do the first one
      !
      ! The k points are distributed among pools: here we collect them
      !
    ENDDO ! first itemp
    IF (iqq == totq) THEN
        !
      ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_el_q', 'Error allocating xkf_all', 1)
      ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_el_q', 'Error allocating etf_all', 1)
      xkf_all(:, :) = zero
      etf_all(:, :) = zero
      !
#if defined(__MPI)
      !
      ! note that poolgather2 works with the doubled grid (k and k+q)
      !
      CALL poolgather2(3,       nkqtotf, nkqf, xkf, xkf_all)
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
      CALL mp_sum(esigmar_all, inter_pool_comm)
      CALL mp_sum(esigmai_all, inter_pool_comm)
      CALL mp_sum(fermicount, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
#else
      !
      xkf_all = xkf
      etf_all = etf
      !
#endif
      DO itemp = 1, nstemp ! second temperature loop to write data
        inv_eptemp = one / gtemp(itemp)
        !
        ! Output electron spectral function here after looping over all q-points
        ! (with their contributions summed in a etc.)
        !
        WRITE(stdout, '(5x, "WARNING: only the eigenstates within the Fermi window are meaningful")')
        !
        ! construct the trace of the spectral function (assume diagonal selfenergy
        ! and constant matrix elements for dipole transitions)
        !
        IF (mpime == ionode_id) THEN
          WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
          filespec = 'specfun.elself.' // trim(adjustl(tp)) // 'K'
          filespecsup = 'specfun_sup.elself.' // trim(adjustl(tp)) // 'K'
          OPEN(UNIT = iospectral, FILE = filespec )
          OPEN(UNIT = iospectral_sup, FILE = filespecsup )
          WRITE(iospectral, '(/2x, a/)') '#Electronic spectral function (meV)'
          WRITE(iospectral_sup, '(/2x, a/)') '#KS eigenenergies + real and im part of electronic self-energy (meV)'
          WRITE(iospectral, '(/2x, a/)') '#K-point    Energy[meV]     A(k,w)[meV^-1]'
          WRITE(iospectral_sup, '(/2x, a/)') '#K-point    Band   e_nk[eV]   w[eV]      Real Sigma[meV]  Im Sigma[meV]'
        ENDIF
        !
        DO ik = 1, nktotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          WRITE(stdout, '(/5x, "ik = ", i5, " coord.: ", 3f12.7, " Temp. : ", f8.3)') ik, xkf_all(:, ikk), &
                                                                                      gtemp(itemp) * ryd2ev / kelvin2eV
          WRITE(stdout, '(5x, a)') REPEAT('-', 67)
          !
          DO iw = 1, nw_specfun
            !
            DO ibnd = 1, nbndfst
              !
              !  the energy of the electron at k
              ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
              !
              a_all(iw, ik, itemp) = a_all(iw, ik, itemp) + ABS(esigmai_all(ibnd, ik, iw, itemp)) / pi / &
                   ((ww(iw) - ekk - esigmar_all(ibnd, ik, iw, itemp))**two + (esigmai_all(ibnd, ik, iw, itemp))**two)
              !
            ENDDO
            !
            WRITE(stdout, 101) ik, ryd2ev * ww(iw), a_all(iw, ik, itemp) / ryd2mev
            !
          ENDDO
          !
          WRITE(stdout, '(5x, a/)') REPEAT('-', 67)
          !
        ENDDO
        !
        DO ik = 1, nktotf
          !
          ! The spectral function should integrate to 1 for each k-point
          specfun_sum = 0.0
          !
          DO iw = 1, nw_specfun
            !
            fermi(iw) = wgauss(-ww(iw) * inv_eptemp, -99)
            !
            specfun_sum = specfun_sum + a_all(iw, ik, itemp) * fermi(iw) * dw
            !
            IF (mpime == ionode_id) WRITE(iospectral, '(2x, i7, 2x, f10.5, 2x, E12.5)') ik, ryd2ev * ww(iw), &
                                                                                  a_all(iw, ik, itemp) / ryd2mev
            !
          ENDDO
          !
          IF (mpime == ionode_id) WRITE(iospectral, '(a)') ' '
          IF (mpime == ionode_id) WRITE(iospectral, '(2x, a, 2x, E12.5)') '# Integrated spectral function ', specfun_sum
          !
        ENDDO
        !
        IF (mpime == ionode_id) CLOSE(iospectral)
        !
        DO ibnd = 1, nbndfst
          DO ik = 1, nktotf
            !
            ikk = 2 * ik - 1
            ikq = ikk + 1
            !
            !  the energy of the electron at k
            ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
            !
            DO iw = 1, nw_specfun
              !
              WRITE(stdout, 102) ik, ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2ev * ww(iw), &
                    ryd2mev * esigmar_all(ibnd, ik, iw, itemp), ryd2mev * esigmai_all(ibnd, ik, iw, itemp)
              !
              IF (mpime == ionode_id) &
              WRITE(iospectral_sup, 102) ik, ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2ev * ww(iw), &
                    ryd2mev * esigmar_all(ibnd, ik, iw, itemp), ryd2mev * esigmai_all(ibnd, ik, iw, itemp)
              !
            ENDDO
            !
          ENDDO
          !
          WRITE(stdout, *) ' '
          !
        ENDDO
        !
        IF (mpime == ionode_id) CLOSE(iospectral_sup)
        !
      ENDDO ! itemp
      DEALLOCATE(xkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_el_q', 'Error deallocating xkf_all', 1)
      DEALLOCATE(etf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_el_q', 'Error deallocating etf_all', 1)
      !
    ENDIF
    !
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    101 FORMAT(5x, 'ik = ', i7, '  w = ', f9.4, ' eV   A(k,w) = ', e12.5, ' meV^-1')
    102 FORMAT(2i9, 2x, f12.4, 2x, f12.4, 2x, f12.4, 2x, f12.4, 2x, f12.4)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE spectral_func_el_q
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE spectral_func_ph_q(iqq, iq, totq)
    !-----------------------------------------------------------------------
    !!
    !! Compute the imaginary part of the phonon self energy due to electron-
    !! phonon interaction in the Migdal approximation. This corresponds to
    !! the phonon linewidth (half width). The phonon frequency is taken into
    !! account in the energy selection rule.
    !!
    !! Use matrix elements, electronic eigenvalues and phonon frequencies
    !! from ep-wannier interpolation.  This routine is similar to the one above
    !! but it is ONLY called from within ephwann_shuffle and calculates
    !! the selfenergy for one phonon at a time.  Much smaller footprint on the disk
    !!
    !-----------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout, ionode_id
    USE io_var,    ONLY : iospectral_sup, iospectral
    USE modes,     ONLY : nmodes
    USE epwcom,    ONLY : nbndsub, fsthick, shortrange, ngaussw, degaussw, &
                          nsmear, delta_smear, eps_acustic, efermi_read, fermi_energy, &
                           nstemp, wmin_specfun, wmax_specfun, nw_specfun
    USE pwcom,     ONLY : nelec, ef
    USE klist_epw, ONLY : isk_dummy
    USE elph2,     ONLY : gtemp, epf17, ibndmax, ibndmin, etf, nbndfst, &
                          wkf, xqf, nkqf, nkf, wf, a_all_ph, efnew
    USE constants_epw, ONLY : kelvin2eV, ryd2mev, ryd2ev, one, two, zero, cone, ci, eps8
    USE constants,     ONLY : pi
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 20) :: tp
    !! String for temperatures
    CHARACTER(LEN = 256) :: filespec
    !! File name for spectral function
    CHARACTER(LEN = 256) :: filespecsup
    !! File name for supporting information
    INTEGER, INTENT(in) :: iqq
    !! Current q-point index from selecq
    INTEGER, INTENT(in) :: iq
    !! Current q-point index
    INTEGER, INTENT(in) :: totq
    !! Total q-points in selecq window
    !
    ! Local variables
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: iw
    !! Counter on frequency for the phonon spectra
    INTEGER :: itemp
    !! Counter on temperature
    INTEGER :: iqq_write
    !! Counter on q-point to write files
    !
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared in Ry^2
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: dosef
    !! Density of state N(Ef)
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: wgkk
    !! Fermi-Dirac occupation factor $f_{nk}(T)$
    REAL(KIND = DP) :: wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
    REAL(KIND = DP) :: fact
    !! Temporary variable to store fact = wgkq - wgkk
    REAL(KIND = DP) :: etmp1
    !! Temporary variable to store etmp1 = ekq - ekk
    REAL(KIND = DP) :: etmp2
    !! Temporary variable to store etmp2 = ekq - ekk - ww
    REAL(KIND = DP) :: weight
    !! Self-energy factor
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse Gaussian for efficiency reasons
    REAL(KIND = DP) :: inv_eptemp
    !! Inverse temperature
    REAL(KIND = DP) :: inv_pi
    !! Inverse pi
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! Function to compute the Density of States at the Fermi level
    REAL(KIND = DP), EXTERNAL :: efermig
    !! Function to compute the Fermi energy
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    REAL(KIND = DP) :: gamma0(nmodes)
   !! Phonon self-energy
    REAL(KIND = DP) :: wq(nmodes)
    !! Phonon frequency on the fine grid
    REAL(KIND = DP) :: inv_wq(nmodes)
    !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
    REAL(KIND = DP) :: dwq(nmodes)
    !! $2\omega_{q\nu}$ defined for efficiency reasons
    REAL(KIND = DP) :: g2_tmp(nmodes)
    !! If the phonon frequency is too small discart g
    REAL(KIND = DP) :: ww(nw_specfun)
    !! Current frequency
    REAL(KIND = DP) :: gammai_all(nw_specfun, nmodes)
    !! Imaginary part of the frequency dependent spectral function
    REAL(KIND = DP) :: gammar_all(nw_specfun, nmodes)
    !!  Real part of the Phonon self-energy (freq. dependent for spectral function)
    !
    dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
    DO iw = 1, nw_specfun
      ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
    ENDDO
    !
    ! Now pre-treat phonon modes for efficiency
    ! Treat phonon frequency and Bose occupation
    wq(:)    = zero
    dwq(:)   = zero
    DO imode = 1, nmodes
      wq(imode) = wf(imode, iq)
      dwq(imode) = two * wq(imode)
      IF (wq(imode) > eps_acustic) THEN
        g2_tmp(imode) = one
        inv_wq(imode) = one / (two * wq(imode))
      ELSE
        g2_tmp(imode) = zero
        inv_wq(imode) = zero
      ENDIF
    ENDDO
    !
    DO itemp = 1, nstemp
      gammar_all(:, :) = zero
      gammai_all(:, :) = zero
      !
      ! Thomas-Fermi screening according to Resta PRB 1977
      ! Here specific case of Diamond
      !eps0   = 5.7
      !rtf    = 2.76
      !qtf    = 1.36
      !qsquared = (xqf(1,iq)**2 + xqf(2,iq)**2 + xqf(3,iq)**2) * tpiba2
      !epstf =  (qtf**2 + qsquared) / (qtf**2/eps0 * sin (sqrt(qsquared)*rtf)/(sqrt(qsquared)*rtf)+qsquared)
      !
      IF (iqq == 1) THEN
        WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
        WRITE(stdout, '(5x, "Phonon Spectral Function Self-Energy in the Migdal Approximation (on the fly)")')
        WRITE(stdout, '(5x, a/)') REPEAT('=', 67)
        !
        IF (fsthick < 1.d3) WRITE(stdout, '(/5x, a, f10.6, a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
        WRITE(stdout, '(/5x, a, f10.6, a)' ) 'Golden Rule strictly enforced with T = ', gtemp(itemp) * ryd2ev, ' eV'
      ENDIF
      !
      ! SP: Multiplication is faster than division ==> Important if called a lot in inner loops
      inv_degaussw = one / degaussw
      inv_eptemp   = one / gtemp(itemp)
      inv_pi       = one / pi
      !
      ! Fermi level and corresponding DOS
      !
      IF (efermi_read) THEN
        !
        ef0 = fermi_energy
        !
      ELSEIF (nsmear > 1) THEN
        !
        ! if some bands are skipped (nbndskip /= 0), nelec has already been recalculated in ephwann_shuffle
        ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk_dummy)
        !
      ELSE !SP: This is added for efficiency reason because the efermig routine is slow
        ef0 = efnew
      ENDIF
      !
      ! N(Ef) in the equation for lambda is the DOS per spin
      dosef = dos_ef(ngaussw, degaussw, ef0, etf, wkf, nkqf, nbndsub)
      dosef = dosef / two
      !
      IF (iqq == 1) THEN
        WRITE(stdout, 100) degaussw * ryd2ev, ngaussw
        WRITE(stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
      ENDIF
      !
      CALL start_clock('PH SPECTRAL-FUNCTION')
      !
      fermicount = 0
      gamma0(:)  = zero
      !
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        ! Here we must have ef, not ef0, to be consistent with ephwann_shuffle
        IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
            (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
          !
          fermicount = fermicount + 1
          !
          DO imode = 1, nmodes
            !
            DO ibnd = 1, nbndfst
              !
              !  the fermi occupation for k
              ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
              wgkk = wgauss(-ekk * inv_eptemp, -99)
              !
              DO jbnd = 1, nbndfst
                !
                !  the fermi occupation for k+q
                ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
                wgkq = wgauss(-ekq * inv_eptemp, -99)
                !
                ! here we take into account the zero-point DSQRT(hbar/2M\omega)
                ! with hbar = 1 and M already contained in the eigenmodes
                ! g2 is Ry^2, wkf must already account for the spin factor
                !
                IF (shortrange .AND. (ABS(xqf(1, iq))> eps8 .OR. ABS(xqf(2, iq))> eps8 &
                   .OR. ABS(xqf(3, iq))> eps8 )) THEN
                  ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary
                  !     number, in which case its square will be a negative number.
                  g2 = REAL((epf17(jbnd, ibnd, imode, ik)**two) * inv_wq(imode) * g2_tmp(imode)) !* epsTF
                ELSE
                  g2 = (ABS(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq(imode) * g2_tmp(imode) !* epsTF
                ENDIF
                !
                ! SP - 03/2019 - Retarded phonon self-energy
                !                See Eq. 145 of RMP 89, 015003 (2017)
                ! \Pi^R = k-point weight * [ [f(E_k+q) - f(E_k)]/ [E_k+q - E_k -w_q - id]
                !                           -[f(E_k+q) - f(E_k)]/ [E_k+q - E_k - id] ]
                ! The second term is gamma0 (static)
                !
                fact = wgkq - wgkk
                etmp1 = ekq - ekk
                weight = wkf(ikk) * fact * REAL(cone / (etmp1 + ci * degaussw))
                !
                gamma0(imode) = gamma0(imode) + weight * g2
                !
                DO iw = 1, nw_specfun
                  !
                  etmp2 = etmp1 - ww(iw)
                  weight = wkf(ikk) * fact * REAL(cone / (etmp2 + ci * degaussw))
                  gammar_all(iw, imode) = gammar_all(iw, imode) + weight * g2
                  !
                  ! Normal implementation
                  !weight = wkf (ikk) * fact * AIMAG(cone / (etmp2 + ci * degaussw))
                  !
                  ! More stable:
                  ! Analytical im. part
                  weight = pi * wkf(ikk) * fact * w0gauss(etmp2 * inv_degaussw, 0) * inv_degaussw
                  !
                  gammai_all(iw, imode) = gammai_all(iw, imode) + weight * g2
                  !
                ENDDO
              ENDDO ! jbnd
            ENDDO   ! ibnd
          ENDDO ! loop on q-modes
        ENDIF ! endif fsthick
      ENDDO ! loop on k
      !
      CALL stop_clock('PH SPECTRAL-FUNCTION')
      !
#if defined(__MPI)
      !
      ! collect contributions from all pools (sum over k-points) this finishes the integral over the BZ  (k)
      !
      CALL mp_sum(gammai_all, inter_pool_comm)
      CALL mp_sum(gammar_all, inter_pool_comm)
      CALL mp_sum(gamma0, inter_pool_comm)
      CALL mp_sum(fermicount, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
#endif
      !
      WRITE(stdout, '(5x, a)')
      !
      IF (iqq == 1 .and. itemp == 1) THEN
        IF (mpime == ionode_id) THEN
          WRITE(tp, "(f8.1)") gtemp(itemp) * ryd2ev / kelvin2eV
          filespecsup = 'specfun_sup.phon'! // trim(adjustl(tp)) // 'K'
          OPEN(UNIT = iospectral_sup, FILE = filespecsup)
          WRITE(iospectral_sup, '(2x, a)') '#Phonon eigenenergies + real and im part of phonon self-energy (meV)'
          WRITE(iospectral_sup, '(2x, a)') '#Q-point    Mode      Temp.[K]       w_q[eV]        w[eV]    &
                                           &Real Sigma(w)[meV]   Real Sigma(w=0)[meV]     Im Sigma(w)[meV]'
        ENDIF
      ENDIF
      !
      ! Write to output file
      WRITE(stdout, '(/5x, a)') 'Real and Imaginary part of the phonon self-energy (omega=0) without gamma0.'
      DO imode = 1, nmodes
        ! Real and Im part of Phonon self-energy at 0 freq.
        WRITE(stdout, 105) imode, ryd2ev * wq(imode), ryd2mev * gammar_all(1, imode), ryd2mev * gammai_all(1, imode)
      ENDDO
      WRITE(stdout, '(5x, a, i8, a, i8)' ) 'Number of (k,k+q) pairs on the Fermi surface: ', fermicount, ' out of ', totq
      !
      ! Write to support files
      DO iw = 1, nw_specfun
        !
        DO imode = 1, nmodes
          !
          !a_all(iw,iq) = a_all(iw,iq) + ABS(gammai_all(imode,iq,iw) ) / pi / &
          !      ( ( ww - wq - gammar_all (imode,iq,iw) + gamma0 (imode))**two + (gammai_all(imode,iq,iw) )**two )
          ! SP: From Eq. 16 of PRB 9, 4733 (1974)
          !    Also in Eq.2 of PRL 119, 017001 (2017).
          a_all_ph(iw, iqq, itemp) = a_all_ph(iw, iqq, itemp) + inv_pi * dwq(imode)**two * ABS(gammai_all(iw, imode)) / &
                              ((ww(iw)**two - wq(imode)**two - dwq(imode) * (gammar_all(iw, imode) - gamma0(imode)))**two + &
                               (dwq(imode) * gammai_all(iw, imode))**two)
          !
          IF (mpime == ionode_id) THEN
            WRITE(iospectral_sup, 102) iq, imode, gtemp(itemp) * ryd2ev / kelvin2eV, ryd2ev * wq(imode), ryd2ev * ww(iw), &
                                       ryd2mev * gammar_all(iw, imode), ryd2mev * gamma0(imode), ryd2mev * gammai_all(iw, imode)
          ENDIF
          !
        ENDDO
        !
!        IF (mpime == ionode_id) THEN
!          WRITE(iospectral, 103) iq, ryd2ev * ww(iw), a_all_ph(iw, iqq) / ryd2mev ! print to file
!        ENDIF
        !
      ENDDO
    ENDDO ! itemp
      !
    IF (iqq == totq) THEN
      IF (mpime == ionode_id) THEN
        DO itemp = 1, nstemp
          WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
          filespec = 'specfun.phon.' // trim(adjustl(tp)) // 'K'
          OPEN(UNIT = iospectral, FILE = filespec)
          WRITE(iospectral, '(/2x, a)') '#Phonon spectral function (meV)'
          WRITE(iospectral, '(/2x, a)') '#Q-point    Energy[eV]     A(q,w)[meV^-1]'
          DO iqq_write = 1, totq
            DO iw = 1, nw_specfun
              WRITE(iospectral, 103) iqq_write, ryd2ev * ww(iw), a_all_ph(iw, iqq_write, itemp) / ryd2mev ! print to file
            ENDDO
          ENDDO
          CLOSE(iospectral)
        ENDDO
        CLOSE(iospectral_sup)
      ENDIF
    ENDIF
    WRITE(stdout, '(5x, a/)') REPEAT('-',67)
    !
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    101 FORMAT(5x, 'DOS =', f10.6,' states/spin/eV/Unit Cell at Ef=', f10.6, ' eV')
    102 FORMAT(2i9, 2x, f8.3, 2x, f12.5, 2x, f12.5, 2x, E22.14, 2x, E22.14, 2x, E22.14)
    103 FORMAT(2x, i7, 2x, f12.5, 2x, E22.14)
    105 FORMAT(5x, 'Omega( ', i3, ' )=', f9.4,' eV   Re[Pi]=', f15.6, ' meV Im[Pi]=', f15.6, ' meV')
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE spectral_func_ph_q
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE spectral_func_pl_q(iqq, iq, totq, first_cycle)
    !-----------------------------------------------------------------------
    !!
    !!  Compute the electron spectral function including the  electron-
    !!  phonon interaction in the Migdal approximation.
    !!
    !!  We take the trace of the spectral function to simulate the photoemission
    !!  intensity. I do not consider the c-axis average for the time being.
    !!  The main approximation is constant dipole matrix element and diagonal
    !!  selfenergy. The diagonality can be checked numerically.
    !!
    !!  Use matrix elements, electronic eigenvalues and phonon frequencies
    !!  from ep-wannier interpolation
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iospectral_sup, iospectral
    USE epwcom,        ONLY : nbndsub, fsthick, ngaussw, degaussw, nw_specfun, &
                              wmin_specfun, wmax_specfun, efermi_read, fermi_energy, &
                              nstemp, nel, meff, epsiheg, restart, restart_step
    USE pwcom,         ONLY : nelec, ef
    USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, nbndfst, wkf, nkf, wqf, xkf, &
                              nkqtotf, xqf, dmef, esigmar_all, esigmai_all, a_all, &
                              gtemp, nktotf, lower_bnd, efnew
    USE constants_epw, ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, ci, eps6
    USE constants,     ONLY : pi
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE cell_base,     ONLY : omega, alat, bg
    USE selfen,        ONLY : get_eps_mahan
    USE io_selfen,     ONLY : spectral_write
    USE poolgathering, ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: first_cycle
    !! Use to determine weather this is the first cycle after restart
    INTEGER, INTENT(in) :: iqq
    !! Q-point index in selecq
    INTEGER, INTENT(in) :: iq
    !! Q-point index
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points in fsthick window
    !
    ! Local variables
    CHARACTER(LEN = 20) :: tp
    !! String for temperatures
    CHARACTER(LEN = 256) :: filespec
    !! File name for spectral function
    CHARACTER(LEN = 256) :: filespecsup
    !! File name for supporting information
    INTEGER :: iw
    !! Counter on the frequency
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared in Ry^2
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: ekk
    !! Eigen energy at k on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekk1
    !! Eigen energy at k on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy at k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: wq
    !! Plasmon frequency
    REAL(KIND = DP) :: wgq
    !! Bose occupation factor $n_{q wpl}(T)$
    REAL(KIND = DP) :: wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
    REAL(KIND = DP) :: fact1
    !! Temporary variable to store $f_{mk+q}(T) + n_{q wpl}(T)$
    REAL(KIND = DP) :: fact2
    !! Temporary variable to store $1 - f_{mk+q}(T) + n_{q wpl}(T)$
    REAL(KIND = DP) :: weight
    !! SE factor
    REAL(KIND = DP) :: inv_eptemp
    !! Inverse of temperature defined for efficiency reasons
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse of degaussw defined for efficiency reasons
    REAL(KIND = DP) :: sq_degaussw
    !! Squared degaussw defined for efficiency reasons
    REAL(KIND = DP) :: g2_tmp
    !! Temporary variable defined for efficiency reasons
    REAL(KIND = DP) :: dw
    !! Spectral frequency increment
    REAL(KIND = DP) :: specfun_sum
    !! Sum of spectral function
    REAL(KIND = DP) :: esigmar0
    !! static SE
    REAL(KIND = DP) :: tpiba_new
    !! 2 \pi / alat
    REAL(KIND = DP) :: kf
    !! Fermi wave-vector
    REAL(KIND = DP) :: vf
    !! Fermi velocity
    REAL(KIND = DP) :: fermiheg
    !! Fermi energy of a homageneous electron gas
    REAL(KIND = DP) :: qnorm
    !! |q|
    REAL(KIND = DP) :: qin
    !! (2 \pi / alat) |q|
    REAL(KIND = DP) :: sq_qin
    !! Squared qin defined for efficiency reasons
    REAL(KIND = DP) :: wpl0
    !! Plasmon frequency
    REAL(KIND = DP) :: eps0
    !! Dielectric function at zero frequency
    REAL(KIND = DP) :: deltaeps
    !!
    REAL(KIND = DP) :: qcut
    !! Cut-off of the maximum wave-vector of plasmon modes (qcut = wpl0 / vf)
    REAL(KIND = DP) :: qtf
    !! Thomas-Fermi screening wave-vector
    REAL(KIND = DP) :: dipole
    !! Dipole
    REAL(KIND = DP) :: rs
    !! Spherical radius used to describe the density of an electron gas
    REAL(KIND = DP) :: degen
    !! Degeneracy of the electron gas
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP) :: q(3)
    !! The q-point in cartesian unit.
    REAL(KIND = DP) :: fermi(nw_specfun)
    !! Spectral function
    REAL(KIND = DP) :: ww(nw_specfun)
    !! Current frequency
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Collect eigenenergies from all pools in parallel case
    !
    COMPLEX(KIND = DP) :: etmp1
    !! Temporary variable to store etmp1 = ekq - wq + ci * degaussw
    COMPLEX(KIND = DP) :: etmp2
    !! Temporary variable to strore etmp2 = ekq + wq + ci * degaussw
    COMPLEX(KIND = DP) :: etmpw1
    !! Temporary variable to store etmpw1 = ww - etmp1
    COMPLEX(KIND = DP) :: etmpw2
    !! Temporary variable to store etmpw1 = ww - etmp2
    COMPLEX(KIND = DP) :: fact
    !! SE factor
    !
    ! loop over temperatures can be introduced
    !
    inv_degaussw = one / degaussw
    sq_degaussw = degaussw * degaussw
    DO itemp = 1, nstemp
      inv_eptemp   = one / gtemp(itemp)
      ! energy range and spacing for spectral function
      !
      dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
      DO iw = 1, nw_specfun
        ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
      ENDDO
      !
      IF (iqq == 1) THEN
        !
        WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
        WRITE(stdout, '(5x, "Electron Spectral Function in the Migdal Approximation")')
        WRITE(stdout, '(5x, a/)') REPEAT('=', 67)
        !
        IF (fsthick < 1.d3) WRITE(stdout, '(/5x, a, f10.6, a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
        WRITE(stdout, '(/5x, a, f10.6, a)' ) 'Golden Rule strictly enforced with T = ', gtemp(itemp) * ryd2ev, ' eV'
        !
      ENDIF
      !
      ! Fermi level
      !
      IF (efermi_read) THEN
        ef0 = fermi_energy
      ELSE
        ef0 = efnew
      ENDIF
      !
      IF (iqq == 1) THEN
        WRITE(stdout, 100) degaussw * ryd2ev, ngaussw
        WRITE(stdout, '(a)') ' '
      ENDIF
      !
      ! SP: Sum rule added to conserve the number of electron.
      IF (iqq == 1) THEN
        WRITE(stdout, '(5x, a)') 'The sum rule to conserve the number of electron is enforced.'
        WRITE(stdout, '(5x, a)') 'The self energy is rescaled so that its real part is zero at the Fermi level.'
        WRITE(stdout, '(5x, a)') 'The sum rule replace the explicit calculation of the Debye-Waller term.'
        WRITE(stdout, '(a)') ' '
      ENDIF
      !
      !nel      =  0.01    ! this should be read from input - # of doping electrons
      !epsiheg  =  12.d0   ! this should be read from input - # dielectric constant at zero doping
      !meff     =  0.25    ! this should be read from input - effective mass
      !
      tpiba_new = two * pi / alat
      degen     = one
      !
      ! Based on Eqs. (5.3)-(5.6) and (5.127) of Mahan 2000.
      !
      ! omega is the unit cell volume in Bohr^3
      rs = (3.d0 / (4.d0 * pi * nel / omega / degen))**(1.d0 / 3.d0) * meff * degen
      kf = (3.d0 * (pi**2.d0) * nel / omega / degen)**(1.d0 / 3.d0)
      vf = (1.d0 / meff) * kf
      !
      ! fermiheg in [Ry] (multiplication by 2 converts from Ha to Ry)
      fermiheg = 2.d0 * (1.d0 / (2.d0 * meff)) * kf**2.d0
      ! qtf in ! [a.u.]
      qtf = DSQRT(6.d0 * pi * nel / omega / degen / (fermiheg / 2.d0))
      ! wpl0 in [Ry] (multiplication by 2 converts from Ha to Ry)
      wpl0 = two * DSQRT(4.d0 * pi * nel / omega / meff / epsiheg)
      wq = wpl0
      !
      q(:) = xqf(:, iq)
      CALL cryst_to_cart(1, q, bg, 1)
      qnorm = DSQRT(q(1)**two + q(2)**two + q(3)**two)
      qin = qnorm * tpiba_new
      sq_qin = qin * qin
      !
      ! qcut in [Ha] (1/2 converts from Ry to Ha)
      qcut = wpl0 / vf / tpiba_new / 2.d0
      !
      !IF (.TRUE.) qcut = qcut / 2.d0 ! renormalize to account for Landau damping
      !
      ! qin should be in atomic units for Mahan formula
      CALL get_eps_mahan(qin, rs, kf, eps0)
      deltaeps = -(1.d0 / (epsiheg + eps0 - 1.d0) - 1.d0 / epsiheg)
      !
      g2_tmp = 4.d0 * pi * (wq * deltaeps / 2.d0) / omega * 2.d0
      !
      IF (iqq == 1) THEN
        WRITE(stdout, '(12x, " nel       = ", E15.6)') nel
        WRITE(stdout, '(12x, " meff      = ", E15.6)') meff
        WRITE(stdout, '(12x, " rs        = ", E15.6)') rs
        WRITE(stdout, '(12x, " kf        = ", E15.6)') kf
        WRITE(stdout, '(12x, " vf        = ", E15.6)') vf
        WRITE(stdout, '(12x, " fermi_en  = ", E15.6)') fermiheg
        WRITE(stdout, '(12x, " qtf       = ", E15.6)') qtf
        WRITE(stdout, '(12x, " wpl       = ", E15.6)') wpl0
        WRITE(stdout, '(12x, " qcut      = ", E15.6)') qcut
        WRITE(stdout, '(12x, " eps0      = ", E15.6)') eps0
        WRITE(stdout, '(12x, " epsiheg   = ", E15.6)') epsiheg
        WRITE(stdout, '(12x, " deltaeps  = ", E15.6)') deltaeps
      ENDIF
      !
      IF (restart) THEN
        ! Make everythin 0 except the range of k-points we are working on
        esigmar_all(:, 1:lower_bnd - 1, :, :) = zero
        esigmar_all(:, lower_bnd + nkf:nktotf, :, :) = zero
        esigmai_all(:, 1:lower_bnd - 1, :, :) = zero
        esigmai_all(:, lower_bnd + nkf:nktotf, :, :) = zero
        !
      ENDIF
      !
      ! In the case of a restart do not add the first step
      IF (first_cycle .and. itemp == nstemp) THEN
        first_cycle = .FALSE.
      ELSE
        IF (qnorm < qcut) THEN
          !
          ! wq is the plasmon frequency
          ! Bose occupation
          wgq = wgauss(-wq * inv_eptemp, -99)
          wgq = wgq / (one - two * wgq)
          !
          ! loop over all k points of the fine mesh
          !
          fermicount = 0
          DO ik = 1, nkf
            !
            ikk = 2 * ik - 1
            ikq = ikk + 1
            !
            ! here we must have ef, not ef0, to be consistent with ephwann_shuffle (but in this case they are the same)
            !
            IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
                (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
              !
              fermicount = fermicount + 1
              !
              DO ibnd = 1, nbndfst
                !
                !  the energy of the electron at k (relative to Ef)
                ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
                !
                DO jbnd = 1, nbndfst
                  !
                  ekk1 = etf(ibndmin - 1 + jbnd, ikk) - ef0
                  ! the energy of the electron at k+q (relative to Ef)
                  ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
                  ! the Fermi occupation at k+q
                  wgkq = wgauss(-ekq * inv_eptemp, -99)
                  !
                  ! Computation of the dipole
                  IF (ibnd == jbnd) THEN
                    IF (qnorm > eps6) THEN
                      dipole = one / sq_qin
                    ELSE
                      dipole = zero
                    ENDIF
                  ELSE
                    IF (ABS(ekq - ekk1) > eps6) THEN
                      dipole = REAL(      dmef(1, ibndmin - 1 + jbnd, ibndmin - 1 + ibnd, ikk) *  &
                                    CONJG(dmef(1, ibndmin - 1 + jbnd, ibndmin - 1 + ibnd, ikk)) / &
                                    ((ekk1 - ekk)**two + sq_degaussw))
                    ELSE
                      dipole = zero
                    ENDIF
                  ENDIF
                  !
                  ! The q^-2 is cancelled by the q->0 limit of the dipole.
                  ! See e.g., pg. 258 of Grosso Parravicini.
                  ! electron-plasmon scattering matrix elements squared
                  g2 = dipole * g2_tmp
                  !
                  fact1 =       wgkq + wgq
                  fact2 = one - wgkq + wgq
                  etmp1 = ekq - wq + ci * degaussw
                  etmp2 = ekq + wq + ci * degaussw
                  !
                  DO iw = 1, nw_specfun
                    !
                    etmpw1 = ww(iw) - etmp1
                    etmpw2 = ww(iw) - etmp2
                    !
                    fact = (fact1 / etmpw1) + (fact2 / etmpw2)
                    !
                    weight = wqf(iq) * REAL(fact)
                    !
                    ! \Re\Sigma [Eq. 3 in Comput. Phys. Commun. 209, 116 (2016)]
                    esigmar_all(ibnd, ik + lower_bnd - 1, iw, itemp) = esigmar_all(ibnd, ik + lower_bnd - 1, iw, itemp) + &
                                                                       g2 * weight
                    !
                    ! SP : Application of the sum rule
                    esigmar0 = - g2 *  wqf(iq) * REAL((fact1 / etmp1) + (fact2 / etmp2))
                    esigmar_all(ibnd, ik + lower_bnd - 1, iw, itemp) = esigmar_all(ibnd, ik + lower_bnd - 1, iw, itemp) - &
                                                                       esigmar0
                    !
                    weight = wqf(iq) * AIMAG(fact)
                    !
                    ! \Im\Sigma [Eq. 3 in Comput. Phys. Commun. 209, 116 (2016)]
                    esigmai_all(ibnd, ik + lower_bnd - 1, iw, itemp) = esigmai_all(ibnd, ik + lower_bnd - 1, iw, itemp) + &
                                                                       g2 * weight
                    !
!                    WRITE(stdout, '(5x, f8.5, f8.5)') g2_tmp, dipole
                  ENDDO
                ENDDO !jbnd
              ENDDO !ibnd
            ENDIF ! endif  fsthick
          ENDDO ! end loop on k
        ENDIF ! endif qnorm
        !
        ! Creation of a restart point
        IF (restart) THEN
          IF (MOD(iqq, restart_step) == 0 .and. itemp == nstemp) THEN
            WRITE(stdout, '(5x, a, i10)' ) 'Creation of a restart point at ', iqq
            CALL mp_sum(esigmar_all, inter_pool_comm)
            CALL mp_sum(esigmai_all, inter_pool_comm)
            CALL mp_sum(fermicount, inter_pool_comm)
            CALL mp_barrier(inter_pool_comm)
            CALL spectral_write(iqq, totq, nktotf, esigmar_all, esigmai_all)
          ENDIF
        ENDIF
      ENDIF ! in case of restart, do not do the first one
    ENDDO ! itemp
    !
    ! The k points are distributed among pools: here we collect them
    ! 
    IF (iqq == totq) THEN
      ! Collect pools and write the spectral function
      !
      ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_pl_q', 'Error allocating xkf_all', 1)
      ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_pl_q', 'Error allocating etf_all', 1)
      xkf_all(:, :) = zero
      etf_all(:, :) = zero
      !
#if defined(__MPI)
      !
      ! Note that poolgather2 works with the doubled grid (k and k+q)
      !
      CALL poolgather2(3, nkqtotf, nkqf, xkf, xkf_all)
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
      CALL mp_sum(esigmar_all, inter_pool_comm)
      CALL mp_sum(esigmai_all, inter_pool_comm)
      CALL mp_sum(fermicount, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
#else
      !
      xkf_all = xkf
      etf_all = etf
      !
#endif
      DO itemp = 1, nstemp
        inv_eptemp = one / gtemp(itemp)
        !
        ! Output electron spectral function here after looping over all q-points (with their contributions summed in a etc.)
        !
        WRITE(stdout, '(5x, "WARNING: only the eigenstates within the Fermi window are meaningful")')
        !
        ! construct the trace of the spectral function (assume diagonal selfenergy
        ! and constant matrix elements for dipole transitions)
        !
        IF (mpime == ionode_id) then
          WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
          filespec = 'specfun.plself.' // trim(adjustl(tp)) // 'K'
          filespecsup = 'specfun_sup.plself.' // trim(adjustl(tp)) // 'K'
          OPEN(UNIT = iospectral, FILE = filespec )
          OPEN(UNIT = iospectral_sup, FILE = filespecsup )
          WRITE(iospectral, '(/2x, a/)') '#Electron-plasmon spectral function (meV)'
          WRITE(iospectral_sup, '(/2x, a/)') '#KS eigenenergies + real and im part of electron-plasmon self-energy (meV)'
          WRITE(iospectral, '(/2x, a/)') '#K-point    Energy[meV]     A(k,w)[meV^-1]'
          WRITE(iospectral_sup, '(/2x, a/)') '#K-point    Band   e_nk[eV]   w[eV]       Real Sigma[meV]  Im Sigma[meV]'
        ENDIF
        !
        DO ik = 1, nktotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          WRITE(stdout, '(/5x, "ik = ", i5, " coord.: ", 3f12.7, " Temp.: ", f8.3 )') ik, xkf_all(:, ikk), &
                                                                                      gtemp(itemp) * ryd2ev / kelvin2eV
          WRITE(stdout, '(5x, a)') REPEAT('-', 67)
          !
          DO iw = 1, nw_specfun
            !
            DO ibnd = 1, nbndfst
              !
              !  the energy of the electron at k
              ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
              !
              a_all(iw, ik, itemp) = a_all(iw, ik, itemp) + ABS(esigmai_all(ibnd, ik, iw, itemp) ) / pi / &
                   ((ww(iw) - ekk - esigmar_all(ibnd, ik, iw, itemp))**two + (esigmai_all(ibnd, ik, iw, itemp))**two)
              !
            ENDDO
            !
            WRITE(stdout, 101) ik, ryd2ev * ww(iw), a_all(iw, ik, itemp) / ryd2mev
            !
          ENDDO
          !
          WRITE(stdout, '(5x, a/)') REPEAT('-', 67)
          !
        ENDDO
        !
        DO ik = 1, nktotf
          !
          ! The spectral function should integrate to 1 for each k-point
          specfun_sum = zero
          !
          DO iw = 1, nw_specfun
            !
            fermi(iw) = wgauss(-ww(iw) * inv_eptemp, -99)
            specfun_sum = specfun_sum + a_all(iw, ik, itemp) * fermi(iw) * dw !/ ryd2mev
            !
           IF (mpime == ionode_id) WRITE(iospectral, '(2x, i7, 2x, f10.5, 2x, E12.5)') ik, ryd2ev * ww(iw), &
                                                                                 a_all(iw, ik, itemp) / ryd2mev
          ENDDO
          !
          IF (mpime == ionode_id) WRITE(iospectral, '(a)') ' '
          IF (mpime == ionode_id) WRITE(iospectral, '(2x, a, 2x, E12.5)') '# Integrated spectral function ', specfun_sum
        ENDDO
        !
        IF (mpime == ionode_id) CLOSE(iospectral)
        !
        DO ibnd = 1, nbndfst
          !
          DO ik = 1, nktotf
            !
            ikk = 2 * ik - 1
            ikq = ikk + 1
            !
            !  the energy of the electron at k
            ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
            !
            DO iw = 1, nw_specfun
              !
              WRITE(stdout, 102) ik, ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2ev * ww(iw), &
                    ryd2mev * esigmar_all(ibnd, ik, iw, itemp), ryd2mev * esigmai_all(ibnd, ik, iw, itemp)
              !
              IF (mpime == ionode_id) &
              WRITE(iospectral_sup, 102) ik, ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2ev * ww(iw), &
                    ryd2mev * esigmar_all(ibnd, ik, iw, itemp), ryd2mev * esigmai_all(ibnd, ik, iw, itemp)
              !
            ENDDO
            !
          ENDDO
          !
          WRITE(stdout, *) ' '
          !
        ENDDO
        !
        IF (mpime == ionode_id) CLOSE(iospectral_sup)
        !
      ENDDO ! itemp
      DEALLOCATE(xkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_pl_q', 'Error deallocating xkf_all', 1)
      DEALLOCATE(etf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('spectral_func_pl_q', 'Error deallocating etf_all', 1)
    ENDIF
    !
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    101 FORMAT(5x, 'ik = ', i7, '  w = ', f9.4, ' eV   A(k,w) = ', e12.5, ' meV^-1')
    102 FORMAT(2i9, 2x, f12.4, 2x, f12.4, 2x, f12.4, 2x, f12.4, 2x, f12.4)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE spectral_func_pl_q
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE a2f_main()
    !-----------------------------------------------------------------------
    !!
    !! Compute the Eliasberg spectral function
    !! in the Migdal approximation.
    !!
    !! If the q-points are not on a uniform grid (i.e. a line)
    !! the function will not be correct
    !!
    !! 02/2009 works in serial on ionode at the moment.  can be parallelized
    !! 03/2009 added transport spectral function -- this involves a v_k dot v_kq term
    !!         in the quantities coming from selfen_phon.f90.  Not fully implemented
    !! 10/2009 the code is transitioning to 'on-the-fly' phonon selfenergies
    !!         and this routine is not currently functional
    !! 10/2015 RM: added calcution of Tc based on Allen-Dynes formula
    !! 09/2019 SP: Cleaning
    !!
    !
    USE kinds,     ONLY : DP
    USE modes,     ONLY : nmodes
    USE cell_base, ONLY : omega
    USE epwcom,    ONLY : degaussq, delta_qsmear, nqsmear, nqstep, nsmear, eps_acustic, &
                          nstemp, delta_smear, degaussw, fsthick, nc
    USE elph2,     ONLY : gtemp, nqtotf, wf, wqf, lambda_all, lambda_v_all
    USE constants_epw, ONLY : ryd2mev, ryd2ev, kelvin2eV, one, two, zero, kelvin2Ry
    USE constants, ONLY : pi
    USE mp,        ONLY : mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    USE io_global, ONLY : ionode_id, stdout
    USE io_var,    ONLY : iua2ffil, iudosfil, iua2ftrfil, iures
    USE io_files,  ONLY : prefix
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 20) :: tp
    CHARACTER(LEN = 256) :: fila2f
    !! File name for Eliashberg spectral function
    CHARACTER(LEN = 256) :: fila2ftr
    !! File name for transport Eliashberg spectral function
    CHARACTER(LEN = 256) :: fildos
    !! File name for phonon density of states
    CHARACTER(LEN = 256) :: filres
    !! File name for resistivity
    !
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: iq
    !! Counter on the q-point index
    INTEGER :: iw
    !! Counter on the frequency
    INTEGER :: ismear
    !! Counter on smearing values (phonons)
    INTEGER :: isig
    !! Counter on smearing values (electrons)
    INTEGER :: i
    !! Counter on mu
    INTEGER :: itemp
    !! Counter on temperature
    INTEGER :: itemprho
    !! Counter on temparture for rho
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: weight
    !! Factor in a2f
    REAL(KIND = DP) :: temp
    !! Temperature
    REAL(KIND = DP) :: n
    !! Carrier density
    REAL(KIND = DP) :: be
    !! Bose-Einstein distribution
    REAL(KIND = DP) :: prefact
    !! Prefactor in resistivity
    REAL(KIND = DP) :: lambda_tot
    !! Total e-ph coupling strength (summation)
    REAL(KIND = DP) :: lambda_tr_tot
    !! Total transport e-ph coupling strength (summation)
    REAL(KIND = DP) :: degaussq0
    !! Phonon smearing
    REAL(KIND = DP) :: inv_degaussq0
    !! Inverse of the smearing for efficiency reasons
    REAL(KIND = DP) :: a2f_tmp
    !! Temporary variable for Eliashberg spectral function
    REAL(KIND = DP) :: a2f_tr_tmp
    !! Temporary variable for transport Eliashberg spectral function
    REAL(KIND = DP) :: om_max
    !! max phonon frequency increased by 10%
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP) :: w0
    !! Current frequency w(imode, iq)
    REAL(KIND = DP) :: l
    !! Temporary variable for e-ph coupling strength
    REAL(KIND = DP) :: l_tr
    !! Temporary variable for transport e-ph coupling strength
    REAL(KIND = DP) :: tc
    !! Critical temperature
    REAL(KIND = DP) :: mu
    !! Coulomb pseudopotential
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    REAL(KIND = DP) :: ww(nqstep)
    !! Current frequency
    REAL(KIND = DP), ALLOCATABLE :: a2f_(:, :)
    !! Eliashberg spectral function for different ismear
    REAL(KIND = DP), ALLOCATABLE :: a2f_tr(:, :)
    !! Transport Eliashberg spectral function for different ismear
    REAL(KIND = DP), ALLOCATABLE :: l_a2f(:)
    !! total e-ph coupling strength (a2f_ integration) for different ismear
    REAL(KIND = DP), ALLOCATABLE :: l_a2f_tr(:)
    !! total transport e-ph coupling strength (a2f_tr integration) for different ismear
    REAL(KIND = DP), ALLOCATABLE :: dosph(:, :)
    !! Phonon density of states for different for different ismear
    REAL(KIND = DP), ALLOCATABLE :: logavg(:)
    !! logavg phonon frequency for different ismear
    REAL(KIND = DP), ALLOCATABLE :: rho(:, :)
    !! Resistivity for different for different ismear
    !
    CALL start_clock('a2F')
    IF (mpime == ionode_id) THEN
      !
      ALLOCATE(a2f_(nqstep, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f_main', 'Error allocating a2f_', 1)
      ALLOCATE(a2f_tr(nqstep, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f_main', 'Error allocating a2f_tr', 1)
      ALLOCATE(dosph(nqstep, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f_main', 'Error allocating dosph', 1)
      ALLOCATE(l_a2f(nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f_main', 'Error allocating l_a2f', 1)
      ALLOCATE(l_a2f_tr(nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f_main', 'Error allocating l_a2f_tr', 1)
      ALLOCATE(logavg(nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f_main', 'Error allocating logavg', 1)
      ! The resitivity is computed for temperature between 0K-1000K by step of 10
      ! This is hardcoded and needs to be changed here if one wants to modify it
      ALLOCATE(rho(100, nqsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('a2f_main', 'Error allocating rho', 1)
      !
      DO itemp = 1, nstemp
        DO isig = 1, nsmear
          !
          WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
          IF (isig < 10) THEN
            WRITE(fila2f,   '(a, a6, i1, a, a)') TRIM(prefix), '.a2f.0', isig, '.', trim(adjustl(tp))
            WRITE(fila2ftr, '(a, a9, i1, a, a)') TRIM(prefix), '.a2f_tr.0', isig, '.', trim(adjustl(tp))
            WRITE(filres,   '(a, a6, i1, a, a)') TRIM(prefix), '.res.0', isig, '.', trim(adjustl(tp))
            WRITE(fildos,   '(a, a8, i1, a, a)') TRIM(prefix), '.phdos.0', isig, '.', trim(adjustl(tp))
          ELSE
            WRITE(fila2f,   '(a, a5, i2, a, a)') TRIM(prefix), '.a2f.', isig, '.', trim(adjustl(tp))
            WRITE(fila2ftr, '(a, a8, i2, a, a)') TRIM(prefix), '.a2f_tr.', isig, '.', trim(adjustl(tp))
            WRITE(filres,   '(a, a5, i2, a, a)') TRIM(prefix), '.res.', isig, '.', trim(adjustl(tp))
            WRITE(fildos,   '(a, a7, i2, a ,a)') TRIM(prefix), '.phdos.', isig, '.', trim(adjustl(tp))
          ENDIF
          OPEN(UNIT = iua2ffil, FILE = fila2f, FORM = 'formatted')
          OPEN(UNIT = iua2ftrfil, FILE = fila2ftr, FORM = 'formatted')
          OPEN(UNIT = iures, FILE = filres, FORM = 'formatted')
          OPEN(UNIT = iudosfil, FILE = fildos, FORM = 'formatted')
          !
          WRITE(stdout, '(/5x, a)') REPEAT('=',67)
          WRITE(stdout, '(5x, "Eliashberg Spectral Function in the Migdal Approximation")')
          WRITE(stdout, '(5x, a/)') REPEAT('=',67)
          !
          om_max = 1.1d0 * MAXVAL(wf(:, :)) ! increase by 10%
          dw = om_max / DBLE(nqstep)
          DO iw = 1, nqstep  !
            ww(iw) = DBLE(iw) * dw
          ENDDO
          !
          lambda_tot    = zero
          l_a2f(:)      = zero
          a2f_(:, :)    = zero
          lambda_tr_tot = zero
          l_a2f_tr(:)   = zero
          a2f_tr(:, :)  = zero
          dosph(:, :)   = zero
          logavg(:)     = zero
          !
          DO ismear = 1, nqsmear
            !
            degaussq0 = degaussq + (ismear - 1) * delta_qsmear
            inv_degaussq0 = one / degaussq0
            !
            DO iw = 1, nqstep  ! loop over points on the a2F(w) graph
              !
              DO iq = 1, nqtotf ! loop over q-points
                DO imode = 1, nmodes ! loop over modes
                  w0 = wf(imode, iq)
                  !
                  IF (w0 > eps_acustic) THEN
                    !
                    l = lambda_all(imode, iq, isig, itemp)
                    IF (lambda_all(imode, iq, isig, itemp) < 0.d0) l = zero ! sanity check
                    !
                    a2f_tmp = wqf(iq) * w0 * l / two
                    !
                    weight = w0gauss((ww(iw) - w0) * inv_degaussq0, 0) * inv_degaussq0
                    a2f_(iw, ismear) = a2f_(iw, ismear) + a2f_tmp * weight
                    dosph(iw, ismear) = dosph(iw, ismear) + wqf(iq) * weight
                    !
                    l_tr = lambda_v_all(imode, iq, isig, itemp)
                    IF (lambda_v_all(imode, iq, isig, itemp) < 0.d0) l_tr = zero !sanity check
                    !
                    a2f_tr_tmp = wqf(iq) * w0 * l_tr / two
                    !
                    a2f_tr(iw, ismear) = a2f_tr(iw, ismear) + a2f_tr_tmp * weight
                    !
                  ENDIF
                ENDDO
              ENDDO
              !
              ! output a2f
              !
              IF (ismear == nqsmear) WRITE(iua2ffil,   '(f12.7, 15f12.7)') ww(iw) * ryd2mev, a2f_(iw, :)
              IF (ismear == nqsmear) WRITE(iua2ftrfil, '(f12.7, 15f12.7)') ww(iw) * ryd2mev, a2f_tr(iw, :)
              IF (ismear == nqsmear) WRITE(iudosfil,   '(f12.7, 15f12.7)') ww(iw) * ryd2mev, dosph(iw, :) / ryd2mev
              !
              ! do the integral 2 int (a2F(w)/w dw)
              !
              l_a2f(ismear) = l_a2f(ismear) + two * a2f_(iw, ismear) / ww(iw) * dw
              l_a2f_tr(ismear) = l_a2f_tr(ismear) + two * a2f_tr(iw, ismear) / ww(iw) * dw
              logavg(ismear) = logavg(ismear) + two *  a2f_(iw, ismear) * LOG(ww(iw)) / ww(iw) * dw
              !
            ENDDO
            !
            logavg(ismear) = EXP(logavg(ismear) / l_a2f(ismear))
            !
          ENDDO
          !
          DO iq = 1, nqtotf ! loop over q-points
            DO imode = 1, nmodes ! loop over modes
              IF (lambda_all(imode, iq, isig, itemp) > 0.d0 .AND. wf(imode, iq) > eps_acustic ) &
                lambda_tot = lambda_tot + wqf(iq) * lambda_all(imode, iq, isig, itemp)
              IF (lambda_v_all(imode, iq, isig, itemp) > 0.d0 .AND. wf(imode, iq) > eps_acustic) &
                lambda_tr_tot = lambda_tr_tot + wqf(iq) * lambda_v_all(imode, iq, isig, itemp)
            ENDDO
          ENDDO
          WRITE(stdout, '(5x, a, f12.7)') "lambda : ", lambda_tot
          WRITE(stdout, '(5x, a, f12.7)') "lambda_tr : ", lambda_tr_tot
          WRITE(stdout, '(a)') " "
          !
          !
          ! Allen-Dynes estimate of Tc for ismear = 1
          !
          WRITE(stdout, '(5x, a, f12.7, a)') "Estimated Allen-Dynes Tc"
          WRITE(stdout, '(a)') " "
          WRITE(stdout, '(5x, a, f12.7, a, f12.7)') "logavg = ", logavg(1), " l_a2f = ", l_a2f(1)
          DO i = 1, 6
            !
            mu = 0.1d0 + 0.02d0 * DBLE(i - 1)
            tc = logavg(1) / 1.2d0 * EXP(-1.04d0 * (1.d0 + l_a2f(1)) / (l_a2f(1) - mu * ( 1.d0 + 0.62d0 * l_a2f(1))))
            ! tc in K
            !
            tc = tc * ryd2ev / kelvin2eV
            !SP: IF Tc is too big, it is not physical
            IF (tc < 1000.0) THEN
              WRITE(stdout, '(5x, a, f6.2, a, f22.12, a)') "mu = ", mu, " Tc = ", tc, " K"
            ENDIF
            !
          ENDDO
          !
          rho(:, :) = zero
          ! Now compute the Resistivity of Metal using the Ziman formula
          ! rho(T,smearing) = 4 * pi * me/(n * e**2 * kb * T) int dw hbar w a2F_tr(w,smearing) n(w,T)(1+n(w,T))
          ! n is the number of electron per unit volume and n(w,T) is the Bose-Einstein distribution
          ! Usually this means "the number of electrons that contribute to the mobility" and so it is typically 8 (full shell)
          ! but not always. You might want to check this.
          !
          n = nc / omega
          WRITE(iures, '(a)') '# Temperature [K]                &
                              Resistivity [micro Ohm cm] for different Phonon smearing (meV)        '
          WRITE(iures, '("#     ", 15f12.7)') ((degaussq + (ismear - 1) * delta_qsmear) * ryd2mev, ismear = 1, nqsmear)
          DO ismear = 1, nqsmear
            DO itemprho = 1, 100 ! Per step of 10K
              temp = itemprho * 10.d0 * kelvin2Ry
              ! omega is the volume of the primitive cell in a.u.
              !
              prefact = 4.d0 * pi / (temp * n)
              DO iw = 1, nqstep  ! loop over points on the a2F(w)
                !
                be = one / (EXP(ww(iw) / temp) - one)
                ! Perform the integral with rectangle.
                rho(itemprho, ismear) = rho(itemprho, ismear) + prefact * ww(iw) * a2f_tr(iw, ismear) * be * (1.d0 + be) * dw
                !
              ENDDO
              ! From a.u. to micro Ohm cm
              ! Conductivity 1 a.u. = 2.2999241E6 S/m
              ! Now to go from Ohm*m to micro Ohm cm we need to multiply by 1E8
              rho(itemprho, ismear) = rho(itemprho, ismear) * 1E8 / 2.2999241E6
              IF (ismear == nqsmear) WRITE (iures, '(i8, 15f12.7)') itemprho * 10, rho(itemprho, :)
            ENDDO
          ENDDO
          CLOSE(iures)
          !
          WRITE(iua2ffil, *) "Integrated el-ph coupling"
          WRITE(iua2ffil, '("  #         ", 15f12.7)') l_a2f(:)
          WRITE(iua2ffil, *) "Phonon smearing (meV)"
          WRITE(iua2ffil, '("  #         ", 15f12.7)') ((degaussq + (ismear - 1) * delta_qsmear) * ryd2mev, ismear = 1, nqsmear)
          WRITE(iua2ffil, '(" Electron smearing (eV)", f12.7)') ((isig - 1) * delta_smear + degaussw) * ryd2ev
          WRITE(iua2ffil, '(" Fermi window (eV)", f12.7)') fsthick * ryd2ev
          WRITE(iua2ffil, '(" Summed el-ph coupling ", f12.7)') lambda_tot
          CLOSE(iua2ffil)
          !
          WRITE(iua2ftrfil, *) "Integrated el-ph coupling"
          WRITE(iua2ftrfil, '("  #         ", 15f12.7)') l_a2f_tr(:)
          WRITE(iua2ftrfil, *) "Phonon smearing (meV)"
          WRITE(iua2ftrfil, '("  #         ", 15f12.7)') ((degaussq + (ismear - 1) * delta_qsmear) * ryd2mev, ismear = 1, nqsmear)
          WRITE(iua2ftrfil, '(" Electron smearing (eV)", f12.7)') ((isig - 1) * delta_smear + degaussw) * ryd2ev
          WRITE(iua2ftrfil, '(" Fermi window (eV)", f12.7)') fsthick * ryd2ev
          WRITE(iua2ftrfil, '(" Summed el-ph coupling ", f12.7)') lambda_tot
          CLOSE(iua2ftrfil)
          !
          CLOSE(iudosfil)
          !
        ENDDO ! isig
        !
      ENDDO ! itemp
      DEALLOCATE(l_a2f, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating l_a2f', 1)
      DEALLOCATE(l_a2f_tr, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating l_a2f_tr', 1)
      DEALLOCATE(a2f_, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating a2f', 1)
      DEALLOCATE(a2f_tr, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating a2f_tr', 1)
      DEALLOCATE(rho, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating rho', 1)
      DEALLOCATE(dosph, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating dosph', 1)
      DEALLOCATE(logavg, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_a2f', 'Error deallocating logavg', 1)
      !
    ENDIF
    !
    CALL stop_clock('a2F')
    CALL print_clock('a2F')
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE a2f_main
    !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  END MODULE spectral_func
  !-----------------------------------------------------------------------

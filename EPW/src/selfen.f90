  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE selfen
  !----------------------------------------------------------------------
  !!
  !! This module contains the various self-energy routines
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE selfen_elec_q(iqq, iq, totq, first_cycle)
    !-----------------------------------------------------------------------
    !!
    !!  Compute the imaginary part of the electron self energy due to electron-
    !!  phonon interaction in the Migdal approximation. This corresponds to
    !!  the electron linewidth (half width). The phonon frequency is taken into
    !!  account in the energy selection rule.
    !!
    !!  Use matrix elements, electronic eigenvalues and phonon frequencies
    !!  from ep-wannier interpolation
    !!
    !!  This routines computes the contribution from phonon iq to all k-points
    !!  The outer loop in ephwann_shuffle.f90 will loop over all iq points
    !!  The contribution from each iq is summed at the end of this subroutine
    !!  for iqq=totq to recover the per-ik electron self energy
    !!
    !!  RM 24/02/2014
    !!  Redefined the size of sigmar_all, sigmai_all, and zi_all within the fermi windwow
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE modes,         ONLY : nmodes
    USE input,         ONLY : nstemp, fsthick, ngaussw, degaussw, &
                              eps_acoustic, efermi_read, fermi_energy, restart, restart_step, &
                              lwfpt, ahc_win_min, ahc_win_max, elecselfen_type, specfun_el, &
                              elecselfen, wmin_specfun, wmax_specfun, nw_specfun, &
                              lfast_kmesh, nqf1, nqf2, nqf3
    USE pwcom,         ONLY : ef
    USE global_var,    ONLY : etf, ibndmin, xqf, eta, nbndfst, &
                              nkf, epf17, wf, wqf, adapt_smearing, &
                              sigmar_all, sigmai_all, sigmai_mode, zi_all, efnew, &
                              nktotf, lower_bnd, gtemp, dwf17, esigmar_all, esigmai_all, &
                              sigmar_dw_all
    USE control_flags, ONLY : iverbosity
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, ryd2ev, one, two, zero, ci, eps6, eps8
    USE constants,     ONLY : pi
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, inter_image_comm
    USE io_selfen,     ONLY : selfen_el_write, selfen_el_write_wfpt, spectral_write
    USE parallelism,   ONLY : poolgather2
    USE utilities,     ONLY : fermi_dirac
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: first_cycle
    !! Use to determine weather this is the first cycle after restart
    INTEGER, INTENT(in) :: iqq
    !! Q-point index from selecq.fmt window
    INTEGER, INTENT(in) :: iq
    !! Q-point index from full grid
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points from the selecq.fmt grid.
    !
    ! Local variables
    CHARACTER(LEN = 10) :: broadening_method
    !! Function to use for the broadened delta function. Loretnzian or Gaussian.
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ik_global
    !! Global k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: iw
    !! Counter on the frequency
    !
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared in Ry^2
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: ekk
    !! Eigen energy at k on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy at k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: etmp1
    !! Temporary variable to store etmp1 = ekk - (ekq - wq)
    REAL(KIND = DP) :: etmp2
    !! Temporary variable to strore etmp2 = ekk - (ekq + wq)
    REAL(KIND = DP) :: etmpw1
    !! Temporary variable to store etmpw1 = ww - (ekq - wq)
    REAL(KIND = DP) :: etmpw2
    !! Temporary variable to store etmpw2 = ww - (ekq + wq)
    REAL(KIND = DP) :: sq_etmp1
    !! Temporary variable to store etmp1^2
    REAL(KIND = DP) :: sq_etmp2
    !! Temporary variable to store etmp2^2
    REAL(KIND = DP) :: wgkq
    !! Fermi-Dirac occupation factor $f_{mk+q}(T)$
    REAL(KIND = DP) :: fact1
    !! Temporary variable to store $f_{mk+q}(T) + n_{q\nu}(T)$
    REAL(KIND = DP) :: fact2
    !! Temporary variable to store $1 - f_{mk+q}(T) + n_{q\nu}(T)$
    REAL(KIND = DP) :: weight
    !! Self-energy factor
    !!$$ N_q \Re(\frac{f_{mk+q}(T) + n_{q\nu}(T)}{\varepsilon_{nk} - \varepsilon_{mk+q} + \omega_{q\nu} - i\delta}) $$
    !!$$ + N_q \Re(\frac{1 - f_{mk+q}(T) + n_{q\nu}(T)}{\varepsilon_{nk} - \varepsilon_{mk+q} - \omega_{q\nu} - i\delta}) $$
    REAL(KIND = DP) :: weight0
    !! Self-energy factor at zero frequency
    REAL(KIND = DP) :: w0g1
    !! Dirac delta at k for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: w0g2
    !! Dirac delta at k+q for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse of degaussw define for efficiency reasons
    REAL(KIND = DP) :: eta_tmp
    !! Temporary variable eta2
    REAL(KIND = DP) :: sq_eta_tmp
    !! Temporary eta2^2
    REAL(KIND = DP) :: inv_eta_tmp
    !! Temporary varialbe inv_eta
    REAL(KIND = DP) :: dw
    !! Frequency intervals
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    REAL(KIND = DP) :: wq(nmodes)
    !! Phonon frequency on the fine grid
    REAL(KIND = DP) :: inv_wq(nmodes)
    !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
    REAL(KIND = DP) :: g2_tmp(nmodes)
    !! If the phonon frequency is too small discart g
    REAL(KIND = DP) :: wgq(nmodes)
    !! Bose occupation factor $n_{q\nu}(T)$
    REAL(KIND = DP) :: eta2(nbndfst, nmodes, nktotf)
    !! Temporary array to store the current smearing eta
    REAL(KIND = DP) :: inv_eta(nbndfst, nmodes, nktotf)
    !! Temporary array to store the inverse of the eta for speed purposes
    REAL(KIND = DP) :: ww(nw_specfun)
    !! Current frequency
    REAL(KIND = DP) :: wqf_loc
    !! Local q-point weight
    COMPLEX(KIND = DP) :: fact
    !! Self-energy factor
    !
    CALL start_clock('selfen_elec_q')
    !
    IF (lwfpt) THEN
      broadening_method = "Lorentzian"
    ELSE
      broadening_method = "Gaussian"
    ENDIF
    !
    ! Weight of the q-points
    IF (lfast_kmesh) THEN
      wqf_loc = 1.0d0 / REAL(nqf1 * nqf2 * nqf3, KIND = DP)
    ELSE
      wqf_loc = wqf(iq)
    ENDIF
    !
    ! energy range and spacing for spectral function
    !
    IF (specfun_el) THEN
      dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
      DO iw = 1, nw_specfun
        ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
      ENDDO
    ENDIF
    !
    ! SP: Define the inverse so that we can efficiently multiply instead of dividing
    inv_degaussw = one / degaussw
    ! To avoid if branching in the loop
    inv_eta(:, :, :) = zero
    IF (adapt_smearing) THEN
      DO ik = 1, nkf
        DO ibnd = 1, nbndfst
          DO imode = 1, nmodes
            inv_eta(ibnd, imode, ik) = one / (DSQRT(two) * eta(imode, ibnd, ik))
            eta2(ibnd, imode, ik) = DSQRT(two) * eta(imode, ibnd, ik)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO ik = 1, nkf
        DO ibnd = 1, nbndfst
          DO imode = 1, nmodes
            inv_eta(ibnd, imode, ik) = inv_degaussw
            eta2(ibnd, imode, ik) = degaussw
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    IF (iqq == 1) THEN
      !
      ! Print calculation information to stdout
      !
      IF (.NOT. adapt_smearing) THEN
        WRITE(stdout, '(5x, a, " Broadening: ", f10.6, " eV, ngauss=", i4)') &
            broadening_method, degaussw * ryd2ev, ngaussw
        WRITE(stdout, '(a)') ' '
      ELSE
        WRITE(stdout, '(5x, a)') "Adaptive broadening"
      ENDIF
      !
      IF (elecselfen) THEN
        WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
        WRITE(stdout, '(5x, "Electron (Imaginary) Self-Energy in the Migdal Approximation")')
        WRITE(stdout, '(5x, a/)') REPEAT('=', 67)
      ENDIF
      !
      IF (specfun_el) THEN
        WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
        WRITE(stdout, '(5x, "Electron Spectral Function in the Migdal Approximation")')
        WRITE(stdout, '(5x, a/)') REPEAT('=', 67)
        !
        IF (lwfpt) THEN
          WRITE(stdout, '(a)') ' '
          WRITE(stdout, '(5x, a)') 'The sum rule to conserve the number of electron is NOT enforced.'
          WRITE(stdout, '(5x, a)') 'The Debye-Waller and upper Fan term are calculated using WFPT.'
          WRITE(stdout, '(a)') ' '
        ELSE
          WRITE(stdout, '(a)') ' '
          WRITE(stdout, '(5x, a)') 'The sum rule to conserve the number of electron is enforced.'
          WRITE(stdout, '(5x, a)') 'The self energy is rescaled so that its real part is zero at the Fermi level.'
          WRITE(stdout, '(5x, a)') 'The sum rule replace the explicit calculation of the Debye-Waller term.'
          WRITE(stdout, '(a)') ' '
        ENDIF
      ENDIF
      !
      IF (fsthick < 1.d3) WRITE(stdout, '(/5x, a, f10.6, a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
    ENDIF
    !
    DO itemp = 1, nstemp ! loop over temperatures
      !
      ! Now pre-treat phonon modes for efficiency
      ! Treat phonon frequency and Bose occupation
      wq(:) = zero
      DO imode = 1, nmodes
        IF (lfast_kmesh) THEN
          wq(imode) = wf(imode, iqq)
        ELSE
          wq(imode) = wf(imode, iq)
        ENDIF
        IF (wq(imode) > eps_acoustic) THEN
          g2_tmp(imode) = one
          wgq(imode)    = fermi_dirac(wq(imode), gtemp(itemp))
          wgq(imode)    = wgq(imode) / (one - two * wgq(imode))
          inv_wq(imode) = one / (two * wq(imode))
        ELSE
          g2_tmp(imode) = zero
          wgq(imode)    = zero
          inv_wq(imode) = zero
        ENDIF
      ENDDO
      !
      IF (iqq == 1) THEN
        !
        WRITE(stdout, '(/5x, a, f10.6, a)') 'Golden Rule strictly enforced with T = ', gtemp(itemp) * ryd2ev, ' eV'
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
      IF (restart) THEN
        ! Make everythin 0 except the range of k-points we are working on
        sigmar_all(:, 1:lower_bnd - 1, :) = zero
        sigmar_all(:, lower_bnd + nkf:nktotf, :) = zero
        sigmai_all(:, 1:lower_bnd - 1, :) = zero
        sigmai_all(:, lower_bnd + nkf:nktotf, :) = zero
        zi_all(:, 1:lower_bnd - 1, :) = zero
        zi_all(:, lower_bnd + nkf:nktotf, :) = zero
        IF (lwfpt) THEN
          sigmar_dw_all(:, 1:lower_bnd - 1, :) = zero
          sigmar_dw_all(:, lower_bnd + nkf:nktotf, :) = zero
        ENDIF
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
        DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          ik_global = ik + lower_bnd - 1
          !
          ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
          ! (but in this case they are the same)
          !
          IF (lwfpt .OR. ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
                          (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick))) THEN
            !
            DO imode = 1, nmodes
              !
              DO ibnd = 1, nbndfst
                !
                ! the energy of the electron at k (relative to Ef)
                ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
                !
                eta_tmp     = eta2(ibnd, imode, ik)
                sq_eta_tmp  = eta_tmp**two
                inv_eta_tmp = inv_eta(ibnd, imode, ik)
                !
                DO jbnd = 1, nbndfst
                  !
                  ! the energy of the electron at k+q (relative to Ef)
                  ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
                  !
                  IF (lwfpt) THEN
                    !
                    ! Skip coupling with oneself or between degenerate states at the same k point
                    IF (ALL(xqf(:, iq) < eps8) .AND. ABS(ekq - ekk) < 2.d-5) CYCLE
                    ! Skip active states outside the ahc window
                    IF (ekq + ef0 < ahc_win_min .OR. ekq + ef0 > ahc_win_max) CYCLE
                    !
                  ENDIF
                  !
                  ! the Fermi occupation at k+q
                  wgkq = fermi_dirac(ekq, gtemp(itemp))
                  !
                  ! here we take into account the zero-point DSQRT(hbar/2M\omega)
                  ! with hbar = 1 and M already contained in the eigenmodes
                  ! g2 is Ry^2, wkf must already account for the spin factor
                  !
                  g2 = (ABS(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq(imode) * g2_tmp(imode)
                  !
                  ! There is a sign error for wq in Eq. 9 of Comp. Phys. Comm. 181, 2140 (2010). - RM
                  ! The sign was corrected according to Eq. (7.282) page 489 from Mahan's book
                  ! (Many-Particle Physics, 3rd edition)
                  !
                  fact1 =       wgkq + wgq(imode)
                  fact2 = one - wgkq + wgq(imode)
                  !
                  IF (elecselfen_type == 'adiabatic') THEN
                    etmp1 = ekk - ekq
                    etmp2 = ekk - ekq
                  ELSEIF (elecselfen_type == 'nonadiabatic') THEN
                    etmp1 = ekk - (ekq - wq(imode))
                    etmp2 = ekk - (ekq + wq(imode))
                  ENDIF
                  !
                  ! Self-energy at the bare electron energy
                  !
                  IF (elecselfen) THEN
                    !
                    ! Real part of self-energy
                    !
                    weight = wqf_loc * REAL(fact1 / (etmp1 - ci * eta_tmp) + fact2 / (etmp2 - ci * eta_tmp))
                    !
                    ! \Re\Sigma [Eq. 3 in Comput. Phys. Commun. 209, 116 (2016)]
                    sigmar_all(ibnd, ik_global, itemp) = sigmar_all(ibnd, ik_global, itemp) + g2 * weight
                    !
                    ! Imaginary part of self-energy
                    !
                    IF (broadening_method == "Lorentzian") THEN
                      ! Logical implementation
                      ! weight = wqf_loc * aimag(                                                  &
                      !         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
                      !           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) )
                      !
                      weight = wqf_loc * (  fact1 * eta_tmp / (etmp1**2 + eta_tmp**2) &
                                          + fact2 * eta_tmp / (etmp2**2 + eta_tmp**2) )
                    ELSEIF (broadening_method == "Gaussian") THEN
                      ! Gaussian broadening of delta function
                      w0g1 = w0gauss(etmp1 * inv_eta_tmp, 0) * inv_eta_tmp
                      w0g2 = w0gauss(etmp2 * inv_eta_tmp, 0) * inv_eta_tmp
                      !
                      weight = pi * wqf_loc * (fact1 * w0g1 + fact2 * w0g2)
                    ELSE
                      CALL errore("selfen_elec_q", "Wrong broadening_method", 1)
                    ENDIF
                    !
                    sigmai_all(ibnd, ik_global, itemp) = sigmai_all(ibnd, ik_global, itemp) + g2 * weight
                    !
                    ! Mode-resolved
                    IF (iverbosity == 3) THEN
                      sigmai_mode(ibnd, imode, ik_global, itemp) = &
                      sigmai_mode(ibnd, imode, ik_global, itemp) + g2 * weight
                    ENDIF
                    !
                    ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
                    !
                    sq_etmp1 = etmp1 * etmp1
                    sq_etmp2 = etmp2 * etmp2
                    !
                    weight = wqf_loc * &
                             (fact1 * (sq_etmp1 - sq_eta_tmp) / (sq_etmp1 + sq_eta_tmp)**two +  &
                              fact2 * (sq_etmp2 - sq_eta_tmp) / (sq_etmp2 + sq_eta_tmp)**two)
                    !
                    zi_all(ibnd, ik_global, itemp) = zi_all(ibnd, ik_global, itemp) + g2 * weight
                    !
                  ENDIF ! elecselfen
                  !
                  ! Frequency-dependent self-energy
                  ! See Eq. 3 in Comput. Phys. Commun. 209, 116 (2016)
                  !
                  IF (specfun_el) THEN
                    !
                    etmp1 = -(ekq - wq(imode))
                    etmp2 = -(ekq + wq(imode))
                    !
                    ! Self-energy weight at the chemical potential (w = 0)
                    weight0 = REAL(fact1 / (etmp1 - ci * degaussw) + fact2 / (etmp2 - ci * degaussw))
                    !
                    DO iw = 1, nw_specfun
                      !
                      etmpw1 = ww(iw) + etmp1
                      etmpw2 = ww(iw) + etmp2
                      !
                      fact = fact1 / (etmpw1 - ci * degaussw) + fact2 / (etmpw2 - ci *  degaussw)
                      !
                      ! \Re\Sigma
                      !
                      weight = REAL(fact)
                      !
                      IF (.NOT. lwfpt) THEN
                        ! SP : Application of the sum rule
                        ! If using WFPT, do not apply the sum rule because the Debye-Waller
                        ! term is explicitly calculated.
                        weight = weight - weight0
                      ENDIF
                      !
                      esigmar_all(ibnd, ik_global, iw, itemp) = &
                      esigmar_all(ibnd, ik_global, iw, itemp) + wqf_loc * g2 * weight
                      !
                      ! \Im\Sigma
                      !
                      weight = AIMAG(fact)
                      !
                      esigmai_all(ibnd, ik_global, iw, itemp) = &
                      esigmai_all(ibnd, ik_global, iw, itemp) + wqf_loc * g2 * weight
                      !
                    ENDDO
                  ENDIF ! specfun_el
                  !
                ENDDO !jbnd
              ENDDO !ibnd
              !
              ! Active space Debye-Waller term
              !
              IF (lwfpt) THEN
                !
                weight = (wgq(imode) + one / two) * inv_wq(imode) * wqf_loc
                !
                IF (elecselfen) THEN
                  DO ibnd = 1, nbndfst
                    sigmar_dw_all(ibnd, ik_global, itemp) = sigmar_dw_all(ibnd, ik_global, itemp) &
                      + weight * REAL(dwf17(ibnd, ibnd, imode, ik))
                  ENDDO ! ibnd
                ENDIF
                !
                IF (specfun_el) THEN
                  DO ibnd = 1, nbndfst
                    DO iw = 1, nw_specfun
                      esigmar_all(ibnd, ik_global, iw, itemp) = &
                      esigmar_all(ibnd, ik_global, iw, itemp) + weight * REAL(dwf17(ibnd, ibnd, imode, ik))
                    ENDDO
                  ENDDO ! ibnd
                ENDIF
                !
              ENDIF
              !
            ENDDO !imode
          ENDIF ! endif  fsthick
        ENDDO ! end loop on k
        !
        ! Creation of a restart point
        IF (restart) THEN
          IF (MOD(iqq, restart_step) == 0 .and. itemp == nstemp) THEN
            WRITE(stdout, '(5x, a, i10)' ) 'Creation of a restart point at ', iqq
            IF (elecselfen) THEN
              CALL mp_sum(sigmar_all, inter_pool_comm)
              CALL mp_sum(sigmai_all, inter_pool_comm)
              CALL mp_sum(zi_all, inter_pool_comm)
              IF (lwfpt) THEN
                CALL mp_sum(sigmar_dw_all, inter_pool_comm)
                CALL selfen_el_write_wfpt(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all, sigmar_dw_all)
              ELSE
                CALL selfen_el_write(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all)
              ENDIF
            ENDIF
            IF (specfun_el) THEN
              CALL mp_sum(esigmar_all, inter_pool_comm)
              CALL mp_sum(esigmai_all, inter_pool_comm)
              CALL spectral_write(iqq, totq, nktotf, esigmar_all, esigmai_all)
            ENDIF
          ENDIF
        ENDIF
      ENDIF ! in case of restart, do not do the first one
    ENDDO ! itemp
    !
    CALL stop_clock('selfen_elec_q')
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE selfen_elec_q
    !-----------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE selfen_elec_print
    !--------------------------------------------------------------------------
    !! Collect self-energy and print them to stdout and file
    !--------------------------------------------------------------------------
    !
    USE, INTRINSIC :: ieee_arithmetic, ONLY: IEEE_VALUE, IEEE_QUIET_NAN
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : ryd2mev, ryd2ev, czero, zero, eps6, one, kelvin2eV
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, inter_image_comm
    USE io_global,     ONLY : stdout, ionode
    USE io_var,        ONLY : linewidth_elself, iuelself_wfpt
    USE input,         ONLY : nbndsub, efermi_read, fermi_energy, nstemp, &
                              ahc_win_min, ahc_win_max, lwfpt
    USE control_flags, ONLY : iverbosity
    USE modes,         ONLY : nmodes
    USE global_var,    ONLY : etf, ibndmin, nkqf, nbndfst, xkf, nkqtotf, gtemp, &
                              sigma_ahc_hdw, sigma_ahc_uf, nktotf, efnew, sigmar_all, &
                              sigmai_all, zi_all, sigmai_mode, sigmar_dw_all
    USE utilities,     ONLY : degenerate_average_real
    USE parallelism,   ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 20) :: tp
    !! string for temperatures
    CHARACTER(LEN = 256) :: fileselfen
    !! file name for self energy
    INTEGER :: ierr
    !! Error status
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! k+q-point index
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: itemp
    !! Counter on temperatures
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: ekk
    !! Eigen energy at k on the fine grid relative to the Fermi level
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Collect eigenenergies from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: sigmar_all_sum(:, :, :)
    !! Real part of the total self-energy, temporarily used for printing
    !
    ! Fermi level
    !
    IF (efermi_read) THEN
      ef0 = fermi_energy
    ELSE
      ef0 = efnew
    ENDIF
    !
    ! The k points are distributed among pools: here we collect them
    !
    ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('selfen_elec_print', 'Error allocating xkf_all', 1)
    ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('selfen_elec_print', 'Error allocating etf_all', 1)
    xkf_all(:, :) = zero
    etf_all(:, :) = zero
    !
    ! note that poolgather2 works with the doubled grid (k and k+q)
    CALL poolgather2(3, nkqtotf, nkqf, xkf, xkf_all)
    CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
    !
    ! Collect k-point distributed to pools.
    !
    CALL mp_sum(sigmar_all, inter_pool_comm)
    CALL mp_sum(sigmai_all, inter_pool_comm)
    CALL mp_sum(zi_all, inter_pool_comm)
    CALL mp_sum(sigmar_all, inter_image_comm)
    CALL mp_sum(sigmai_all, inter_image_comm)
    CALL mp_sum(zi_all, inter_image_comm)
    IF (lwfpt) THEN
      CALL mp_sum(sigmar_dw_all, inter_pool_comm)
      CALL mp_sum(sigma_ahc_hdw, inter_pool_comm)
      CALL mp_sum(sigma_ahc_uf, inter_pool_comm)
    ENDIF
    IF (iverbosity == 3) THEN 
      CALL mp_sum(sigmai_mode, inter_pool_comm)
      CALL mp_sum(sigmai_mode, inter_image_comm)
    ENDIF
    !
    ! Average over degenerate eigenstates
    !
    WRITE(stdout, '(5x,"Average over degenerate eigenstates is performed")')
    CALL degenerate_average_real(sigmar_all, etf_all)
    CALL degenerate_average_real(sigmai_all, etf_all)
    CALL degenerate_average_real(zi_all, etf_all)
    IF (lwfpt) THEN
      CALL degenerate_average_real(sigmar_dw_all, etf_all)
      CALL degenerate_average_real(sigma_ahc_hdw, etf_all)
      CALL degenerate_average_real(sigma_ahc_uf, etf_all)
    ENDIF
    !
    ALLOCATE(sigmar_all_sum(nbndfst, nktotf, nstemp), STAT=ierr)
    IF (ierr /=0) CALL errore('selfen_elec_print', 'Error allocating sigmar_all_sum', 1)
    !
    sigmar_all_sum(:,:,:) = 0.d0
    !
    IF (lwfpt) THEN
      sigmar_all_sum = sigmar_all + sigmar_dw_all + sigma_ahc_hdw + sigma_ahc_uf
    ELSE
      sigmar_all_sum = sigmar_all
    ENDIF
    !
    IF (lwfpt) THEN
      WRITE(stdout, '(5x, "Electron Self-Energy using Wannier function perturbation theory")')
      !
      ! Upper Fan self-energy is valid only if the energy is inside the AHC window.
      ! Set it to NaN otherwise.
      !
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          ikk = 2 * ik - 1
          DO ibnd = 1, nbndfst
            !
            ekk = etf_all(ibndmin - 1 + ibnd, ikk)
            IF (ekk < ahc_win_min .OR. ahc_win_max < ekk) THEN
              ! sigma_ahc_uf(ibnd, ik, itemp) = IEEE_VALUE(sigma_ahc_uf(ibnd, ik, itemp), IEEE_QUIET_NAN)
              ! QE testcode cannot deal with NaNs. So, I set the values to zero instead of NaN.
              ! FIXME: Delete the next two lines once gitlab.com/QEF/q-e/-/issues/623 is fixed.
              sigma_ahc_uf(ibnd, ik, itemp) = zero
              sigmar_all_sum(ibnd, ik, itemp) = zero
            ENDIF
            !
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    ! Output electron self-energy
    !
    DO itemp = 1, nstemp
      WRITE(stdout, '(5x, a, f8.3, a)') "Temperature: ", gtemp(itemp) * ryd2ev / kelvin2eV, "K"
      !
      ! Output electron SE here after looping over all q-points (with their contributions
      ! summed in sigmar_all, etc.)
      !
      WRITE(stdout, '(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
      !
      IF (ionode) THEN
        ! Write to file
        WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
        fileselfen = 'linewidth.elself.' // trim(adjustl(tp)) // 'K'
        OPEN(UNIT = linewidth_elself, FILE = fileselfen)
        WRITE(linewidth_elself, '(a)') '# Electron linewidth = 2*Im(Sigma) (meV)'
        IF (iverbosity == 3) THEN
          WRITE(linewidth_elself, '(a)') '#      ik       ibnd                 E(ibnd)      imode          Im(Sigma)(meV)'
        ELSE
          WRITE(linewidth_elself, '(a)') '#      ik       ibnd                 E(ibnd)      Im(Sigma)(meV)'
        ENDIF
        !
        DO ik = 1, nktotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          WRITE(stdout, '(/5x, "ik = ", i7," coord.: ", 3f12.7)') ik, xkf_all(:, ikk)
          WRITE(stdout, '(5x, a)') REPEAT('-', 67)
          !
          DO ibnd = 1, nbndfst
            !
            ! note that ekk does not depend on q
            ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
            !
            ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
            zi_all(ibnd, ik, itemp) = one / (one + zi_all(ibnd, ik, itemp))
            !
            WRITE(stdout, 102) ibndmin - 1 + ibnd, &
              ryd2ev * ekk, &
              ryd2mev * sigmar_all_sum(ibnd, ik, itemp), &
              ryd2mev * sigmai_all(ibnd,ik, itemp), &
              zi_all(ibnd, ik, itemp), &
              one / zi_all(ibnd, ik, itemp) - one
            !
            IF (iverbosity == 3) THEN
              DO imode = 1, nmodes
                WRITE(linewidth_elself, '(i9, 2x)', ADVANCE = 'no') ik
                WRITE(linewidth_elself, '(i9, 2x)', ADVANCE = 'no') ibndmin - 1 + ibnd
                WRITE(linewidth_elself, '(E22.14, 2x)', ADVANCE = 'no') ryd2ev * ekk
                WRITE(linewidth_elself, '(i9, 2x)', ADVANCE = 'no') imode
                WRITE(linewidth_elself, '(E22.14, 2x)') ryd2mev * sigmai_mode(ibnd, imode, ik, itemp)
              ENDDO
            ELSE
              WRITE(linewidth_elself, '(i9, 2x)', ADVANCE = 'no') ik
              WRITE(linewidth_elself, '(i9, 2x)', ADVANCE = 'no') ibndmin - 1 + ibnd
              WRITE(linewidth_elself, '(E22.14, 2x)', ADVANCE = 'no') ryd2ev * ekk
              WRITE(linewidth_elself, '(E22.14, 2x)') ryd2mev * sigmai_all(ibnd, ik, itemp)
            ENDIF
            !
          ENDDO
          WRITE(stdout, '(5x, a/)') REPEAT('-', 67)
          !
        ENDDO
        CLOSE(linewidth_elself)
      ENDIF ! inode
      !
      ! Print self-energy and Z factor to stdout
      !
      DO ibnd = 1, nbndfst
        DO ik = 1, nktotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          ! note that ekk does not depend on q
          ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
          !
          ! calculate Z = 1 / (1 - \frac{\partial\Sigma}{\partial\omega})
          !zi_all(ibnd,ik) = one / (one + zi_all(ibnd,ik))
          !
          WRITE(stdout, '(2i9, 5E22.14)') ik, ibndmin - 1 + ibnd, &
            ryd2ev * ekk,&
            ryd2mev * sigmar_all_sum(ibnd, ik, itemp), &
            ryd2mev * sigmai_all(ibnd, ik, itemp), &
            zi_all(ibnd, ik, itemp), &
            one / zi_all(ibnd, ik, itemp) - one
          !
        ENDDO
        !
        WRITE(stdout, '(a)') '  '
        !
      ENDDO
      !
      ! Print WFPT output self-energy to stdout
      !
      IF (lwfpt) THEN
        !
        WRITE(stdout, '(a)') ''
        WRITE(stdout, '(5x,a)') 'Full decomposition of the Allen-Heine-Cardona self-energy into the Fan/Debye-Waller'
        WRITE(stdout, '(5x,a)') 'and active-space/rest-space contributions is written to file elself_wfpt_sup.#K'
        WRITE(stdout, '(a)') ''
        WRITE(stdout, '(a)') ''
        !
        WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
        fileselfen = 'elself_wfpt_sup.' // trim(adjustl(tp)) // 'K'
        OPEN(UNIT = iuelself_wfpt, FILE = fileselfen)
        !
        WRITE(iuelself_wfpt, '(a)') '# Electron self-energy (meV) in the Allen-Heine-Cardona&
                          & formalism at T = ' // trim(adjustl(tp)) // 'K'
        WRITE(iuelself_wfpt, '(a)') '#   ik  ibnd      E_nk (eV) Re[Active_Fan]      Active_DW&
                          &       Rest_Fan        Rest_DW Im[Active_Fan] (meV)'
        !
        DO ibnd = 1, nbndfst
          DO ik = 1, nktotf
            !
            ikk = 2 * ik - 1
            !
            ekk = etf_all(ibndmin - 1 + ibnd, ikk)
            !
            WRITE(iuelself_wfpt, '(2i6,6E22.14)') ik, ibndmin - 1 + ibnd, &
              ryd2ev * ekk, &
              ryd2mev * sigmar_all(ibnd, ik, itemp), &
              ryd2mev * sigmar_dw_all(ibnd, ik, itemp), &
              ryd2mev * sigma_ahc_uf(ibnd, ik, itemp), &
              ryd2mev * sigma_ahc_hdw(ibnd, ik, itemp), &
              ryd2mev * sigmai_all(ibnd, ik, itemp)
          ENDDO ! ik
          !
          WRITE(iuelself_wfpt, '(a)') '  '
          !
        ENDDO ! ibnd
        !
        CLOSE(iuelself_wfpt)
        !
      ENDIF ! lwfpt
      !
    ENDDO ! itemp
    !
    DEALLOCATE(sigmar_all_sum, STAT = ierr)
    IF (ierr /= 0) CALL errore('selfen_elec_print', 'Error deallocating sigmar_all_sum', 1)
    DEALLOCATE(xkf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('selfen_elec_print', 'Error deallocating xkf_all', 1)
    DEALLOCATE(etf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('selfen_elec_print', 'Error deallocating etf_all', 1)
    !
    102 FORMAT(5x, 'E( ', i3, ' )=', f12.6, ' eV   Re[Sigma]=', E22.14, ' meV Im[Sigma]=', &
               E22.14, ' meV     Z=', E22.14, ' lam=', E22.14)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE selfen_elec_print
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE selfen_phon_q(iqq, iq, totq)
    !-----------------------------------------------------------------------
    !!
    !! Compute the imaginary part of the phonon self energy due to electron-
    !! phonon interaction in the Migdal approximation and dynamically renormalize the 
    !! adiabatic frequencies. The imaginary part corresponds to
    !! the phonon linewidth (half width). The phonon frequency is taken into
    !! account in the energy selection rule.
    !!
    !! Use matrix elements, electronic eigenvalues and phonon frequencies
    !! from ep-wannier interpolation.  This routine is similar to the one above
    !! but it is ONLY called from within ephwann_shuffle and calculates
    !! the selfenergy for one phonon at a time.  Much smaller footprint on
    !! the disk
    !!
    !! RM 24/02/2014
    !! redefined the size of coskkq, vkk, vkq within the fermi windwow
    !! cleaned up the subroutine
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE modes,         ONLY : nmodes
    USE input,         ONLY : nbndsub, fsthick, efermi_read, fermi_energy,  &
                              nstemp, ngaussw, degaussw, shortrange,        &
                              nsmear, delta_smear, eps_acoustic, specfun_ph, &
                              delta_approx, vme, lfast_kmesh, isk_dummy, &
                              wmax_specfun, wmin_specfun, nw_specfun, &
                              phonselfen
    USE pwcom,         ONLY : nelec, ef
    USE global_var,    ONLY : epf17, ibndmin, etf, wkf, xqf, wqf, nkqf,  &
                              nkf, wf, xqf, lambda_all, lambda_v_all,    &
                              vmef, gamma_all, gamma_v_all, efnew, nbndfst, &
                              gtemp, nktotf, adapt_smearing, pi_0, gammai_all, &
                              pir_all
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, ryd2ev, one, two, zero, eps4,&
                              eps6, eps8, pi, ci, cone
    USE utilities,     ONLY : fermi_dirac
    USE klist,         ONLY : degauss, ngauss 
   !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point index from the selecq
    INTEGER, INTENT(in) :: iq
    !! Current q-point index
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points in selecq.fmt
    !
    ! Local variables
    !
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: jmode
    !! Counter on mode
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: ismear
    !! Number of smearing values for the Gaussian function
    INTEGER :: n
    !! Counter on number of mode degeneracies
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: iw
    !! Counter on energy for the spectral function
    !
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared in Ry^2
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: dosef
    !! Density of state N(Ef)
    REAL(KIND = DP) :: ekk
    !! Eigen energy at k on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy at k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: etmp1
    !! Temporary variable to store etmp1 = ekq - ekk
    REAL(KIND = DP) :: etmp2
    !! Temporary variable to store etmp2 = ekq - ekk - ww
    REAL(KIND = DP) :: wgkk
    !! Fermi-Dirac occupation factor $f_{nk}(T)$
    REAL(KIND = DP) :: wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
    REAL(KIND = DP) :: weight
    !! Imaginary or real part of the phonon self-energy factor
    !!$$ \pi N_q \Im(\frac{f_{nk}(T) - f_{mk+q(T)}}{\varepsilon_{nk}-\varepsilon_{mk+q}-\omega_{q\nu}+i\delta}) $$
    !! In practice the imaginary is performed with a delta Dirac
    REAL(KIND = DP) :: weight0
    !! static part of the phonon self-energy factor
    REAL(KIND = DP) :: w0g1
    !! Dirac delta at k for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: w0g2
    !! Dirac delta at k+q for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: degaussw0
    !! degaussw0 = (ismear-1) * delta_smear + degaussw
    REAL(KIND = DP) :: inv_degaussw0
    !! Inverse degaussw0 for efficiency reasons
    REAL(KIND = DP) :: inv_degauss
    !! inverse smearing from the scf calculation
    REAL(KIND = DP) :: ww(nw_specfun)
    !! Current frequency
    REAL(KIND = DP) :: dw
    !! frequency increment
    REAL(KIND = DP) :: lambda_tot
    !! Integrated lambda function
    REAL(KIND = DP) :: lambda_tr_tot
    !! Integrated transport lambda function
    REAL(KIND = DP) :: tmp1
    !! Temporary value of lambda for av.
    REAL(KIND = DP) :: tmp2
    !! Temporary value of lambda_v for av.
    REAL(KIND = DP) :: tmp3
    !! Temporary value of lambda_v for av.
    REAL(KIND = DP) :: tmp4
    !! Temporary value of lambda_v for av.
    REAL(KIND = DP) :: DDOT
    !! Dot product function
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! Function to compute the Density of States at the Fermi level
    REAL(KIND = DP), EXTERNAL :: efermig
    !! Return the fermi energy
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    REAL(KIND = DP) :: fact
    !! difference between occupations: fnk - fmk+q
    REAL(KIND = DP) :: fact_scf
    !! difference between occupations: fnk - fmk+q, but calculated for the smearing type and value from the ph.x calculation
    REAL(KIND = DP) :: wq(nmodes)
    !! Phonon frequency on the fine grid
    REAL(KIND = DP) :: wq_dyn(nmodes)
    !! One-shot nonadiabatic phonon frequency on the fine grid (on-shell)
    REAL(KIND = DP) :: wq_dyn_sc(nmodes)
    !! self-consistent nonadiabatic phonon frequency
    REAL(KIND = DP) :: inv_wq(nmodes)
    !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
    REAL(KIND = DP) :: g2_tmp(nmodes)
    !! If the phonon frequency is too small discart g
    REAL(KIND = DP) :: gamma(nmodes)
    !! Gamma is the imaginary part of the phonon self-energy
    REAL(KIND = DP) :: gamma_v(nmodes)
    !! Gamma is the imaginary part of the phonon self-energy multiplied by (1-coskkq)
    REAL(KIND = DP) :: pi_r(nmodes)
    !! dynamical real part of the phonon self-energy
    REAL(KIND = DP) :: lambda_tmp(nmodes)
    !! Temporary value of lambda for av.
    REAL(KIND = DP) :: lambda_v_tmp(nmodes)
    !! Temporary value of lambda v for av.
    REAL(KIND = DP) :: gamma_tmp(nmodes)
    !! Temporary value of gamma for av.
    REAL(KIND = DP) :: gamma_v_tmp(nmodes)
    !! Temporary value of gamma v for av.
    REAL(KIND = DP) :: vkk(3, nbndfst)
    !! Electronic velocity $v_{nk}$
    REAL(KIND = DP) :: vkq(3, nbndfst)
    !! Electronic velocity $v_{nk+q}$
    REAL(KIND = DP) :: coskkq(nbndfst, nbndfst)
    !! $$(v_k \cdot v_{k+q}) / |v_k|^2$$
    !
    CHARACTER(LEN=256) :: ngauss_         
    !! smearing type used for the ph.x calculation, read from an xml file
    !
    !
    IF (adapt_smearing) CALL errore('selfen_phon_q', 'adapt_smearing cannot be used with phonon self-energy', 1)
    !
    IF (specfun_ph) THEN
      dw = (wmax_specfun - wmin_specfun) / DBLE(nw_specfun - 1)
      DO iw = 1, nw_specfun
        ww(iw) = wmin_specfun + DBLE(iw - 1) * dw
      ENDDO
      pir_all(:, :, :, :) = zero
      gammai_all(:, :, :, :) = zero
    ENDIF
    !
    !    
    DO itemp = 1, nstemp
      IF (iq == 1) THEN
        IF (phonselfen) THEN
          WRITE(stdout, '(/5x, a)') REPEAT('=',67)
          WRITE(stdout, '(5x, "Phonon (Imaginary) Self-Energy in the Migdal Approximation")')
          WRITE(stdout, '(5x, a/)') REPEAT('=',67)
          !
          IF (ngauss == 0) THEN
            ngauss_ = 'Gaussian'
          ELSEIF (ngauss == 1) THEN
            ngauss_ =  'Methfessel-Paxton'
          ELSEIF (ngauss == -1) THEN
            ngauss_ = 'Marzari-Vanderbilt'
          ELSEIF (ngauss == -99) THEN
            ngauss_ = 'Fermi-Dirac'
          ENDIF
        IF (fsthick < 1.d3 ) WRITE(stdout, '(/5x, a, f10.6, a)' ) &
                 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
          WRITE(stdout, '(/5x, a, f10.6, a)' ) 'Golden Rule strictly enforced with T = ', gtemp(itemp) * ryd2ev, ' eV'
          WRITE(stdout, '(/5x, a)' ) 'On- and Off-shell dynamical frequency renormalization will be performed using:'
          WRITE(stdout, '(/5x, a)' ) ' Omega^2 = omega^2 + 2omega(Re[Pi(Omega,T_low)]- Pi[0,T_high]).'
          WRITE(stdout, '(5x, a, a)') 'Smearing type for the static self-energy: ', trim(ngauss_)
          WRITE(stdout, '(5x, a, f10.6, a)') 'smearing value for the static self-energy: T_high = ', degauss,' Ry.'
        ENDIF
      ENDIF
      !
      !
      ! Now pre-treat phonon modes for efficiency
      ! Treat phonon frequency and Bose occupation
      wq(:) = zero
      wq_dyn(:) = zero
      wq_dyn_sc(:) = zero
      DO imode = 1, nmodes
        IF (lfast_kmesh) THEN
          wq(imode) = wf(imode, iqq)
        ELSE
          wq(imode) = wf(imode, iq)
        ENDIF
        IF (wq(imode) > eps_acoustic) THEN
          g2_tmp(imode) = one
          inv_wq(imode) = one / (two * wq(imode))
        ELSE
          g2_tmp(imode) = zero
          inv_wq(imode) = zero
        ENDIF
      ENDDO
      !
      inv_degauss = one / degauss
      !
      DO ismear = 1, nsmear
        !
        degaussw0 = (ismear - 1) * delta_smear + degaussw
        !
        ! SP: Multiplication is faster than division ==> Important if called a lot
        !     in inner loops
        inv_degaussw0 = one / degaussw0
        !
        ! Fermi level and corresponding DOS
        !
        IF (efermi_read) THEN
          ef0 = fermi_energy
        ELSEIF (nsmear > 1) THEN
          !
          ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw0, ngaussw, 0, isk_dummy)
          ! if some bands are skipped (nbndskip /= 0), nelec has already been
          ! recalculated in ephwann_shuffle
          !
        ELSE !SP: This is added for efficiency reason because the efermig routine is slow
          ef0 = efnew
        ENDIF
        !
        dosef = dos_ef(ngaussw, degaussw0, ef0, etf, wkf, nkqf, nbndsub)
        !  N(Ef) in the equation for lambda is the DOS per spin
        dosef = dosef / two
        !
        IF (iq == 1) THEN
          WRITE (stdout, 100) degaussw0 * ryd2ev, ngaussw
          WRITE (stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
        ENDIF
        !
        CALL start_clock('PH SELF-ENERGY')
        !
        fermicount = 0
        wgkk = zero
        w0g1 = zero
        gamma(:)   = zero
        gamma_v(:) = zero
        pi_0(:)  = zero
        pi_r(:) = zero
        !
        DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          coskkq = zero
          ! coskkq = (vk dot vkq) / |vk|^2  appears in Grimvall 8.20
          ! this is different from :   coskkq = (vk dot vkq) / |vk||vkq|
          ! In principle the only coskkq contributing to lambda_tr are both near the
          ! Fermi surface and the magnitudes will not differ greatly between vk and vkq
          ! we may implement the approximation to the angle between k and k+q
          ! vectors also listed in Grimvall
          !
          DO ibnd = 1, nbndfst
            DO jbnd = 1, nbndfst
              !
              ! vmef is in units of Ryd * bohr
              !
              vkk(:, ibnd) = REAL(vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
              vkq(:, jbnd) = REAL(vmef(:, ibndmin - 1 + jbnd, ibndmin - 1 + jbnd, ikq))
              IF (ABS(vkk(1, ibnd)**two + vkk(2, ibnd)**two + vkk(3, ibnd)**two) > eps4) &
                coskkq(ibnd, jbnd) = DDOT(3, vkk(:, ibnd), 1, vkq(:, jbnd), 1) / &
                                     DDOT(3, vkk(:, ibnd), 1, vkk(:, ibnd), 1)
            ENDDO
          ENDDO
          !
          ! Here we must have ef, not ef0, to be consistent with ephwann_shuffle
          IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
              (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
            !
            fermicount = fermicount + 1
            DO imode = 1, nmodes
              !
              DO ibnd = 1, nbndfst
                !
                !  the fermi occupation for k
                ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
                !
                IF (delta_approx) THEN
                  w0g1 = w0gauss(ekk * inv_degaussw0, 0) * inv_degaussw0
                ENDIF
                wgkk = fermi_dirac(ekk, gtemp(itemp))
                !
                DO jbnd = 1, nbndfst
                  !
                  !  the fermi occupation for k+q
                  ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
                  etmp1 = ekq - ekk
                  !
                  ! here we take into account the zero-point DSQRT(hbar/2M\omega)
                  ! with hbar = 1 and M already contained in the eigenmodes
                  ! g2 is Ry^2, wkf must already account for the spin factor
                  !
                  IF (shortrange .AND. (ABS(xqf(1, iq)) > eps8 .OR. ABS(xqf(2, iq)) > eps8 &
                                        .OR. ABS(xqf(3, iq)) > eps8)) THEN
                    ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary
                    !     number, in which case its square will be a negative number.
                    g2 = REAL((epf17(jbnd, ibnd, imode, ik)**two) * inv_wq(imode) * g2_tmp(imode))
                  ELSE
                    g2 = (ABS(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq(imode) * g2_tmp(imode)
                  ENDIF
                  !
                  wgkq = fermi_dirac(ekq, gtemp(itemp))
                  fact = wgkq - wgkk
                  !
                  IF (delta_approx) THEN
                    !
                    w0g2 = w0gauss(ekq * inv_degaussw0, 0) * inv_degaussw0
                    ! the expression below is positive-definite, but also an
                    ! approximation which neglects some fine features
                    weight = pi * wq(imode) * wkf(ikk) * w0g1 * w0g2
                    !
                  ELSE
                    !
                    ! = k-point weight * [f(E_k) - f(E_k+q)] / [E_k+q - E_k - w_q + id]
                    ! This is the imaginary part of the phonon self-energy, sans
                    ! the matrix elements [Eq. 4 in Comput. Phys. Commun. 209, 116 (2016)]
                    !
                    ! weight = wkf (ikk) * (wgkk - wgkq) * AIMAG(cone / (ekq - ekk - wq - ci * degaussw0))
                    !
                    ! SP: The expression below is the imag part of phonon self-energy,
                    ! sans matrix elements [Eq. 9 in Comput. Phys. Commun. 209, 116 (2016)]
                    !  = pi * k-point weight * [f(E_k) - f(E_k+q)] * delta[E_k+q - E_k - w_q]
                    !
                    weight = -pi * wkf(ikk) * fact * w0gauss((etmp1 - wq(imode)) * inv_degaussw0, 0) * inv_degaussw0
                    !
                  ENDIF
                  !
                  gamma(imode)   = gamma(imode)   + weight * g2
                  gamma_v(imode) = gamma_v(imode) + weight * g2 * (1.0d0 - coskkq(ibnd, jbnd))
                  !
                  ! calculate the static part of the phonon self-energy from the ph.x calculation
                  ! in the limit of a vanishing denominator, the static self-energy becomes a derivative 
                  ! (f(ekq) - f(ekk) / (ekq - ekk) -> -(df/de)_{e=ekk}
                  !
                  fact_scf = wgauss(-ekq * inv_degauss, ngauss) - wgauss(-ekk * inv_degauss, ngauss)
                  fact = wgkq - wgkk
                  !
                  IF (abs(etmp1) < 1.0D-5 ) THEN
                    weight0 = -wkf (ikk) * inv_degauss * w0gauss(ekk * inv_degauss, ngauss)
                  ELSE
                    weight0 = wkf(ikk) * fact_scf / etmp1
                  ENDIF
                  !
                  pi_0(imode)    =  pi_0(imode)   + weight0 * g2
                  ! calculate the real part of the dynamical phonon self-energy          
                  weight = wkf(ikk) * fact * REAL(cone / (etmp1 -  wq(imode) - ci * degaussw0))
                  pi_r(imode)    = pi_r(imode)   + weight * g2
                  !
                  IF (specfun_ph) THEN
                    DO iw = 1, nw_specfun
                      !
                      ! below we calculate the full phonon self-energy, its real and imaginary parts
                      ! 
                      etmp2 = etmp1 - ww(iw)
                      weight = wkf(ikk) * fact * REAL(cone / (etmp2 - ci * degaussw0))
                      !
                      pir_all(iw, imode, itemp, ismear) = pir_all(iw, imode, itemp, ismear) + weight * g2
                      !
                      ! Normal implementation
                      ! weight = wkf (ikk) * fact * AIMAG(cone / (etmp2 - ci * degaussw))
                      !
                      ! More stable:
                      ! Analytical im. part
                      weight = pi * wkf(ikk) * fact * w0gauss(etmp2 * inv_degaussw0, 0) * inv_degaussw0
                      !
                      gammai_all(iw, imode, itemp, ismear) = gammai_all(iw, imode, itemp, ismear) + weight * g2
                      !
                    ENDDO ! iw
                  ENDIF 
                ENDDO ! jbnd
              ENDDO   ! ibnd
            ENDDO ! loop on q-modes
          ENDIF ! endif fsthick
        ENDDO ! loop on k
        !
        CALL stop_clock('PH SELF-ENERGY')
        !
        ! collect contributions from all pools (sum over k-points)
        ! this finishes the integral over the BZ  (k)
        !
	CALL mp_sum(pi_r, inter_pool_comm)
	CALL mp_sum(pi_0, inter_pool_comm)
        CALL mp_sum(gamma, inter_pool_comm)
        CALL mp_sum(gamma_v, inter_pool_comm)
        CALL mp_sum(fermicount, inter_pool_comm)
        CALL mp_barrier(inter_pool_comm)
        !
        IF (phonselfen) THEN
          ! An average over degenerate phonon-mode is performed.
          DO imode = 1, nmodes
            ! calculate the one-shot dynamical correction to the adiabatic frequency.
            wq_dyn(imode) = wq(imode)**2 + 2*wq(imode) * (pi_r(imode) - pi_0(imode))
            wq_dyn(imode) =  SQRT(ABS(wq_dyn(imode))) * wq_dyn(imode) / ABS(wq_dyn(imode))
            !
            ! If there are any dynamical corrections to the DFPT frequency, calculate also the self-consistent (off-shell) correction 
            IF (abs(wq_dyn(imode) - wq(imode)) > 1.0D-5) THEN
              CALL omega_dyn_sc(iqq, iq, imode, degaussw0, gtemp(itemp), pi_0(imode), ef0, ef, wq(imode), wq_dyn_sc(imode))
            ELSE
              wq_dyn_sc(imode) = wq_dyn(imode)
            ENDIF
            ! 
            n = 0
            tmp1 = zero
            tmp2 = zero
            tmp3 = zero
            tmp4 = zero
            DO jmode = 1, nmodes
              IF (ABS(wq(imode) - wq(jmode)) < eps6) THEN
                n = n + 1
                IF (wq(jmode) > eps_acoustic) THEN
                  tmp1 =  tmp1 + gamma(jmode)   / pi / wq(imode)**two / dosef
                  tmp2 =  tmp2 + gamma_v(jmode) / pi / wq(imode)**two / dosef
                ENDIF
                tmp3 =  tmp3 + gamma(jmode)
                tmp4 =  tmp4 + gamma_v(jmode)
              ENDIF
            ENDDO ! jbnd
            lambda_tmp(imode)   = tmp1 / FLOAT(n)
            lambda_v_tmp(imode) = tmp2 / FLOAT(n)
            gamma_tmp(imode)    = tmp3 / FLOAT(n)
            gamma_v_tmp(imode)  = tmp4 / FLOAT(n)
          ENDDO
          lambda_all(:, iq, ismear, itemp)   = lambda_tmp(:)
          lambda_v_all(:, iq, ismear, itemp) = lambda_v_tmp(:)
          gamma_all(:, iq, ismear, itemp)    = gamma_tmp(:)
          gamma_v_all(:, iq, ismear, itemp)  = gamma_v_tmp(:)
          lambda_tot    = SUM(lambda_all(:, iq, ismear, itemp))
          lambda_tr_tot = SUM(lambda_v_all(:, iq, ismear, itemp))
          !
          WRITE(stdout, '(/5x, "ismear = ",i5," iq = ",i7," coord.: ", 3f9.5, " wt: ", f9.5, " Temp: ", f8.3, "K")') ismear, iq, &
                                                                        xqf(:, iq), wqf(iq), gtemp(itemp) * ryd2ev / kelvin2eV
          WRITE(stdout, '(5x, a)') REPEAT('-', 67)
          !
          DO imode = 1, nmodes
            !
            WRITE(stdout, 102) imode, lambda_all(imode, iq, ismear, itemp), &
                               ryd2mev * gamma_all(imode, iq, ismear, itemp), &
                               ryd2mev * wq(imode), ryd2mev * wq_dyn(imode), &
                               ryd2mev * wq_dyn_sc(imode)
            WRITE(stdout, 104) imode, lambda_v_all(imode, iq, ismear, itemp), &
                               ryd2mev * gamma_v_all(imode, iq, ismear, itemp),&
                               ryd2mev * wq(imode), ryd2mev * wq_dyn(imode), &
                               ryd2mev * wq_dyn_sc(imode)
          !
          ENDDO
          !
          WRITE(stdout, 103) lambda_tot
          WRITE(stdout, 105) lambda_tr_tot
          !
          WRITE(stdout, '(5x, a/)') REPEAT('-', 67)
          WRITE(stdout, '(/5x, a, i8, a, i8/)' ) 'Number of (k,k+q) pairs on the Fermi surface: ', fermicount, ' out of ', nktotf
          !
        ENDIF
      ENDDO !smears
      !
    ENDDO ! itemp
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    101 FORMAT(5x, 'DOS =', f10.6, ' states/spin/eV/Unit Cell at Ef=', f10.6, ' eV')
    102 FORMAT(5x, 'lambda___( ', i3, ' )=', f15.6, '   gamma___=', f15.6, ' meV', '  omega=', f12.4, ' meV', &
                   '   on_shell_Omega=', f12.4, ' meV'  , '  off_shell_Omega=', f12.4, ' meV')
    103 FORMAT(5x, 'lambda___( tot )=', f15.6)
    104 FORMAT(5x, 'lambda_tr( ',i3,' )=', f15.6, '   gamma_tr=', f15.6, ' meV', '  omega=', f12.4, ' meV' , &
                   '   on_shell_Omega=', f12.4, ' meV', '  off_shell_Omega=', f12.4, ' meV')
    105 FORMAT(5x, 'lambda_tr( tot )=', f15.6)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE selfen_phon_q
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE omega_dyn_sc(iqq, iq, imode, degaussw0, eptemp,pi_0, ef0, ef, wq, wq_dyn_sc)
    !-----------------------------------------------------------------------
    !! This routine uses the bisection method to calculate the off-shell dynamical correction to the DFPT frequency. 
    !! The self-consistent loop, determines its renormalization using the equation: Omega**2 = omega_DFPT**2 + 2*omega_DFPT*(Pi(Omega) - Pi(0))
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE modes,         ONLY : nmodes
    USE input,         ONLY : nbndsub, fsthick, shortrange, eps_acoustic
    USE global_var,    ONLY : epf17, ibndmin, etf, wkf, xqf, wqf, nkqf,  &
                              nkf, wf, xqf, nbndfst, nktotf
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE ep_constants,  ONLY : one, two, zero, eps4, eps6, eps8, pi, ci, cone, ryd2mev
    USE utilities,     ONLY : fermi_dirac
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point index from the selecq
    INTEGER, INTENT(in) :: iq
    !! Current q-point index
    ! Local variables
    !
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER, INTENT(in) :: imode
    !! Counter on mode
    INTEGER :: iter
    !! Counter on iterations 
    INTEGER :: max_iter
    !! maximum number of iterations
    !
    REAL(KIND = DP), INTENT(in) :: eptemp
    !! temperature from the epw input file
    REAL(KIND = DP), INTENT(in) :: ef0
    !! Fermi energy level
    REAL(KIND = DP), INTENT(in) ::  ef
    !! fermi energy consistent with the ephwann shuffle routine
    REAL(KIND = DP), INTENT(in) :: degaussw0
    !! smearing fromt he input file
    REAL(KIND = DP), INTENT(out) :: wq_dyn_sc
    !! off-shell dynamically renormalized frequency
    REAL(KIND = DP), INTENT(in) :: pi_0
    !! static part of the phonon self-energy
    REAL(KIND = DP), INTENT(in) :: wq
    !! adiabatic fine-grid frequency
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared in Ry^2
    REAL(KIND = DP) :: ekk
    !! Eigen energy at k on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy at k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: wgkk
    !! Fermi-Dirac occupation factor $f_{nk}(T)$
    REAL(KIND = DP) :: wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
    REAL(KIND = DP) :: weight
    !! Imaginary and real  part of the phonhon self-energy factor
    REAL(KIND = DP) :: inv_wq(nmodes)
    !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    REAL(KIND = DP) :: g2_tmp(nmodes)
    !! If the phonon frequency is too small discart g
    REAL(KIND = DP) ::  pi_a
    !! dynamical real part of the phonon self-energy for omega = a
    REAL(KIND = DP) ::  pi_b
    !! dynamical real part of the phonon self-energy  for omega = b
    REAL(KIND = DP) ::  pi_c
    !! dynamical real part of the phonon self-energy for omega = c
    REAL(KIND = DP) ::  b
    !! upper bound on the dynamically renormalized freq. for the bisection method
    REAL(KIND = DP) ::  a
    !! lower bound on the dynamically renormalized freq. for the bisection method
    REAL(KIND = DP) ::  c
    !! (a+b)/2
    REAL(KIND = DP) ::  ga
    !! dynamically renormalized frequency minus the argument of the ph.self-energy: Omega(a) - a 
    REAL(KIND = DP) ::  gb
    !! dynamically renormalized frequency minus the argument of the ph.self-energy: Omega(b) - b 
    REAL(KIND = DP) ::  gc
    !! dynamically renormalized frequency minus the argument of the ph.self-energy: Omega(c) -c
    REAL(KIND = DP) ::  tol
    !! the threshold for the bisection method convergence
    !   
    wq_dyn_sc = zero
    tol = 1.0D-5
    IF (wq > eps_acoustic) THEN
      g2_tmp(imode) = one
      inv_wq(imode) = one / (two * wq)
    ELSE
      g2_tmp(imode) = zero
      inv_wq(imode) = zero
    ENDIF
    !
    !
    max_iter = 200
    CALL start_clock ('SC-OMEGA')
    !the lower bound of the off-shell frequency
    a = wq * 0.5
    ! the upper bound for the off-shell frequency
    b = wq * 3.0
    !
    DO iter = 1, max_iter
      wgkk = zero
      wgkq = zero
      pi_a = zero
      pi_b = zero
      pi_c = zero
      c = (a + b) * 0.5
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
            (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
          !
          !
          DO ibnd = 1, nbndfst
            !
            !  the fermi occupation for k
            ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
            !
            wgkk = fermi_dirac(ekk, eptemp)
            !
            DO jbnd = 1, nbndfst
              !
              ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
              !
              ! here we take into account the zero-point DSQRT(hbar/2M\omega)
              ! with hbar = 1 and M already contained in the eigenmodes
              ! g2 is Ry^2, wkf must already account for the spin factor
              !
              IF (shortrange .AND. (ABS(xqf(1, iq)) > eps8 .OR. ABS(xqf(2, iq)) > eps8 &
                                       .OR. ABS(xqf(3, iq)) > eps8)) THEN
                ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary
                !     number, in which case its square will be a negative number.
                g2 = REAL((epf17(jbnd, ibnd, imode, ik)**two) * inv_wq(imode) * g2_tmp(imode))
              ELSE
                g2 = (ABS(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq(imode) * g2_tmp(imode)
              ENDIF
              !
              wgkq = fermi_dirac(ekq, eptemp)
              ! calculate the real part of the dynamical phonon self-energy          
              weight = wkf(ikk) * (wgkq - wgkk) * REAL(cone / (ekq - ekk - a - ci * degaussw0))
              pi_a    = pi_a   + weight * g2
              weight = wkf(ikk) * (wgkq - wgkk) * REAL(cone / (ekq - ekk - b - ci * degaussw0))
              pi_b    = pi_b   + weight * g2
              weight = wkf(ikk) * (wgkq - wgkk) * REAL(cone / (ekq - ekk - c - ci * degaussw0))
              pi_c    = pi_c   + weight * g2
              !
            ENDDO ! jbnd
          ENDDO   ! ibnd
        ENDIF ! endif fsthick
      ENDDO ! loop on k
      !
      CALL mp_sum(pi_a, inter_pool_comm)
      CALL mp_sum(pi_b, inter_pool_comm)
      CALL mp_sum(pi_c, inter_pool_comm)
      !Update the freq and calculate f(omega) and check if f(omega) =  Omega(omega) - omega = 0
      ga = wq**two + two * wq * (pi_a - pi_0)
      ga =  SQRT(ABS(ga)) * ga / ABS(ga) - a
      gb = wq**two + two * wq * (pi_b - pi_0)
      gb =  SQRT(ABS(gb)) * gb / ABS(gb) - b
      gc = wq**two + two * wq * (pi_c - pi_0)
      gc =  SQRT(ABS(gc)) * gc / ABS(gc)- c   
      IF (ABS(gc) < tol .OR. (b - a) * 0.5 < tol )  THEN
        wq_dyn_sc = c 
        EXIT
      ELSEIF (( ga * gc) < 0) THEN
        b = c 
      ELSE
        a = c
      ENDIF 
    ENDDO
    CALL stop_clock ('SC-OMEGA')
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE omega_dyn_sc
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE selfen_pl_q(iqq, iq, totq, first_cycle)
    !-----------------------------------------------------------------------
    !!
    !!  Compute the imaginary part of the electron self energy due to electron-
    !!  plasmon interaction.
    !!
    !!  The coupling coefficients have been evaluated analytically employing a
    !!  Lindhard function model for the dielectric function contribution due to
    !!  the extrinsic carriers.
    !!
    !!  There are 3 parameters that the users should provide in the input:
    !!    - DOS effective mass;
    !!    - carrier concentration (Only for doped semiconductors, it shouldn't be used for insulators);
    !!    - epsilon_infinity (e.g, from exp. or from RPA).
    !!
    !!  F. Caruso and S. Ponce - 2017
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_var,        ONLY : linewidth_elself
    USE input,         ONLY : nbndsub, fsthick, ngaussw, efermi_read, &
                              fermi_energy, degaussw, nel, meff, epsiheg, &
                              restart, restart_step, nstemp
    USE pwcom,         ONLY : ef
    USE global_var,    ONLY : etf, ibndmin, nkqf, xqf, vmef, adapt_smearing, &
                              nkf, wqf, xkf, nkqtotf, efnew, nbndfst, nktotf,  &
                              gtemp, sigmar_all, sigmai_all, zi_all, lower_bnd
    USE ep_constants,  ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, ci, eps6, eps8
    USE ep_constants,  ONLY : pi
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE cell_base,     ONLY : omega, alat, bg
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE io_selfen,     ONLY : selfen_el_write
    USE parallelism,   ONLY : poolgather2
    USE utilities,     ONLY : fermi_dirac
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: first_cycle
    !! Use to determine weather this is the first cycle after restart
    INTEGER, INTENT(in) :: iqq
    !! Q-index from the selected q
    INTEGER, INTENT(in) :: iq
    !! Q-index from the global q
    INTEGER, INTENT(in) :: totq
    !! Number of q-points in selecq window
    !
    ! Local varialbes
    CHARACTER(LEN = 20) :: tp
    !! String for temperatures
    CHARACTER(LEN = 256) :: fileselfen
    !! File name for self energy
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
    INTEGER :: n
    !! Integer for the degenerate average over eigenstates
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
    REAL(KIND = DP) :: etmp1
    !! Temporary variable to store etmp1 = ekk - (ekq - wq)
    REAL(KIND = DP) :: etmp2
    !! Temporary variable to strore etmp2 = ekk - (ekq + wq)
    REAL(KIND = DP) :: sq_etmp1
    !! Temporary variable to store etmp1^2
    REAL(KIND = DP) :: sq_etmp2
    !! Temporary variable to store etmp2^2
    REAL(KIND = DP) :: wgq
    !! Bose occupation factor $n_{q wpl}(T)$
    REAL(KIND = DP) :: wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
    REAL(KIND = DP) :: fact1
    !! Temporary variable to store $f_{mk+q}(T) + n_{q wpl}(T)$
    REAL(KIND = DP) :: fact2
    !! Temporary variable to store $1 - f_{mk+q}(T) + n_{q wpl}(T)$
    REAL(KIND = DP) :: weight
    !! SE factors
    REAL(KIND = DP) :: w0g1
    !! Dirac delta at k for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: w0g2
    !! Dirac delta at k+q for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse of degaussw defined for efficiency reasons
    REAL(KIND = DP) :: sq_degaussw
    !! Squared degaussw defined for efficiency reasons
    REAL(KIND = DP) :: g2_tmp
    !! Temporary variable defined for efficiency reasons
    REAL(KIND = DP) :: tmp1
    !! Temporary variable to store real part of Sigma for the degenerate average
    REAL(KIND = DP) :: tmp2
    !! Temporary variable to store imag part of Sigma for the degenerate average
    REAL(KIND = DP) :: tmp3
    !! Temporary variable to store Z for the degenerate average
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
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    REAL(KIND = DP) :: q(3)
    !! The q-point in cartesian unit.
    REAL(KIND = DP) :: sigmar_tmp(nbndfst)
    !! Temporary array to store the real-part of Sigma
    REAL(KIND = DP) :: sigmai_tmp(nbndfst)
    !! Temporary array to store the imag-part of Sigma
    REAL(KIND = DP) :: zi_tmp(nbndfst)
    !! Temporary array to store the Z
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Collect eigenenergies from all pools in parallel case
    !
    IF (adapt_smearing) CALL errore('selfen_pl_q', 'adapt_smearing cannot be used with plasmon self-energy', 1)
    !
    !
    DO itemp = 1, nstemp
      ! SP: Define the inverse so that we can efficiently multiply instead of dividing
      inv_degaussw = one / degaussw
      sq_degaussw  = degaussw * degaussw
      !
      IF (iqq == 1) THEN
        !
        WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
        WRITE(stdout, '(5x, "Electron-plasmon Self-Energy in the Migdal Approximation")')
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
        WRITE (stdout, 100) degaussw * ryd2ev, ngaussw
        WRITE (stdout,'(a)') ' '
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
        sigmar_all(:, 1:lower_bnd - 1, :) = zero
        sigmar_all(:, lower_bnd + nkf:nktotf, :) = zero
        sigmai_all(:, 1:lower_bnd - 1, :) = zero
        sigmai_all(:, lower_bnd + nkf:nktotf, :) = zero
        zi_all(:, 1:lower_bnd - 1, :) = zero
        zi_all(:, lower_bnd + nkf:nktotf, :) = zero
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
          ! Bose occupation factor
          wgq = fermi_dirac(wq, gtemp(itemp))
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
            ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
            ! (but in this case they are the same)
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
                  wgkq = fermi_dirac(ekq, gtemp(itemp))
                  !
                  ! Computation of the dipole
                  IF (ibnd == jbnd) THEN
                    IF (qnorm > eps8) THEN
                      dipole = one / sq_qin
                    ELSE
                      dipole = zero
                    ENDIF
                  ELSE
                    IF (ABS(ekk - ekk1) > eps8) THEN
                      ! TODO: Check the expression to confirm that division by 2 is correct.
                      dipole = REAL(      vmef(1, ibndmin - 1 + jbnd, ibndmin - 1 + ibnd, ikk) / 2.d0 *  &
                                    CONJG(vmef(1, ibndmin - 1 + jbnd, ibndmin - 1 + ibnd, ikk) / 2.d0) / &
                                    ((ekk1 - ekk)**two + sq_degaussw))
                    ELSE
                      dipole = zero
                    ENDIF
                  ENDIF
                  !
                  IF (ABS(dipole * sq_qin) > 1.d0) THEN
                    dipole = one / sq_qin
                  ENDIF
                  !
                  ! The q^-2 is cancelled by the q->0 limit of the dipole.
                  ! See e.g., pg. 258 of Grosso Parravicini.
                  ! electron-plasmon scattering matrix elements squared
                  g2 = dipole * g2_tmp
                  !
                  fact1 =       wgkq + wgq
                  fact2 = one - wgkq + wgq
                  etmp1 = ekk - (ekq - wq)
                  etmp2 = ekk - (ekq + wq)
                  sq_etmp1 = etmp1 * etmp1
                  sq_etmp2 = etmp2 * etmp2
                  !
                  weight = wqf(iq) * REAL(fact1 / (etmp1 - ci * degaussw) + fact2 / (etmp2 - ci * degaussw))
                  !
                  ! \Re\Sigma [Eq. 1 in PRB 94, 115208 (2016)]
                  sigmar_all(ibnd, ik + lower_bnd - 1, itemp) = sigmar_all(ibnd, ik + lower_bnd - 1, itemp) + g2 * weight
                  !
                  ! Delta implementation
                  w0g1 = w0gauss(etmp1 * inv_degaussw, 0) * inv_degaussw
                  w0g2 = w0gauss(etmp2 * inv_degaussw, 0) * inv_degaussw
                  !
                  weight = pi * wqf(iq) * (fact1 * w0g1 + fact2 * w0g2)
                  !
                  ! \Im\Sigma using delta approx. [Eq. 1 in PRB 94, 115208 (2016)]
                  sigmai_all(ibnd, ik + lower_bnd - 1, itemp) = sigmai_all(ibnd, ik + lower_bnd - 1, itemp) + g2 * weight
                  !
                  ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
                  !
                  weight = wqf(iq) * &
                          (fact1 * (sq_etmp1 - sq_degaussw) / (sq_etmp1 + sq_degaussw)**two +  &
                           fact2 * (sq_etmp2 - sq_degaussw) / (sq_etmp2 + sq_degaussw)**two)
                  !
                  zi_all(ibnd, ik + lower_bnd - 1, itemp) = zi_all(ibnd, ik + lower_bnd - 1, itemp) + g2 * weight
                  !
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
            CALL mp_sum(sigmar_all, inter_pool_comm)
            CALL mp_sum(sigmai_all, inter_pool_comm)
            CALL mp_sum(zi_all, inter_pool_comm)
            CALL mp_sum(fermicount, inter_pool_comm)
            CALL mp_barrier(inter_pool_comm)
            CALL selfen_el_write(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all)
          ENDIF
        ENDIF
      ENDIF ! in case of restart, do not do the first one
      !
    ENDDO ! itemp
    !
    ! The k points are distributed among pools: here we collect them
    !
    IF (iqq == totq) THEN
      !
      ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_pl_q', 'Error allocating xkf_all', 1)
      ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_pl_q', 'Error allocating etf_all', 1)
      xkf_all(:, :) = zero
      etf_all(:, :) = zero
      !
#if defined(__MPI)
      !
      ! note that poolgather2 works with the doubled grid (k and k+q)
      !
      CALL poolgather2(3, nkqtotf, nkqf, xkf, xkf_all)
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
      CALL mp_sum(sigmar_all, inter_pool_comm)
      CALL mp_sum(sigmai_all, inter_pool_comm)
      CALL mp_sum(zi_all, inter_pool_comm)
      CALL mp_sum(fermicount, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
#else
      !
      xkf_all = xkf
      etf_all = etf
      !
#endif
      !
      DO itemp = 1, nstemp
        !
        ! Average over degenerate eigenstates:
        WRITE(stdout, '(5x, "Average over degenerate eigenstates is performed")')
        WRITE(stdout, '(5x, a, f8.3, a)') "Temperature: ", gtemp(itemp) * ryd2ev / kelvin2eV, "K"
        !
        DO ik = 1, nktotf
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          DO ibnd = 1, nbndfst
            ekk = etf_all(ibndmin - 1 + ibnd, ikk)
            n = 0
            tmp1 = zero
            tmp2 = zero
            tmp3 = zero
            DO jbnd = 1, nbndfst
              ekk1 = etf_all(ibndmin - 1 + jbnd, ikk)
              IF (ABS(ekk1 - ekk) < eps6) THEN
                n = n + 1
                tmp1 = tmp1 + sigmar_all(jbnd, ik, itemp)
                tmp2 = tmp2 + sigmai_all(jbnd, ik, itemp)
                tmp3 = tmp3 + zi_all(jbnd, ik, itemp)
              ENDIF
              !
            ENDDO ! jbnd
            sigmar_tmp(ibnd) = tmp1 / FLOAT(n)
            sigmai_tmp(ibnd) = tmp2 / FLOAT(n)
            zi_tmp(ibnd)     = tmp3 / FLOAT(n)
            !
          ENDDO ! ibnd
          sigmar_all(:, ik, itemp) = sigmar_tmp(:)
          sigmai_all(:, ik, itemp) = sigmai_tmp(:)
          zi_all(:, ik, itemp)     = zi_tmp(:)
          !
        ENDDO ! nktotf
        !
        ! Output plasmon SE here after looping over all q-points (with their contributions summed in sigmar_all, etc.)
        !
        WRITE(stdout, '(5x, "WARNING: only the eigenstates within the Fermi window are meaningful")')
        !
        IF (mpime == ionode_id) THEN
          ! Write to file
          WRITE(tp, "(f8.3)") gtemp(itemp) * ryd2ev / kelvin2eV
          fileselfen = 'linewidth.plself.' // trim(adjustl(tp)) // 'K'
          OPEN(UNIT = linewidth_elself, FILE = fileselfen)
          WRITE(linewidth_elself, '(a)') '# Electron lifetime (meV)'
          WRITE(linewidth_elself, '(a)') '#      ik       ibnd                 E(ibnd)      Im(Sigma)(meV)'
          !
          DO ik = 1, nktotf
            !
            ikk = 2 * ik - 1
            ikq = ikk + 1
            !
            WRITE(stdout, '(/5x, "ik = ", i7, " coord.: ", 3f12.7)') ik, xkf_all(:, ikk)
            WRITE(stdout, '(5x, a)') REPEAT('-', 67)
            !
            DO ibnd = 1, nbndfst
              !
              ! note that ekk does not depend on q
              ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
              !
              ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
              zi_all(ibnd, ik, itemp) = one / (one + zi_all(ibnd, ik, itemp))
              !
              WRITE(stdout, 102) ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2mev * sigmar_all(ibnd, ik, itemp), &
                                 ryd2mev * sigmai_all(ibnd, ik, itemp), zi_all(ibnd, ik, itemp), &
                                 one / zi_all(ibnd, ik, itemp) - one
              WRITE(linewidth_elself, '(i9, 2x)', ADVANCE = 'no') ik
              WRITE(linewidth_elself, '(i9, 2x)', ADVANCE = 'no') ibndmin - 1 + ibnd
              WRITE(linewidth_elself, '(E22.14, 2x)', ADVANCE = 'no') ryd2ev * ekk
              WRITE(linewidth_elself, '(E22.14, 2x)') ryd2mev * sigmai_all(ibnd, ik, itemp)
              !
            ENDDO
            WRITE(stdout, '(5x, a/)') REPEAT('-', 67)
          ENDDO
          CLOSE(linewidth_elself)
        ENDIF
        !
      ENDDO ! itemp
      !
      DEALLOCATE(xkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_pl_q', 'Error deallocating xkf_all', 1)
      DEALLOCATE(etf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_pl_q', 'Error deallocating etf_all', 1)
      !
    ENDIF
    !
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    102 FORMAT(5x, 'E( ', i3, ' )=', f9.4, ' eV   Re[Sigma]=', f15.6, ' meV Im[Sigma]=', &
               f15.6, ' meV     Z=', f15.6, ' lam=', f15.6)
    !
    RETURN
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE selfen_pl_q
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE get_eps_mahan(q, rs, kf, eps0)
    !--------------------------------------------------------------------------
    !!
    !! Based on Eq. 5.166 of Mahan 2000.
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : eps6, eps10
    USE ep_constants,  ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(out) :: eps0
    !! Output dielectric function at zero frequency
    REAL(KIND = DP), INTENT(in) ::  q
    !! Norm of q wave-vector
    REAL(KIND = DP), INTENT(in) :: rs
    !! Spherical radius used to describe the density of an electron gas
    REAL(KIND = DP), INTENT(in) :: kf
    !! Fermi wave-vector
    !
    !Local variable
    REAL(KIND = DP) :: x
    !! Temporary variable for q / (2 kf)
    REAL(KIND = DP) :: alpha
    !!Temporary variable
    !
    alpha = (4.d0 / (9.d0 * pi))**(1.d0/3.d0)
    !
    IF (ABS(q) > eps10) THEN
      x    = q / (2.d0 * kf)
      eps0 = 1.d0 + (1.d0 - x**2.d0) / (2.d0 * x) * LOG(ABS((1.d0 + x)/(1.d0 - x)))
      eps0 = 1.d0 + alpha * rs * eps0 / (2.d0 * pi * (x**2.d0))
    ELSE
      x    = (q + eps6) / 2.d0 / kf
      eps0 = 1.d0 + (1.d0 - x**2.d0) / (2.d0 * x) * LOG(ABS((1.d0 + x) / (1.d0 - x)))
      eps0 = 1.d0 + alpha * rs / 2.d0 / pi / x**2.d0 * eps0
    ENDIF
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE get_eps_mahan
    !--------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE nesting_fn_q(iqq, iq)
    !-----------------------------------------------------------------------
    !!
    !! Compute the imaginary part of the phonon self energy due to electron-
    !! phonon interaction in the Migdal approximation. This corresponds to
    !! the phonon linewidth (half width). The phonon frequency is taken into
    !! account in the energy selection rule.
    !!
    !! Use matrix elements, electronic eigenvalues and phonon frequencies
    !! from ep-wannier interpolation.
    !!
    !-----------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE input,     ONLY : nbndsub, fsthick, ngaussw, degaussw, &
                          nsmear, delta_smear, efermi_read, fermi_energy
    USE pwcom,     ONLY : ef
    USE global_var,ONLY : ibndmin, etf, wkf, xqf, wqf, nkqf, nktotf, &
                          nkf, xqf, nbndfst, efnew
    USE ep_constants,  ONLY : ryd2ev, zero, one, two
    USE mp,        ONLY : mp_barrier, mp_sum
    USE mp_global, ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! Current q-point index from selecq
    INTEGER, INTENT(in) :: iq
    !! Current q-point index
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
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: ismear
    !! Counter on smearing values
    !
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: weight
    !! Imaginary part of the phonhon self-energy factor, sans e-ph matrix elements
    REAL(KIND = DP) :: dosef
    !! Density of state N(Ef)
    REAL(KIND = DP) :: w0g1
    !! Dirac delta at k for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: w0g2
    !! Dirac delta at k+q for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: degaussw0
    !! degaussw0 = (ismear-1) * delta_smear + degaussw
    REAL(KIND = DP) :: inv_degaussw0
    !! Inverse degaussw0 for efficiency reasons
    REAL(KIND = DP) :: gamma
    !! Nesting function
    REAL(KIND = DP) :: dos_ef
    !! Function returning the density of states at the Fermi level
    REAL(KIND = DP) :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    !
    IF (iqq == 1) THEN
      WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
      WRITE(stdout, '(5x, "Nesting Function in the double delta approx")')
      WRITE(stdout, '(5x, a/)') REPEAT('=', 67)
      !
      IF (fsthick < 1.d3) WRITE(stdout, '(/5x, a, f10.6, a)' ) &
        'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
    ENDIF
    !
    ! SP: The Gamma function needs to be put to 0 for each q
    gamma = zero
    !
    ! Here we loop on smearing values
    DO ismear = 1, nsmear
      !
      degaussw0 = (ismear - 1) * delta_smear + degaussw
      inv_degaussw0 = one / degaussw0
      !
      ! Fermi level and corresponding DOS
      !
      !   Note that the weights of k+q points must be set to zero here
      !   no spin-polarized calculation here
      IF (efermi_read) THEN
        ef0 = fermi_energy
      ELSE
        ef0 = efnew
      ENDIF
      !
      dosef = dos_ef(ngaussw, degaussw0, ef0, etf, wkf, nkqf, nbndsub)
      !  N(Ef) in the equation for lambda is the DOS per spin
      dosef = dosef / two
      !
      IF (iqq == 1) THEN
        WRITE(stdout, 100) degaussw0 * ryd2ev, ngaussw
        WRITE(stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
      ENDIF
      !
      !
      CALL start_clock('nesting')
      !
      fermicount = 0
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
        IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
            (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) then
          !
          fermicount = fermicount + 1
          !
          DO ibnd = 1, nbndfst
            !
            ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
            w0g1 = w0gauss(ekk * inv_degaussw0, 0) * inv_degaussw0
            !
            DO jbnd = 1, nbndfst
              !
              ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
              w0g2 = w0gauss(ekq *inv_degaussw0, 0) * inv_degaussw0
              !
              ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
              ! This is the imaginary part of the phonon self-energy, sans the matrix elements
              !
              ! weight = wkf (ikk) * (wgkk - wgkq) * &
              !      aimag ( cone / ( ekq - ekk  - ci * degaussw ) )
              !
              ! the below expression is positive-definite, but also an approximation
              ! which neglects some fine features
              !
              weight = wkf(ikk) * w0g1 * w0g2
              !
              gamma  = gamma  + weight
              !
            ENDDO ! jbnd
          ENDDO ! ibnd
        ENDIF ! endif fsthick
      ENDDO ! loop on k
      !
      ! collect contributions from all pools (sum over k-points)
      ! this finishes the integral over the BZ  (k)
      !
      CALL mp_sum(gamma, inter_pool_comm)
      CALL mp_sum(fermicount, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
      WRITE(stdout, '(/5x, "iq = ",i5," coord.: ", 3f9.5, " wt: ", f9.5)') iq, xqf(:, iq) , wqf(iq)
      WRITE(stdout, '(5x, a)') REPEAT('-', 67)
      !
      WRITE(stdout, 102) gamma
      WRITE(stdout, '(5x,a/)') REPEAT('-', 67)
      !
      WRITE(stdout, '(/5x, a, i8, a, i8/)') &
        'Number of (k,k+q) pairs on the Fermi surface: ', fermicount, ' out of ', nktotf
      !
      CALL stop_clock('nesting')
    ENDDO !smears
    !
100 FORMAT(5x, 'Gaussian Broadening: ', f7.3,' eV, ngauss=', i4)
101 FORMAT(5x, 'DOS =', f10.6, ' states/spin/eV/Unit Cell at Ef=', f10.6, ' eV')
102 FORMAT(5x, 'Nesting function (q)=', E15.6, ' [Adimensional]')
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE nesting_fn_q
    !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  END MODULE selfen
  !-----------------------------------------------------------------------

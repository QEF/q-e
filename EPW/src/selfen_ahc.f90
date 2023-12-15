  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE selfen_ahc
  !----------------------------------------------------------------------
  !!
  !! This module contains the various self-energy routines
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !--------------------------------------------------------------------------
    SUBROUTINE selfen_elec_ahc_active(iqq, iq)
    !--------------------------------------------------------------------------
    !!
    !!  Compute the electron self-energy due to electron-phonon interaction
    !!  using the Allen-Heine-Cardona (AHC) theory.
    !!
    !!  Use e-ph matrix elements, Debye-Waller matrix, upper Fan matrix,
    !!  electronic eigenvalues, and phonon frequencies from Wannier
    !!  interpolation.
    !!
    !!  This routine computes the contribution from the coarse phonon iq
    !!  point to all k-points. The outer loop in ephwann_shuffle.f90 over the
    !!  coarse q-grid will loop over all iq points.
    !!  The contribution from each iq is summed at the end of this subroutine
    !!  to recover the per-ik electron self energy
    !!
    !!  In contrast to selfen_elec_q, here this routine computes both the real
    !!  and imaginary part of the self-energy. This can be done because the
    !!  Debye-Waller and upper Fan matrices are interpolated.
    !!  Also, the outer loop loops only over the coarse q points. This is done
    !!  because the sum of the upper Fan and Debye-Waller self-energy converge
    !!  much faster than the lowe Fan self-energy. Accordingly, the phonon
    !!  momentum interpolation of upper Fan matrix is not implemented.
    !!
    !!  Implemented by Jae-Mo Lihm.
    !!
    !!  [1] J.-M. Lihm and C.-H. Park, PRB 101, 121102(R) (2020)
    !!  [2] J.-M. Lihm and C.-H. Park, PRX 11, 041053 (2021)
    !!
    !! TODO: restart
    !!
    !--------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE constants_epw, ONLY : ryd2mev, ryd2ev, one, two, zero, ci, eps6, eps8
    USE modes,         ONLY : nmodes
    USE epwcom,        ONLY : nbndsub, nstemp, ngaussw, degaussw, &
                              eps_acustic, efermi_read, fermi_energy, &
                              ahc_win_min, ahc_win_max, wfpt
    USE elph2,         ONLY : etf, ibndmin, eta, nbndfst, gtemp, &
                              nkf, epf17, wf, wqf, adapt_smearing, &
                              sigma_ahc_act, efnew, nktotf, lower_bnd, &
                              sigma_ahc_ldw, dwf17, xqf
    USE io_selfen,     ONLY : selfen_el_write
    !
    IMPLICIT NONE
    !
    ! LOGICAL, INTENT(inout) :: first_cycle
    ! !! Use to determine weather this is the first cycle after restart
    INTEGER, INTENT(in) :: iqq
    !! Q-point index from selecq.fmt window
    INTEGER, INTENT(in) :: iq
    !! Q-point index from full grid
    ! INTEGER, INTENT(in) :: totq
    ! !! Total number of q-points from the selecq.fmt grid.
    !
    ! Local variables
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ik_gl
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
    !
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
    REAL(KIND = DP) :: wgkq
    !! Fermi-Dirac occupation factor $f_{mk+q}(T)$
    REAL(KIND = DP) :: fact1
    !! Temporary variable to store $f_{mk+q}(T) + n_{q\nu}(T)$
    REAL(KIND = DP) :: fact2
    !! Temporary variable to store $1 - f_{mk+q}(T) + n_{q\nu}(T)$
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse of degaussw define for efficiency reasons
    REAL(KIND = DP) :: inv_eptemp
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP) :: eta_tmp
    !! Temporary variable eta2
    REAL(KIND = DP) :: sq_eta_tmp
    !! Temporary eta2^2
    REAL(KIND = DP) :: inv_eta_tmp
    !! Temporary varialbe inv_eta
    REAL(KIND = DP) :: coeff
    !! Self-energy coefficient
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
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
    COMPLEX(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared in Ry^2
    COMPLEX(KIND = DP) :: weight
    !! Self-energy factor
    !!$$ (\frac{f_{mk+q}(T) + n_{q\nu}(T)}{\varepsilon_{nk} - \varepsilon_{mk+q} + \omega_{q\nu} - i\delta}) $$
    !!$$ + (\frac{1 - f_{mk+q}(T) + n_{q\nu}(T)}{\varepsilon_{nk} - \varepsilon_{mk+q} - \omega_{q\nu} - i\delta}) $$
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
      !
      WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
      WRITE(stdout, '(5x, "Electron Self-Energy using the Allen-Heine-Cardona Theory")')
      WRITE(stdout, '(5x, " - Active part")')
      WRITE(stdout, '(5x, a)') REPEAT('=', 67)
      !
      WRITE(stdout, '(5x, a, f10.3, a, f10.3, a)' ) &
        'AHC active state window: ', ahc_win_min, ' eV to ', ahc_win_max, ' eV'
      !
      WRITE(stdout, '(5x, a, f12.6, a)') 'Fermi energy: ', ef0 * ryd2ev, ' eV'
      !
      IF (adapt_smearing) THEN
        WRITE(stdout, '(5x, a)') 'Use adaptive smearing'
      ELSE
        WRITE(stdout, '(5x, a, f10.6, a, i4)') &
          'Gaussian Broadening: ', degaussw * ryd2ev, ' eV, ngauss=', ngaussw
      ENDIF
      !
      WRITE(stdout, '(a)') ''
      !
    ENDIF ! iqq == 1
    !
    ! SP: Define the inverse so that we can efficiently multiply instead of dividing
    inv_degaussw = one / degaussw
    !
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
    DO itemp = 1, nstemp ! loop over temperatures
      inv_eptemp = one / gtemp(itemp)
      !
      ! Now pre-treat phonon modes for efficiency
      ! Treat phonon frequency and Bose occupation
      wq(:) = zero
      DO imode = 1, nmodes
        wq(imode) = wf(imode, iq)
        IF (wq(imode) > eps_acustic) THEN
          g2_tmp(imode) = one
          wgq(imode)    = wgauss(-wq(imode) * inv_eptemp, -99)
          wgq(imode)    = wgq(imode) / (one - two * wgq(imode))
          inv_wq(imode) = one / (two * wq(imode))
        ELSE
          g2_tmp(imode) = zero
          wgq(imode)    = zero
          inv_wq(imode) = zero
        ENDIF
      ENDDO
      !
      ! IF (restart) THEN
      !   ! Make everythin 0 except the range of k-points we are working on
      !   sigmar_all(:, 1:lower_bnd - 1) = zero
      !   sigmar_all(:, lower_bnd + nkf:nktotf) = zero
      !   sigmai_all(:, 1:lower_bnd - 1) = zero
      !   sigmai_all(:, lower_bnd + nkf:nktotf) = zero
      !   zi_all(:, 1:lower_bnd - 1) = zero
      !   zi_all(:, lower_bnd + nkf:nktotf) = zero
      !   !
      ! ENDIF
      !
      ! ! In the case of a restart do not add the first step
      ! IF (first_cycle) THEN
      !   first_cycle = .FALSE.
      ! ELSE
        !
        ! loop over all k points of the fine mesh
        !
        DO ik = 1, nkf
          !
          ik_gl = ik + lower_bnd - 1
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
          ! (but in this case they are the same)
          !
          ! fsthick is not used for AHC
          !
          DO imode = 1, nmodes
            DO ibnd = 1, nbndfst
              !
              ! the energy of the electron at k (relative to Ef)
              ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
              !
              eta_tmp     = eta2(ibnd, imode, ik)
              sq_eta_tmp  = eta_tmp**two
              inv_eta_tmp = inv_eta(ibnd, imode, ik)
              !
              DO jbnd = 1, nbndsub
                !
                ! the energy of the electron at k+q
                ekq = etf(jbnd, ikq)
                !
                ! Skip active states outside the ahc window
                IF (ekq < ahc_win_min / ryd2ev .OR. ekq > ahc_win_max / ryd2ev) CYCLE
                !
                ! Set energy relative to Fermi energy
                ekq = ekq - ef0
                !
                ! the Fermi occupation at k+q
                wgkq = wgauss(-ekq * inv_eptemp, -99)
                !
                ! Skip coupling with oneself or between degenerate states at the same k point
                IF (ALL(xqf(:, iq) < eps8) .AND. ABS(ekq - ekk) < 2.d-5) CYCLE
                !
                g2 = CONJG(epf17(jbnd, ibnd, imode, ik)) * epf17(jbnd, ibnd, imode, ik) &
                   * inv_wq(imode) * g2_tmp(imode)
                !
                ! There is a sign error for wq in Eq. 9 of Comp. Phys. Comm. 181, 2140 (2010). - RM
                ! The sign was corrected according to Eq. (7.282) page 489 from Mahan's book
                ! (Many-Particle Physics, 3rd edition)
                !
                fact1 =       wgkq + wgq(imode)
                fact2 = one - wgkq + wgq(imode)
                etmp1 = ekk - (ekq - wq(imode))
                etmp2 = ekk - (ekq + wq(imode))
                !
                weight = wqf(iq) * (fact1 / (etmp1 - ci * eta_tmp) + fact2 / (etmp2 - ci * eta_tmp))
                !
                ! \Re\Sigma [Eq. 3 in Comput. Phys. Commun. 209, 116 (2016)]
                sigma_ahc_act(ibnd, ik_gl, itemp) = sigma_ahc_act(ibnd, ik_gl, itemp) + g2 * weight
                !
              ENDDO !jbnd
              !
            ENDDO !ibnd
            !
            ! Active space Debye-Waller term
            !
            ! coeff = (n + 0.5) / (2 * omega) / nqc
            coeff = (wgq(imode) + one / two) * inv_wq(imode) * wqf(iq)
            !
            DO ibnd = 1, nbndfst
              !
              ! Debye-Waller
              IF (wfpt) THEN
                sigma_ahc_ldw(ibnd, ik_gl, itemp) = sigma_ahc_ldw(ibnd, ik_gl, itemp) &
                + coeff * REAL(dwf17(ibnd, ibnd, imode, ik))
              ENDIF
              !
            ENDDO ! ibnd
            !
          ENDDO !imode
          !
        ENDDO ! end loop on k
      ENDDO ! itemp
      !
      ! ! Creation of a restart point
      ! IF (restart) THEN
      !   IF (MOD(iqq, restart_step) == 0) THEN
      !     WRITE(stdout, '(5x, a, i10)' ) 'Creation of a restart point at ', iqq
      !     CALL mp_sum(sigmar_all, inter_pool_comm)
      !     CALL mp_sum(sigmai_all, inter_pool_comm)
      !     CALL mp_sum(zi_all, inter_pool_comm)
      !     CALL mp_barrier(inter_pool_comm)
      !     CALL selfen_el_write(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all)
      !   ENDIF
      ! ENDIF
    ! ENDIF ! in case of restart, do not do the first one
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE selfen_elec_ahc_active
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE selfen_elec_ahc_static(iq, nqc, wf_coarse)
    !--------------------------------------------------------------------------
    !!
    !!  Compute the electron self-energy due to electron-phonon interaction
    !!  using the Allen-Heine-Cardona (AHC) theory.
    !!
    !!  Use e-ph matrix elements, Debye-Waller matrix, upper Fan matrix,
    !!  electronic eigenvalues, and phonon frequencies from Wannier
    !!  interpolation.
    !!
    !!  This routine computes the contribution from the coarse phonon iq
    !!  point to all k-points. The outer loop in ephwann_shuffle.f90 over the
    !!  coarse q-grid will loop over all iq points.
    !!  The contribution from each iq is summed at the end of this subroutine
    !!  to recover the per-ik electron self energy
    !!
    !!  In contrast to selfen_elec_q, here this routine computes both the real
    !!  and imaginary part of the self-energy. This can be done because the
    !!  Debye-Waller and upper Fan matrices are interpolated.
    !!  Also, the outer loop loops only over the coarse q points. This is done
    !!  because the sum of the upper Fan and Debye-Waller self-energy converge
    !!  much faster than the lowe Fan self-energy. Accordingly, the phonon
    !!  momentum interpolation of upper Fan matrix is not implemented.
    !!
    !! TODO: restart
    !!
    !--------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : one, two, zero, czero, ryd2ev
    USE mp,            ONLY : mp_barrier, mp_sum
    USE io_global,     ONLY : stdout
    USE modes,         ONLY : nmodes
    USE epwcom,        ONLY : nbndsub, nstemp, eps_acustic, ahc_win_min, &
                              ahc_win_max, wfpt
    USE elph2,         ONLY : etf, ibndmin, nbndfst, nkf, epf17, gtemp, &
                              sigma_ahc_hdw, sigma_ahc_uf, lower_bnd, sthf17, dwf17
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! Q-point index from full grid
    INTEGER, INTENT(in) :: nqc
    !! Total number of coarse q-points
    REAL(KIND = DP), INTENT(in) :: wf_coarse(nmodes)
    !! Phonon frequency at coarse q-points
    !
    ! Local variables
    LOGICAL :: skip_mode(nmodes)
    !! If the phonon frequency is too small or negative, skip mode.
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ik_gl
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
    REAL(KIND = DP) :: coeff
    !! Self-energy coefficient
    REAL(KIND = DP) :: ekk
    !! Eigen energy at k on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy at k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: inv_eptemp
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    REAL(KIND = DP) :: wq(nmodes)
    !! Phonon frequency on the fine grid
    REAL(KIND = DP) :: inv_wq(nmodes)
    !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
    REAL(KIND = DP) :: wgq(nmodes)
    !! Bose occupation factor $n_{q\nu}(T)$
    COMPLEX(KIND = DP) :: sthmat(nbndfst, nmodes)
    !! Sternheimer matrix at each k-point
    !
    IF (.NOT. wfpt) CALL errore("selfen_elec_ahc_static", &
      "To compute static AHC self-energies, wfpt must be true", 1)
    !
    IF (iq == 1) THEN
      !
      WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
      WRITE(stdout, '(5x, "Electron Self-Energy using the Allen-Heine-Cardona Theory")')
      WRITE(stdout, '(5x, " - Static part")')
      WRITE(stdout, '(5x, a)') REPEAT('=', 67)
      !
      ! Write ahc window
      WRITE(stdout, '(5x, a, f10.3, a, f10.3, a)' ) &
        'AHC active state window: ', ahc_win_min, ' eV to ', ahc_win_max, ' eV'
      WRITE(stdout, *)
      !
    ENDIF
    !
    ! loop over all k points of the fine mesh
    !
    DO ik = 1, nkf
      !
      ik_gl = ik + lower_bnd - 1
      ikk = 2 * ik - 1
      ikq = ikk + 1
      !
      sthmat(:, :) = czero
      !
      ! Sternheimer contribution to the Upper Fan matrix.
      !
      DO imode = 1, nmodes
        DO ibnd = 1, nbndfst
          sthmat(ibnd, imode) = sthf17(ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, imode, ik)
        ENDDO
      ENDDO
      !
      ! Add inactive Wannier band contribution to upper Fan matrix.
      !
      DO imode = 1, nmodes
        DO ibnd = 1, nbndfst
          !
          ! the energy of the electron at k
          ekk = etf(ibndmin - 1 + ibnd, ikk)
          !
          ! JML : If ekk etf(ibndmin - 1 + ibnd, ikk) outside
          !
          DO jbnd = 1, nbndsub
            !
            ! the energy of the electron at k+q
            ekq = etf(jbnd, ikq)
            !
            ! Skip active states inside the ahc window
            IF (ahc_win_min / ryd2ev <= ekq .AND. ekq <= ahc_win_max / ryd2ev) CYCLE
            !
            sthmat(ibnd, imode) = sthmat(ibnd, imode) &
            + CONJG(epf17(jbnd, ibnd, imode, ik)) * epf17(jbnd, ibnd, imode, ik) / (ekk - ekq)
            !
          ENDDO ! jbnd
        ENDDO ! ibnd
      ENDDO ! imode
      !
      ! Calculate self-energy
      !
      DO itemp = 1, nstemp ! loop over temperatures
        !
        inv_eptemp = one / gtemp(itemp)
        !
        wq(:) = zero
        skip_mode(:) = .FALSE.
        !
        DO imode = 1, nmodes
          wq(imode) = wf_coarse(imode)
          IF (wq(imode) > eps_acustic) THEN
            skip_mode(imode) = .FALSE.
            wgq(imode)    = wgauss(-wq(imode) * inv_eptemp, -99)
            wgq(imode)    = wgq(imode) / (one - two * wgq(imode))
            inv_wq(imode) = one / (two * wq(imode))
          ELSE
            skip_mode(imode) = .TRUE.
            wgq(imode)    = zero
            inv_wq(imode) = zero
          ENDIF
          !
          IF (skip_mode(imode)) CYCLE
          !
          ! coeff = (n + 0.5) / (2 * omega) / nqc
          coeff = (wgq(imode) + one / two) * inv_wq(imode) / REAL(nqc, DP)
          !
          DO ibnd = 1, nbndfst
            !
            ekk = etf(ibndmin - 1 + ibnd, ikk)
            !
            ! Upper Fan
            sigma_ahc_uf(ibnd, ik_gl, itemp) = sigma_ahc_uf(ibnd, ik_gl, itemp) &
              + coeff * REAL((sthmat(ibnd, imode) + CONJG(sthmat(ibnd, imode))))
            !
            ! Debye-Waller
            sigma_ahc_hdw(ibnd, ik_gl, itemp) = sigma_ahc_hdw(ibnd, ik_gl, itemp) &
            + coeff * REAL(dwf17(ibnd, ibnd, imode, ik))
            !
          ENDDO ! ibnd
          !
        ENDDO ! imode
      ENDDO ! itemp
      !
    ENDDO ! ik
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE selfen_elec_ahc_static
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE compute_dw_truncated(first_q, ik, nrr, cuf, etf, etf_ks, cfac, &
                                    dims, uf_ph, epmatf, dwf)
    !--------------------------------------------------------------------------
    !! Compute DW term truncated to the active space at ik.
    !--------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : czero, ryd2ev, ci, eps8
    USE modes,         ONLY : nmodes
    USE epwcom,        ONLY : nbndsub, ahc_win_min, ahc_win_max
    USE elph2,         ONLY : ibndmin, nbndfst
    USE wan2bloch,     ONLY : dmewan2bloch, dwwan2blochp
    USE wfpt_mod,      ONLY : dwmatf_trunc
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: first_q
    !! Use to determine weather this is the first q point
    INTEGER, INTENT(in) :: ik
    !! kpoint index
    INTEGER, INTENT(in) :: nrr
    !! Number of real-space grid points
    INTEGER, INTENT(in) :: dims
    !! Is equal to the number of Wannier function if use_ws == .TRUE.
    !! Is equal to 1 otherwise.
    REAL(KIND = DP), INTENT(in) :: etf(nbndsub)
    !! Eigenenergies on the fine grid
    REAL(KIND = DP), INTENT(in) :: etf_ks(nbndsub)
    !! Kohn-Sham eigenvalues
    COMPLEX(KIND = DP), INTENT(IN) :: cuf(nbndsub, nbndsub)
    !! Rotation matrix U^\dagger, fine mesh
    COMPLEX(KIND = DP), INTENT(in) :: cfac(nrr, dims, dims)
    !! Exponential factor
    COMPLEX(KIND = DP), INTENT(IN) :: epmatf(nbndsub, nbndsub, nmodes)
    !! e-p matrix in smooth Bloch basis, fine mesh, phonon Catresian basis
    COMPLEX(KIND = DP), INTENT(in) :: uf_ph(nmodes, nmodes)
    !! phonon eigenmode
    COMPLEX(KIND = DP), INTENT(inout) :: dwf(nbndfst, nbndfst, 3, nmodes)
    !! Truncated Debye-Waller matrix element in electron Bloch basis, phonon Cartesian basis
    !
    ! Local variables
    INTEGER :: ib
    !! Band index
    INTEGER :: jb
    !! Band index
    INTEGER :: pb
    !! Band index
    INTEGER :: idir
    !! Cartesian direction
    INTEGER :: imode
    !! Mode index
    INTEGER :: ib_full
    !! Index on the full band space
    INTEGER :: jb_full
    !! Index on the full band space
    REAL(KIND = DP) :: ekq
    !! Eigen energy at k+q on the fine grid relative to the Fermi level
    COMPLEX(KIND = DP) :: pmef(3, nbndsub, nbndsub)
    !! momentum matrix elements on the fine mesh
    !
    ! Debye-Waller matrix element in phonon Cartesian basis is q-independent.
    ! So, compute it only once for the first q point.
    !
    IF (first_q) THEN
      !
      CALL dmewan2bloch(nbndsub, nrr, cuf, pmef, etf, etf_ks, cfac, dims, use_momentum = .TRUE.)
      !
      DO pb = 1, nbndsub
        !
        ! Skip states outside the AHC window
        ekq = etf(pb)
        IF (ekq < ahc_win_min / ryd2ev .OR. ekq > ahc_win_max / ryd2ev) CYCLE
        !
        DO ib = 1, nbndfst
          DO jb = 1, nbndfst
            ib_full = ib + ibndmin - 1
            jb_full = jb + ibndmin - 1
            !
            DO idir = 1, 3
              DO imode = 1, nmodes
                dwmatf_trunc(ib, jb, idir, imode, ik) &
                = dwmatf_trunc(ib, jb, idir, imode, ik) &
                + ci * CONJG(epmatf(pb, ib_full, imode)) * pmef(idir, pb, jb) &
                - ci * CONJG(pmef(idir, pb, ib)) * epmatf(pb, jb_full, imode)
              ENDDO
            ENDDO
            !
          ENDDO
        ENDDO
      ENDDO
    ENDIF ! first_q
    !
    ! Rotate dwmatf (phonon Cartesian) to dwmatf2 (phonon eigenmode)
    !
    CALL dwwan2blochp(nbndfst, dwmatf_trunc(:, :, :, :, ik), uf_ph, dwf, nmodes)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE compute_dw_truncated
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE ahc_run_static_wfpt(nrr_k, nrr_q, nrr_g, irvec_k, irvec_q, irvec_g, &
                                   ndegen_k, ndegen_q, ndegen_g, dims, dims2, rws, nrws)
    !--------------------------------------------------------------------------
    !!
    !! Allen-Heine-Cardona theory for the real part of electron self-energy
    !!
    !! The Debye-Waller and upper Fan terms are interpolated only for k-points,
    !! not for q-points. So, we calculate them here, after the q-point loop.
    !!
    ! -------------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE mp_world,         ONLY : world_comm
    USE mp,               ONLY : mp_bcast
    USE io_global,        ONLY : ionode_id, stdout, ionode
    USE io_files,         ONLY : prefix, diropn
    USE cell_base,        ONLY : at
    USE ions_base,        ONLY : amass, ityp
    USE modes,            ONLY : nmodes
    USE constants_epw,    ONLY : eps8, czero, twopi, ci, cone, zero
    USE epwcom,           ONLY : nbndsub, lifc, nqc1, nqc2, nqc3, use_ws, eig_read
    USE io_var,           ONLY : iusthwe, iudgwe, iuxqc
    USE elph2,            ONLY : dwmatwe, dgmatwe, sthmatwe, dwf17, sthf17, dgf17,   &
                                 epf17, nkf, nbndfst, xkf,etf, etf_ks, chw_ks, chw,  &
                                 ibndmin, ibndmax
    USE wfpt_mod,         ONLY : dwmatf2, dwmatf_trunc
    USE wan2bloch,        ONLY : dmewan2bloch, hamwan2bloch, dynwan2bloch,           &
                                 ephwan2blochp, ephwan2bloch, vmewan2bloch,          &
                                 dynifc2blochf, vmewan2blochp, rrwan2bloch,          &
                                 dwwan2blochp, sthwan2blochp, dgwan2blochp
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nrr_k
    !! Number of WS points for electrons
    INTEGER, INTENT(IN) :: nrr_q
    !! Number of WS points for phonons
    INTEGER, INTENT(IN) :: nrr_g
    !! Number of WS points for electron-phonons
    INTEGER, INTENT(IN) :: nrws
    !! Number of real-space Wigner-Seitz
    INTEGER, INTENT(IN) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(IN) :: dims2
    !! Dims is either nat if use_ws or 1 if not
    INTEGER, INTENT(IN) :: irvec_k(3, nrr_k)
    !! INTEGER components of the ir-th Wigner-Seitz grid point in the basis
    !! of the lattice vectors for electrons
    INTEGER, INTENT(IN) :: irvec_q(3, nrr_q)
    !! INTEGER components of the ir-th Wigner-Seitz grid point for phonons
    INTEGER, INTENT(IN) :: irvec_g(3, nrr_g)
    !! INTEGER components of the ir-th Wigner-Seitz grid point for electron-phonon
    INTEGER, INTENT(IN) :: ndegen_k(nrr_k, dims, dims)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, INTENT(IN) :: ndegen_q(nrr_q, dims2, dims2)
    !! Wigner-Seitz weights for the phonon grid that depend on atomic positions $R + \tau(nb) - \tau(na)$
    INTEGER, INTENT(IN) :: ndegen_g(dims, nrr_g, dims2)
    !! Wigner-Seitz weights for the electron-phonon grid that depend on
    !! atomic positions $R - \tau(na)$
    REAL(KIND = DP), INTENT(IN) :: rws(0:3, nrws)
    !! Real-space wigner-Seitz vectors
    !
    CHARACTER(LEN = 256) :: filint
    !! Name of the file to write/read
    LOGICAL :: exst
    !! If the file exist
    INTEGER :: imode, nu, mu
    !! Counter on modes
    INTEGER :: na
    !! Counter on atoms
    INTEGER :: ik, ikk, ikq
    !! Counter on k points
    INTEGER :: ir
    !! Counter on R vectors
    INTEGER :: iw, iw2
    !! Counter on Wannier functions
    INTEGER :: ibnd, jbnd
    !! Counter on bands
    INTEGER :: iq
    !! Counter on coarse q-point grid
    INTEGER :: nqc_ahc
    !! Coarse q-point grid size
    INTEGER :: lrsthmatw
    !! record length while reading file
    INTEGER :: lrdgmatw
    !! record length while reading file
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: xxq(3)
    !! Current q-point
    REAL(KIND = DP) :: xkk(3)
    !! Current k-point on the fine grid
    REAL(KIND = DP) :: xkq2(3)
    !! Current k+q point on the fine grid
    REAL(KIND = DP), ALLOCATABLE :: wf_coarse(:, :)
    !! Phonon frequency at coarse q-points
    REAL(KIND = DP), ALLOCATABLE :: w2(:)
    !! Interpolated phonon frequency
    REAL(KIND = DP), ALLOCATABLE :: xqc_loc(:, :)
    !! Coarse q-point grid
    REAL(KIND = DP), ALLOCATABLE :: irvec_r(:, :)
    !! Wigner-Size supercell vectors, store in real instead of integer
    REAL(KIND = DP), ALLOCATABLE :: rdotk(:)
    !! $r\cdot k$
    REAL(KIND = DP), ALLOCATABLE :: rdotk2(:)
    !! $r\cdot k$
    COMPLEX(KIND = DP), ALLOCATABLE :: cufkk(:, :)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP), ALLOCATABLE :: cufkq(:, :)
    !! Rotation matrix, fine mesh, points k+q
    COMPLEX(KIND = DP), ALLOCATABLE :: uf(:, :)
    !! Rotation matrix for phonons
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatwef(:, :, :, :)
    !! e-p matrix  in el wannier - fine Bloch phonon grid
    COMPLEX(KIND = DP), ALLOCATABLE :: dwmatf(:, :, :, :)
    !! Debye-Waller matrix in smooth Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: sthmatf(:, :, :, :)
    !! Sternheimer matrix in smooth Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: dgmatf(:, :, :)
    !! Hopping correction matrix in smooth Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: sthmatf2(:, :, :)
    !! Sternheimer matrix in smooth Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: dgmatf2(:, :, :)
    !! Hopping correction matrix in smooth Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: cfac(:, :, :)
    !! Used to store $e^{2\pi r \cdot k}$ exponential
    COMPLEX(KIND = DP), ALLOCATABLE :: cfacq(:, :, :)
    !! Used to store $e^{2\pi r \cdot k+q}$ exponential
    COMPLEX(KIND = DP), ALLOCATABLE :: rrf(:, :, :)
    !! Interpolation position matrix elements on the fine mesh (ipol, nbnd, nbnd)
    COMPLEX(KIND = DP), ALLOCATABLE :: tmp(:, :, :)
    !! Overlap k and k+q
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatf(:, :, :)
    !! e-p matrix  in smooth Bloch basis, fine mesh
    !
    CALL start_clock('ep-int-ahc')
    !
    nqc_ahc = nqc1 * nqc2 * nqc3
    !
    ALLOCATE(rrf(3, nbndsub, nbndsub), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating rrf', 1)
    rrf(:, :, :) = czero
    ALLOCATE(cfac(nrr_k, dims, dims), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating cfac', 1)
    ALLOCATE(cfacq(nrr_k, dims, dims), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating cfacq', 1)
    ALLOCATE(rdotk(nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating rdotk', 1)
    ALLOCATE(rdotk2(nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating rdotk2', 1)
    ! This is simply because dgemv take only real number (not integer)
    ALLOCATE(irvec_r(3, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating irvec_r', 1)
    irvec_r = REAL(irvec_k, KIND = DP)
    ALLOCATE(tmp(nbndfst, nbndfst, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating tmp', 1)
    ALLOCATE(epmatf(nbndsub, nbndsub, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating epmatf', 1)
    !
    ! Zeroing everything - initialization is important !
    cfac(:, :, :)  = czero
    cfacq(:, :, :) = czero
    rdotk(:)       = zero
    rdotk2(:)      = zero
    !
    ALLOCATE(epmatwef(nbndsub, nbndsub, nrr_k, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating epmatwef', 1)
    ALLOCATE(w2(nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating w2', 1)
    ALLOCATE(cufkk(nbndsub, nbndsub), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating cufkk', 1)
    ALLOCATE(cufkq(nbndsub, nbndsub), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating cufkq', 1)
    ALLOCATE(uf(nmodes, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating uf', 1)
    ALLOCATE(xqc_loc(3, nqc_ahc), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating xqc_loc', 1)
    ALLOCATE(wf_coarse(nmodes, nqc_ahc), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating wf_coarse', 1)
    ALLOCATE(dwmatf(nbndsub, nbndsub, 3, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating dwmatf', 1)
    ALLOCATE(sthmatwe(nbndsub, nbndsub, nrr_k, nmodes, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating sthmatwe', 1)
    ALLOCATE(dgmatwe(nbndsub, nbndsub, nrr_k, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating dgmatwe', 1)
    ALLOCATE(sthmatf(nbndsub, nbndsub, nmodes, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating sthmatf', 1)
    ALLOCATE(dgmatf(nbndsub, nbndsub, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating dgmatf', 1)
    ALLOCATE(sthmatf2(nbndsub, nbndsub, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating sthmatf2', 1)
    ALLOCATE(dgmatf2(nbndsub, nbndsub, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating dgmatf2', 1)
    ALLOCATE(sthf17(nbndfst, nbndfst, nmodes, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating sthf17', 1)
    ALLOCATE(dgf17(nbndfst, nbndfst, nmodes, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating dgf17', 1)
    !
    ! Read xqc from file
    !
    IF (ionode) THEN
      CALL diropn(iuxqc, 'xqc', 3 * nqc_ahc, exst)
      CALL davcio(xqc_loc, 3 * nqc_ahc, iuxqc, 1, -1)
      CLOSE(iuxqc)
      ! Transform to crystal coordinates
      CALL cryst_to_cart(nqc_ahc, xqc_loc, at, -1)
    ENDIF
    CALL mp_bcast(xqc_loc, ionode_id, world_comm)
    !
    ! Open the prefix.sthmatwe file
    IF (ionode) THEN
      lrsthmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
      filint = TRIM(prefix)//'.sthmatwe'
      CALL diropn(iusthwe, 'sthmatwe', lrsthmatw, exst)
      !
      lrdgmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
      filint = TRIM(prefix)//'.dgmatwe'
      CALL diropn(iudgwe, 'dgmatwe', lrdgmatw, exst)
    ENDIF
    !
    ! Loop over the coarse q-points.
    !
    DO iq = 1, nqc_ahc
      !
      xxq = xqc_loc(:, iq)
      !
      epf17(:, :, :, :) = czero
      dwf17(:, :, :, :) = czero
      sthf17(:, :, :, :) = czero
      dgf17(:, :, :, :) = czero
      !
      ! ------------------------------------------------------
      ! dynamical matrix : Wannier -> Bloch
      ! ------------------------------------------------------
      !
      IF (.NOT. lifc) THEN
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq, uf, w2)
      ELSE
        CALL dynifc2blochf(nmodes, rws, nrws, xxq, uf, w2)
      ENDIF
      !
      ! ...then take into account the mass factors and square-root the frequencies...
      !
      DO nu = 1, nmodes
        !
        ! wf are the interpolated eigenfrequencies
        ! (omega on fine grid)
        !
        IF (w2(nu) > -eps8) THEN
          wf_coarse(nu, iq) =  DSQRT(ABS(w2(nu)))
        ELSE
          wf_coarse(nu, iq) = -DSQRT(ABS(w2(nu)))
        ENDIF
        !
        DO mu = 1, nmodes
          na = (mu - 1) / 3 + 1
          uf(mu, nu) = uf(mu, nu) / DSQRT(amass(ityp(na)))
        ENDDO
      ENDDO ! nu
      !
      ! --------------------------------------------------------------
      ! epmat : Wannier el and Wannier ph -> Wannier el and Bloch ph
      ! --------------------------------------------------------------
      !
      CALL ephwan2blochp(nmodes, xxq, irvec_g, ndegen_g, nrr_g, uf, epmatwef, nbndsub, nrr_k, dims, dims2)
      !
      ! Read sthmatwe and dgmatwe from file
      !
      sthmatwe = czero
      !
      IF (ionode) THEN
        DO imode = 1, nmodes
          CALL davcio(sthmatwe(:, :, :, :, imode), lrsthmatw, iusthwe, &
            imode + (iq - 1) * nmodes, -1)
        ENDDO
        CALL davcio(dgmatwe, lrdgmatw, iudgwe, iq, -1)
      ENDIF
      !
      CALL mp_bcast(dgmatwe, ionode_id, world_comm)
      CALL mp_bcast(sthmatwe, ionode_id, world_comm)
      !
      ! This is a loop over k blocks in the pool (size of the local k-set)
      !
      DO ik = 1, nkf
        !
        ! xkf is assumed to be in crys coord
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        xkk = xkf(:, ikk)
        xkq2 = xkk + xxq
        !
        cufkk(:, :) = czero
        cufkq(:, :) = czero
        !
        CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1)
        CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkq2, 1, 0.0_DP, rdotk2, 1)
        !
        IF (use_ws) THEN
          DO iw = 1, dims
            DO iw2 = 1, dims
              DO ir = 1, nrr_k
                IF (ndegen_k(ir, iw2, iw) > 0) THEN
                  cfac(ir, iw2, iw)  = EXP(ci * rdotk(ir))  / ndegen_k(ir, iw2, iw)
                  cfacq(ir, iw2, iw) = EXP(ci * rdotk2(ir)) / ndegen_k(ir, iw2, iw)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ELSE
          cfac(:, 1, 1)  = EXP(ci * rdotk(:))  / ndegen_k(:, 1, 1)
          cfacq(:, 1, 1) = EXP(ci * rdotk2(:)) / ndegen_k(:, 1, 1)
        ENDIF
        !
        ! ------------------------------------------------------
        ! hamiltonian : Wannier -> Bloch
        ! ------------------------------------------------------
        !
        ! Kohn-Sham first, then get the rotation matricies for following interp.
        IF (eig_read) THEN
          CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf_ks(:, ikk), chw_ks, cfac, dims)
          CALL hamwan2bloch(nbndsub, nrr_k, cufkq, etf_ks(:, ikq), chw_ks, cfacq, dims)
        ENDIF
        !
        CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf(:, ikk), chw, cfac, dims)
        CALL hamwan2bloch(nbndsub, nrr_k, cufkq, etf(:, ikq), chw, cfacq, dims)
        !
        ! --------------------------------------------------------------
        ! epmat : Wannier el and Bloch ph -> Bloch el and Bloch ph
        ! --------------------------------------------------------------
        !
        CALL rrwan2bloch(nbndsub, nrr_k, cfac, dims, rrf)
        CALL ephwan2bloch(nbndsub, nrr_k, epmatwef, cufkk, cufkq, epmatf, nmodes, cfac, dims, xxq, rrf)
        !
        ! Store epmatf in memory
        !
        DO jbnd = ibndmin, ibndmax
          DO ibnd = ibndmin, ibndmax
            epf17(ibnd - ibndmin + 1, jbnd - ibndmin + 1, :, ik) = epmatf(ibnd, jbnd, :)
          ENDDO
        ENDDO
        !
        ! Now do the eigenvector rotation: epmatf(j) = sum_i eptmp(i) * uf(i,j)
        tmp(:, :, :) = epf17(:, :, :, ik)
        CALL ZGEMM('n', 'n', nbndfst * nbndfst, nmodes, nmodes, cone, tmp(:,:, :), &
                   nbndfst * nbndfst, uf, nmodes, czero, epf17(:, :, :,ik), nbndfst * nbndfst)
        !
        ! -------------------------------------------------------------------
        ! dwmat : Wannier el, phonon Cartesian -> Bloch el, phonon eigenmode
        ! -------------------------------------------------------------------
        ! We calculate Debye-Waller only at q = Gamma.
        !
        dwmatf2(:, :, :) = czero
        !
        ! dwmatf (Debye-Waller in electron Bloch, phonon Cartesian basis)
        ! is independent of q.
        ! To reuse it, we need to store it for all k-points. (Not implemented)
        !
        ! Rotate dwmatwe from electron Wannier to Bloch basis
        ! We need to use cufkk and cufkk, not cufkk and cufkq
        ! nmodes must be set to 1 because we loop over idir and imode
        CALL ephwan2bloch(nbndsub, nrr_k, dwmatwe, cufkk, cufkk, dwmatf, &
            3*nmodes, cfac, dims, xxq, rrf, longrange = .FALSE.)
        !
        ! Subtract the active space term
        dwmatf = dwmatf - dwmatf_trunc(:, :, :, :, ik)
        !
        ! Rotate dwmatf (phonon Cartesian) to dwmatf2 (phonon eigenmode)
        CALL dwwan2blochp(nbndsub, dwmatf, uf, dwmatf2, nmodes)
        !
        ! Slim down to states inside the fsthick window.
        !
        DO jbnd = ibndmin, ibndmax
          DO ibnd = ibndmin, ibndmax
            dwf17(ibnd - ibndmin + 1, jbnd - ibndmin + 1, :, ik) = dwmatf2(ibnd, jbnd, :)
          ENDDO
        ENDDO
        !
        ! -------------------------------------------------------------------
        ! sthmat : el Wannier, phonon Cartesian -> el Bloch, phonon eigenmode
        ! -------------------------------------------------------------------
        !
        ! Rotate electron Wannier -> Bloch
        ! We need to use cufkk and cufkk, not cufkk and cufkq
        !
        sthmatf = czero
        CALL ephwan2bloch(nbndsub, nrr_k, sthmatwe, cufkk, cufkk, sthmatf, &
            nmodes**2, cfac, dims, xxq, rrf, longrange=.FALSE.)
        !
        ! Rotate sthmatf (phonon Cartesian) -> sthmatf2 (phonon eigenmode)
        !
        CALL sthwan2blochp(nbndsub, sthmatf, uf, sthmatf2, nmodes)
        !
        ! -------------------------------------------------------------------
        ! dgmat : el Wannier, phonon Cartesian -> el Bloch, phonon eigenmode
        ! -------------------------------------------------------------------
        !
        dgmatf = czero
        CALL ephwan2bloch(nbndsub, nrr_k, dgmatwe, cufkk, cufkq, dgmatf, &
            nmodes, cfac, dims, xxq, rrf, longrange=.FALSE.)
        !
        ! Rotate dgmatf (phonon Cartesian) -> dgf (phonon eigenmode)
        CALL dgwan2blochp(nbndsub, dgmatf, uf, dgmatf2, nmodes)
        !
        !
        ! Slim down to states inside the fsthick window.
        !
        DO jbnd = ibndmin, ibndmax
          DO ibnd = ibndmin, ibndmax
            sthf17(ibnd - ibndmin + 1, jbnd - ibndmin + 1, :, ik) = sthmatf2(ibnd, jbnd, :)
            dgf17(ibnd - ibndmin + 1, jbnd - ibndmin + 1, :, ik) = dgmatf2(ibnd, jbnd, :)
          ENDDO
        ENDDO
        !
      ENDDO ! ik
      !
      CALL selfen_elec_ahc_static(iq, nqc_ahc, wf_coarse(:, iq))
      !
    ENDDO ! iq
    !
    IF (ionode) THEN
      CLOSE(iusthwe, STATUS = 'keep')
      CLOSE(iudgwe, STATUS = 'keep')
    ENDIF
    !
    DEALLOCATE(wf_coarse, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating wf_coarse', 1)
    DEALLOCATE(xqc_loc, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating xqc_loc', 1)
    DEALLOCATE(dwmatwe, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating epmatwp', 1)
    DEALLOCATE(dwf17, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating dwf17', 1)
    DEALLOCATE(dwmatf, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating dwmatf', 1)
    DEALLOCATE(dwmatf2, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating dwmatf2', 1)
    DEALLOCATE(dgmatwe, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating dgmatwe', 1)
    DEALLOCATE(sthmatwe, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating sthmatwe', 1)
    DEALLOCATE(dgmatf, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating dgmatf', 1)
    DEALLOCATE(sthmatf, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating sthmatf', 1)
    DEALLOCATE(dgmatf2, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating dgmatf2', 1)
    DEALLOCATE(sthmatf2, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating sthmatf2', 1)
    DEALLOCATE(sthf17, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating sthf17', 1)
    DEALLOCATE(dgf17, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating dgf17', 1)
    !
    CALL stop_clock('ep-int-ahc')
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE ahc_run_static_wfpt
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE selfen_print
    !--------------------------------------------------------------------------
    !! Collect self-energy and print them to stdout and file
    !--------------------------------------------------------------------------
    USE, INTRINSIC :: ieee_arithmetic, ONLY: IEEE_VALUE, IEEE_QUIET_NAN
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : ryd2mev, ryd2ev, czero, zero, eps6, one, kelvin2eV
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE io_global,     ONLY : stdout, ionode
    USE epwcom,        ONLY : nbndsub, efermi_read, fermi_energy, nstemp, &
                              ahc_win_min, ahc_win_max
    USE elph2,         ONLY : etf, ibndmin, nkqf, nbndfst, xkf, nkqtotf, gtemp, &
                              sigma_ahc_act, sigma_ahc_ldw, sigma_ahc_hdw,      &
                              sigma_ahc_uf, nktotf, efnew
    USE poolgathering, ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ibnd
    !! Counter on bands at k
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
    REAL(KIND = DP), ALLOCATABLE :: tmp_sigmar_all(:, :, :)
    !! Real part of the total self-energy, temporarily used for printing
    REAL(KIND = DP), ALLOCATABLE :: tmp_sigmai_all(:, :, :)
    !! Imaginary part of the total self-energy, temporarily used for printing
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
    IF (ierr /= 0) CALL errore('selfen_print', 'Error allocating xkf_all', 1)
    ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('selfen_print', 'Error allocating etf_all', 1)
    xkf_all(:, :) = zero
    etf_all(:, :) = zero
    !
#if defined(__MPI)
    !
    ! Collect k-point distributed to pools.
    ! note that poolgather2 works with the doubled grid (k and k+q)
    !
    CALL poolgather2(3, nkqtotf, nkqf, xkf, xkf_all)
    CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
    CALL mp_sum(sigma_ahc_act, inter_pool_comm)
    CALL mp_sum(sigma_ahc_ldw, inter_pool_comm)
    CALL mp_sum(sigma_ahc_hdw, inter_pool_comm)
    CALL mp_sum(sigma_ahc_uf, inter_pool_comm)
    !
#else
    !
    xkf_all = xkf
    etf_all = etf
    !
#endif
    !
    ! Average over degenerate eigenstates:
    WRITE(stdout, '(5x,"Average over degenerate eigenstates is performed")')
    !
    CALL degenerate_average_cmplx(sigma_ahc_act, etf_all)
    CALL degenerate_average_real(sigma_ahc_ldw, etf_all)
    CALL degenerate_average_real(sigma_ahc_hdw, etf_all)
    CALL degenerate_average_real(sigma_ahc_uf, etf_all)
    !
    ! Output electron self-energy
    !
    IF (ionode) THEN
      !
      WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
      WRITE(stdout, '(5x, "Electron Self-Energy using the Allen-Heine-Cardona Theory")')
      WRITE(stdout, '(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
      WRITE(stdout, '(5x, a)') REPEAT('=', 67)
      !
      ! Write to file
      !
      ALLOCATE(tmp_sigmar_all(nbndfst, nktotf, nstemp))
      ALLOCATE(tmp_sigmai_all(nbndfst, nktotf, nstemp))
      tmp_sigmar_all = REAL(sigma_ahc_act, DP) + sigma_ahc_ldw + sigma_ahc_hdw + sigma_ahc_uf
      tmp_sigmai_all = IMAG(sigma_ahc_act)
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
            IF (ekk < ahc_win_min / ryd2ev .OR. ahc_win_max / ryd2ev < ekk) THEN
              sigma_ahc_uf(ibnd, ik, itemp) = IEEE_VALUE(sigma_ahc_uf(ibnd, ik, itemp), IEEE_QUIET_NAN)
              ! QE testcode cannot deal with NaNs. So, I set the values to zero instead of NaN.
              ! FIXME: Delete the next two lines once gitlab.com/QEF/q-e/-/issues/623 is fixed.
              sigma_ahc_uf(ibnd, ik, itemp) = zero
              tmp_sigmar_all(ibnd, ik, itemp) = zero
            ENDIF
            !
            !
          ENDDO
        ENDDO
      ENDDO
      !
      DO itemp = 1, nstemp
        !
        WRITE(stdout, '(5x, a, f8.3, a)') "Temperature: ", gtemp(itemp) * ryd2ev / kelvin2eV, "K"
        !
        DO ik = 1, nktotf
          !
          ikk = 2 * ik - 1
          !
          WRITE(stdout, '(/5x, "ik = ", i7," coord.: ", 3f12.7)') ik, xkf_all(:, ikk)
          WRITE(stdout, '(5x, a)') REPEAT('-', 67)
          !
          DO ibnd = 1, nbndfst
            !
            ! note that ekk does not depend on q
            ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0
            !
            WRITE(stdout, 103) ibndmin - 1 + ibnd, ryd2ev * ekk, &
              ryd2mev * tmp_sigmar_all(ibnd, ik, itemp), &
              ryd2mev * tmp_sigmai_all(ibnd, ik, itemp)
            !
          ENDDO
          WRITE(stdout, '(5x, a/)') REPEAT('-', 67)
          !
        ENDDO
        !
        WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
        WRITE(stdout, '(5x, a)') "ik ibnd E_nk [eV] Re[Sigma_active] &
                                 &Im[Sigma_active] Sigma_low_DW Sigma_high_DW Sigma_upperfan [meV]"
        WRITE(stdout, '(5x, a)') REPEAT('=', 67)
        !
        DO ibnd = 1, nbndfst
          DO ik = 1, nktotf
            !
            ikk = 2 * ik - 1
            !
            ekk = etf_all(ibndmin - 1 + ibnd, ikk)
            !
            WRITE(stdout, '(i9, i6, f12.6, 5f14.6)') &
              ik, ibndmin - 1 + ibnd, &
              ryd2ev * ekk, &
              ryd2mev * REAL(sigma_ahc_act(ibnd, ik, itemp)), &
              ryd2mev * AIMAG(sigma_ahc_act(ibnd, ik, itemp)), &
              ryd2mev * sigma_ahc_ldw(ibnd, ik, itemp), &
              ryd2mev * sigma_ahc_hdw(ibnd, ik, itemp), &
              ryd2mev * sigma_ahc_uf(ibnd, ik, itemp)
          ENDDO ! ik
          !
          WRITE(stdout, '(a)') '  '
          !
        ENDDO ! ibnd
        !
      ENDDO ! itemp
      !
      DEALLOCATE(tmp_sigmar_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_print', 'Error deallocating tmp_sigmar_all', 1)
      DEALLOCATE(tmp_sigmai_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_print', 'Error deallocating tmp_sigmai_all', 1)
      !
    ENDIF
    !
    DEALLOCATE(xkf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('selfen_print', 'Error deallocating xkf_all', 1)
    DEALLOCATE(etf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('selfen_print', 'Error deallocating etf_all', 1)
    !
103 FORMAT(5x, 'E( ', i3, ' )=', f9.4, ' eV   Re[Sigma]=', f15.6, ' meV Im[Sigma]=', &
           f15.6, ' meV')
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE selfen_print
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE degenerate_average_cmplx(array, etf_all)
    !--------------------------------------------------------------------------
    !! Average a complex-valued array over degenerate states
    !--------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : czero, eps6
    USE epwcom,        ONLY : nbndsub, nstemp
    USE elph2,         ONLY : ibndmin, nbndfst, nkqtotf, nktotf
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(inout) :: array(nbndfst, nktotf, nstemp)
    !! Quantity to be averaged
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndsub, nkqtotf)
    !! Collected eigenenergies
    !
    INTEGER :: n
    !! Integer for the degenerate average over eigenstates
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: itemp
    !! Counter on temperatures
    REAL(KIND = DP) :: ekk
    !! Eigen energy at k on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekk1
    !! Temporary variable to the eigenenergies at k for the degenerate average
    COMPLEX(KIND = DP) :: tmp
    !! Temporary variable to store quantity for the degenerate average
    COMPLEX(KIND = DP) :: array_tmp(nbndfst)
    !! Temporary array to store the Sigma
    !
    DO itemp = 1, nstemp
      DO ik = 1, nktotf
        !
        ikk = 2 * ik - 1
        !
        DO ibnd = 1, nbndfst
          ekk = etf_all(ibndmin - 1 + ibnd, ikk)
          n = 0
          tmp = czero
          DO jbnd = 1, nbndfst
            ekk1 = etf_all(ibndmin - 1 + jbnd, ikk)
            IF (ABS(ekk1 - ekk) < eps6) THEN
              n    = n + 1
              tmp = tmp + array(jbnd, ik, itemp)
            ENDIF
            !
          ENDDO ! jbnd
          array_tmp(ibnd) = tmp / REAL(n, DP)
          !
        ENDDO ! ibnd
        array(:, ik, itemp) = array_tmp(:)
        !
      ENDDO ! nktotf
    ENDDO ! itemp
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE degenerate_average_cmplx
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE degenerate_average_real(array, etf_all)
    !--------------------------------------------------------------------------
    !! Average a real-valued array over degenerate states
    !--------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : zero, eps6
    USE epwcom,        ONLY : nbndsub, nstemp
    USE elph2,         ONLY : ibndmin, nbndfst, nkqtotf, nktotf
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(inout) :: array(nbndfst, nktotf, nstemp)
    !! Quantity to be averaged
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndsub, nkqtotf)
    !! Collected eigenenergies
    !
    INTEGER :: n
    !! Integer for the degenerate average over eigenstates
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: itemp
    !! Counter on temperatures
    REAL(KIND = DP) :: ekk
    !! Eigen energy at k on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekk1
    !! Temporary variable to the eigenenergies at k for the degenerate average
    REAL(KIND = DP) :: tmp
    !! Temporary variable to store quantity for the degenerate average
    REAL(KIND = DP) :: array_tmp(nbndfst)
    !! Temporary array to store the Sigma
    !
    DO itemp = 1, nstemp
      DO ik = 1, nktotf
        !
        ikk = 2 * ik - 1
        !
        DO ibnd = 1, nbndfst
          ekk = etf_all(ibndmin - 1 + ibnd, ikk)
          n = 0
          tmp = zero
          DO jbnd = 1, nbndfst
            ekk1 = etf_all(ibndmin - 1 + jbnd, ikk)
            IF (ABS(ekk1 - ekk) < eps6) THEN
              n    = n + 1
              tmp = tmp + array(jbnd, ik, itemp)
            ENDIF
            !
          ENDDO ! jbnd
          array_tmp(ibnd) = tmp / REAL(n, DP)
          !
        ENDDO ! ibnd
        array(:, ik, itemp) = array_tmp(:)
        !
      ENDDO ! nktotf
    ENDDO ! itemp
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE degenerate_average_real
    !--------------------------------------------------------------------------
    !
  !-----------------------------------------------------------------------
  END MODULE selfen_ahc
  !-----------------------------------------------------------------------

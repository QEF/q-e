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
    USE io_var,        ONLY : linewidth_elself
    USE modes,         ONLY : nmodes
    USE epwcom,        ONLY : nstemp, nbndsub, shortrange, fsthick, ngaussw, degaussw, &
                              eps_acustic, efermi_read, fermi_energy, restart, restart_step
    USE pwcom,         ONLY : ef
    USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, eta, nbndfst, &
                              nkf, epf17, wf, wqf, xkf, nkqtotf, adapt_smearing, &
                              sigmar_all, sigmai_all, sigmai_mode, zi_all, efnew, &
                              nktotf, lower_bnd, gtemp
    USE control_flags, ONLY : iverbosity
    USE constants_epw, ONLY : kelvin2eV, ryd2mev, ryd2ev, one, two, zero, ci, eps6, eps8
    USE constants,     ONLY : pi
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE io_selfen,     ONLY : selfen_el_write
    USE poolgathering, ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 20) :: tp
    !! string for temperatures
    CHARACTER(LEN = 256) :: fileselfen
    !! file name for self energy
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
    INTEGER :: n
    !! Integer for the degenerate average over eigenstates
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
    !! Temporary variable to the eigenenergies at k for the degenerate average
    REAL(KIND = DP) :: ekq
    !! Eigen energy at k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: etmp1
    !! Temporary variable to store etmp1 = ekk - (ekq - wq)
    REAL(KIND = DP) :: etmp2
    !! Temporary variable to strore etmp2 = ekk - (ekq + wq)
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
    REAL(KIND = DP) :: w0g1
    !! Dirac delta at k for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: w0g2
    !! Dirac delta at k+q for the imaginary part of $\Sigma$
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
    REAL(KIND = DP) :: tmp1
    !! Temporary variable to store real part of Sigma for the degenerate average
    REAL(KIND = DP) :: tmp2
    !! Temporary variable to store imag part of Sigma for the degenerate average
    REAL(KIND = DP) :: tmp3
    !! Temporary variable to store Z for the degenerate average
    REAL(KIND = DP) :: sigmar_tmp(nbndfst)
    !! Temporary array to store the real-part of Sigma
    REAL(KIND = DP) :: sigmai_tmp(nbndfst)
    !! Temporary array to store the imag-part of Sigma
    REAL(KIND = DP) :: zi_tmp(nbndfst)
    !! Temporary array to store the Z
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
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! Collect k-point coordinate from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Collect eigenenergies from all pools in parallel case
    !
    ! SP: Define the inverse so that we can efficiently multiply instead of dividing
    inv_degaussw = one /degaussw
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
    DO itemp = 1, nstemp ! loop over temperatures
      inv_eptemp   = one / gtemp(itemp)
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
      IF (iqq == 1) THEN
        !
        WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
        WRITE(stdout, '(5x, "Electron (Imaginary) Self-Energy in the Migdal Approximation")')
        WRITE(stdout, '(5x, a/)') REPEAT('=', 67)
        !
        IF (fsthick < 1.d3) WRITE(stdout, '(/5x, a, f10.6, a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
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
      IF ((iqq == 1) .AND. .NOT. adapt_smearing) THEN
        WRITE(stdout, 100) degaussw * ryd2ev, ngaussw
        WRITE(stdout, '(a)') ' '
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
                eta_tmp     = eta2(ibnd, imode, ik)
                sq_eta_tmp  = eta_tmp**two
                inv_eta_tmp = inv_eta(ibnd, imode, ik)
                !
                DO jbnd = 1, nbndfst
                  !
                  ! the energy of the electron at k+q (relative to Ef)
                  ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
                  ! the Fermi occupation at k+q
                  wgkq = wgauss(-ekq * inv_eptemp, -99)
                  !
                  ! here we take into account the zero-point DSQRT(hbar/2M\omega)
                  ! with hbar = 1 and M already contained in the eigenmodes
                  ! g2 is Ry^2, wkf must already account for the spin factor
                  !
                  ! SP: Shortrange is disabled for efficiency reasons
                  !IF (shortrange .AND. ( ABS(xqf(1, iq))> eps8 .OR. ABS(xqf(2, iq))> eps8 &
                  !   .OR. ABS(xqf(3, iq))> eps8 )) THEN
                  !  ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary
                  !  !     number, in which case its square will be a negative number.
                  !  g2 = REAL((epf17(jbnd, ibnd, imode, ik)**two) * inv_wq(imode) * g2_tmp(imode))
                  !ELSE
                  !  g2 = (ABS(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq(imode) * g2_tmp(imode)
                  !ENDIF
                  !
                  g2 = (ABS(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq(imode) * g2_tmp(imode)
                  !
                  ! There is a sign error for wq in Eq. 9 of Comp. Phys. Comm. 181, 2140 (2010). - RM
                  ! The sign was corrected according to Eq. (7.282) page 489 from Mahan's book
                  ! (Many-Particle Physics, 3rd edition)
                  !
                  fact1 =       wgkq + wgq(imode)
                  fact2 = one - wgkq + wgq(imode)
                  etmp1 = ekk - (ekq - wq(imode))
                  etmp2 = ekk - (ekq + wq(imode))
                  sq_etmp1 = etmp1 * etmp1
                  sq_etmp2 = etmp2 * etmp2
                  !
                  weight = wqf(iq) * REAL(fact1 / (etmp1 - ci * eta_tmp) + fact2 / (etmp2 - ci * eta_tmp))
                  !
                  ! \Re\Sigma [Eq. 3 in Comput. Phys. Commun. 209, 116 (2016)]
                  sigmar_all(ibnd, ik + lower_bnd - 1, itemp) = sigmar_all(ibnd, ik + lower_bnd - 1, itemp) + g2 * weight
                  !
                  ! Logical implementation
                  ! weight = wqf(iq) * aimag(                                                  &
                  !         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
                  !           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) )
                  !
                  ! Delta implementation
                  w0g1 = w0gauss(etmp1 * inv_eta_tmp, 0) * inv_eta_tmp
                  w0g2 = w0gauss(etmp2 * inv_eta_tmp, 0) * inv_eta_tmp
                  !
                  weight = pi * wqf(iq) * (fact1 * w0g1 + fact2 * w0g2)
                  !
                  ! \Im\Sigma using delta approx. [Eq. 8 in Comput. Phys. Commun. 209, 116 (2016)]
                  sigmai_all(ibnd, ik + lower_bnd - 1, itemp) = sigmai_all(ibnd, ik + lower_bnd - 1, itemp) + g2 * weight
                  !
                  ! Mode-resolved
                  IF (iverbosity == 3) THEN
                    sigmai_mode(ibnd, imode, ik + lower_bnd - 1, itemp) = sigmai_mode(ibnd, imode, ik + lower_bnd - 1, itemp) + &
                            g2 * weight
                  ENDIF
                  !
                  ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
                  !
                  weight = wqf(iq) * &
                           (fact1 * (sq_etmp1 - sq_eta_tmp) / (sq_etmp1 + sq_eta_tmp)**two +  &
                            fact2 * (sq_etmp2 - sq_eta_tmp) / (sq_etmp2 + sq_eta_tmp)**two)
                  !
                  zi_all(ibnd, ik + lower_bnd - 1, itemp) = zi_all(ibnd, ik + lower_bnd - 1, itemp) + g2 * weight
                  !
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
            CALL mp_sum(sigmar_all, inter_pool_comm)
            CALL mp_sum(sigmai_all, inter_pool_comm)
            CALL mp_sum(zi_all, inter_pool_comm)
            CALL mp_sum(fermicount, inter_pool_comm)
            CALL mp_barrier(inter_pool_comm)
            CALL selfen_el_write(iqq, totq, nktotf, sigmar_all, sigmai_all, zi_all)
          ENDIF
        ENDIF
      ENDIF ! in case of restart, do not do the first one
    ENDDO ! itemp
    !
    ! The k points are distributed among pools: here we collect them
    !
    IF (iqq == totq) THEN
      !
      ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_elec_q', 'Error allocating xkf_all', 1)
      ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_elec_q', 'Error allocating etf_all', 1)
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
      IF (iverbosity == 3) CALL mp_sum(sigmai_mode, inter_pool_comm)
      CALL mp_sum(zi_all, inter_pool_comm)
      CALL mp_sum(fermicount, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
#else
      !
      xkf_all = xkf
      etf_all = etf
      !
#endif
      !
      DO itemp = 1, nstemp
        ! Average over degenerate eigenstates:
        WRITE(stdout, '(5x,"Average over degenerate eigenstates is performed")')
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
                n    = n + 1
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
        ! Output electron SE here after looping over all q-points (with their contributions
        ! summed in sigmar_all, etc.)
        !
        WRITE(stdout, '(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
        !
        IF (mpime == ionode_id) THEN
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
              WRITE(stdout, 102) ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2mev * sigmar_all(ibnd, ik, itemp), &
                                 ryd2mev * sigmai_all(ibnd,ik, itemp), zi_all(ibnd, ik, itemp), &
                                 one / zi_all(ibnd, ik, itemp) - one
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
        ENDIF
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
            WRITE(stdout, '(2i9, 5f12.4)') ik, ibndmin - 1 + ibnd, ryd2ev * ekk, ryd2mev * sigmar_all(ibnd, ik, itemp), &
                           ryd2mev * sigmai_all(ibnd, ik, itemp), zi_all(ibnd, ik, itemp), &
                           one / zi_all(ibnd, ik, itemp) - one
            !
          ENDDO
          !
          WRITE(stdout, '(a)') '  '
          !
        ENDDO
        !
      ENDDO ! itemp
      DEALLOCATE(xkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_elec_q', 'Error deallocating xkf_all', 1)
      DEALLOCATE(etf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('selfen_elec_q', 'Error deallocating etf_all', 1)
      !
    ENDIF
    !
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    102 FORMAT(5x, 'E( ', i3, ' )=', f9.4, ' eV   Re[Sigma]=', f15.6, ' meV Im[Sigma]=', &
               f15.6, ' meV     Z=', f15.6, ' lam=', f15.6)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE selfen_elec_q
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE selfen_phon_q(iqq, iq, totq)
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
    !! the selfenergy for one phonon at a time.  Much smaller footprint on
    !! the disk
    !!
    !! RM 24/02/2014
    !! redefined the size of coskkq, vkk, vkq within the fermi windwow
    !! cleaned up the subroutine
    !!
    !-----------------------------------------------------------------------
    USE kinds,      ONLY : DP
    USE io_global,  ONLY : stdout
    USE modes,      ONLY : nmodes
    USE epwcom,     ONLY : nbndsub, fsthick, efermi_read, fermi_energy,  &
                           nstemp, ngaussw, degaussw, shortrange,        &
                           nsmear, delta_smear, eps_acustic, specfun_ph, &
                           delta_approx, vme
    USE pwcom,      ONLY : nelec, ef
    USE klist_epw,  ONLY : isk_dummy
    USE elph2,      ONLY : epf17, ibndmax, ibndmin, etf, wkf, xqf, wqf, nkqf,   &
                           nkf, wf, nkqtotf, xqf, lambda_all, lambda_v_all,     &
                           dmef, vmef, gamma_all, gamma_v_all, efnew, nbndfst, &
                           gtemp, nktotf, adapt_smearing
    USE mp,         ONLY : mp_barrier, mp_sum
    USE mp_global,  ONLY : inter_pool_comm
    USE constants_epw, ONLY : kelvin2eV, ryd2mev, ryd2ev, one, two, zero, eps4, eps6, eps8
    USE constants,  ONLY : pi
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
    CHARACTER(LEN = 20) :: tp
    !! String for temperatures
    CHARACTER(LEN = 256) :: filephself
    !! File name for phonon self energy
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
    !! Counter on temperatuers
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
    REAL(KIND = DP) :: wgkk
    !! Fermi-Dirac occupation factor $f_{nk}(T)$
    REAL(KIND = DP) :: wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
    REAL(KIND = DP) :: weight
    !! Imaginary part of the phonhon self-energy factor
    !!$$ \pi N_q \Im(\frac{f_{nk}(T) - f_{mk+q(T)}}{\varepsilon_{nk}-\varepsilon_{mk+q}-\omega_{q\nu}+i\delta}) $$
    !! In practice the imaginary is performed with a delta Dirac
    REAL(KIND = DP) :: w0g1
    !! Dirac delta at k for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: w0g2
    !! Dirac delta at k+q for the imaginary part of $\Sigma$
    REAL(KIND = DP) :: degaussw0
    !! degaussw0 = (ismear-1) * delta_smear + degaussw
    REAL(KIND = DP) :: inv_degaussw0
    !! Inverse degaussw0 for efficiency reasons
    REAL(KIND = DP) :: inv_eptemp
    !! Inverse of temperature define for efficiency reasons
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
    REAL(KIND = DP) :: wq(nmodes)
    !! Phonon frequency on the fine grid
    REAL(KIND = DP) :: inv_wq(nmodes)
    !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
    REAL(KIND = DP) :: g2_tmp(nmodes)
    !! If the phonon frequency is too small discart g
    REAL(KIND = DP) :: gamma(nmodes)
    !! Gamma is the imaginary part of the phonon self-energy
    REAL(KIND = DP) :: gamma_v(nmodes)
    !! Gamma is the imaginary part of the phonon self-energy multiplied by (1-coskkq)
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
    IF (adapt_smearing) CALL errore('selfen_phon_q', 'adapt_smearing cannot be used with phonon self-energy', 1)
    !
    !
    DO itemp = 1, nstemp
      IF (iq == 1) THEN
        WRITE(stdout, '(/5x, a)') REPEAT('=',67)
        WRITE(stdout, '(5x, "Phonon (Imaginary) Self-Energy in the Migdal Approximation")')
        WRITE(stdout, '(5x, a/)') REPEAT('=',67)
        !
        IF (fsthick < 1.d3 ) WRITE(stdout, '(/5x, a, f10.6, a)' ) &
             'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
        WRITE(stdout, '(/5x, a, f10.6, a)' ) 'Golden Rule strictly enforced with T = ', gtemp(itemp) * ryd2ev, ' eV'
        !
      ENDIF
      !
      !
      ! Now pre-treat phonon modes for efficiency
      ! Treat phonon frequency and Bose occupation
      wq(:) = zero
      DO imode = 1, nmodes
        wq(imode) = wf(imode, iq)
        IF (wq(imode) > eps_acustic) THEN
          g2_tmp(imode) = one
          inv_wq(imode) = one / (two * wq(imode))
        ELSE
          g2_tmp(imode) = zero
          inv_wq(imode) = zero
        ENDIF
      ENDDO
      !
      DO ismear = 1, nsmear
        !
        degaussw0 = (ismear - 1) * delta_smear + degaussw
        !
        ! SP: Multiplication is faster than division ==> Important if called a lot
        !     in inner loops
        inv_degaussw0 = one / degaussw0
        inv_eptemp   = one / gtemp(itemp)
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
          IF (vme == 'wannier') THEN
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
          ELSE
            DO ibnd = 1, nbndfst
              DO jbnd = 1, nbndfst
                !
                ! v_(k,i) = 1/m <ki|p|ki> = dmef (:, i,i,k)
                ! 1/m  = 2 in Rydberg atomic units
                !
                vkk(:, ibnd) = REAL(dmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
                vkq(:, jbnd) = REAL(dmef(:, ibndmin - 1 + jbnd, ibndmin - 1 + jbnd, ikq))
                IF (ABS(vkk(1, ibnd)**two + vkk(2, ibnd)**two + vkk(3, ibnd)**two) > eps4) &
                  coskkq(ibnd, jbnd) = DDOT(3, vkk(:, ibnd), 1, vkq(:, jbnd), 1)  / &
                                       DDOT(3, vkk(:, ibnd), 1, vkk(:, ibnd), 1)
              ENDDO
            ENDDO
          ENDIF
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
                ELSE
                  wgkk = wgauss(-ekk * inv_eptemp, -99)
                ENDIF
                !
                DO jbnd = 1, nbndfst
                  !
                  !  the fermi occupation for k+q
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
                  IF (delta_approx) THEN
                    !
                    w0g2 = w0gauss(ekq * inv_degaussw0, 0) * inv_degaussw0
                    ! the expression below is positive-definite, but also an
                    ! approximation which neglects some fine features
                    weight = pi * wq(imode) * wkf(ikk) * w0g1 * w0g2
                    !
                  ELSE
                    !
                    wgkq = wgauss(-ekq * inv_eptemp, -99)
                    !
                    ! = k-point weight * [f(E_k) - f(E_k+q)] / [E_k+q - E_k - w_q + id]
                    ! This is the imaginary part of the phonon self-energy, sans
                    ! the matrix elements [Eq. 4 in Comput. Phys. Commun. 209, 116 (2016)]
                    !
                    !weight = wkf (ikk) * (wgkk - wgkq) * AIMAG(cone / (ekq - ekk - wq - ci * degaussw0))
                    !
                    ! SP: The expression below is the imag part of phonon self-energy,
                    ! sans matrix elements [Eq. 9 in Comput. Phys. Commun. 209, 116 (2016)]
                    !  = pi * k-point weight * [f(E_k) - f(E_k+q)] * delta[E_k+q - E_k - w_q]
                    !
                    weight = pi * wkf(ikk) * (wgkk - wgkq) * w0gauss((ekq - ekk - wq(imode)) * inv_degaussw0, 0) * inv_degaussw0
                    !
                  ENDIF
                  !
                  gamma(imode)   = gamma(imode)   + weight * g2
                  gamma_v(imode) = gamma_v(imode) + weight * g2 * (1.0d0 - coskkq(ibnd, jbnd))
                  !
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
        CALL mp_sum(gamma, inter_pool_comm)
        CALL mp_sum(gamma_v, inter_pool_comm)
        CALL mp_sum(fermicount, inter_pool_comm)
        CALL mp_barrier(inter_pool_comm)
        !
        ! An average over degenerate phonon-mode is performed.
        DO imode = 1, nmodes
          n = 0
          tmp1 = zero
          tmp2 = zero
          tmp3 = zero
          tmp4 = zero
          DO jmode = 1, nmodes
            IF (ABS(wq(imode) - wq(jmode)) < eps6) THEN
              n = n + 1
              IF (wq(jmode) > eps_acustic) THEN
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
                             ryd2mev * gamma_all(imode, iq, ismear, itemp), ryd2mev * wq(imode)
          WRITE(stdout, 104) imode, lambda_v_all(imode, iq, ismear, itemp), &
                             ryd2mev * gamma_v_all(imode, iq, ismear, itemp), ryd2mev * wq(imode)
          !
        ENDDO
        !
        WRITE(stdout, 103) lambda_tot
        WRITE(stdout, 105) lambda_tr_tot
        !
        IF (.NOT. specfun_ph) THEN
          WRITE(stdout, '(5x, a/)') REPEAT('-', 67)
          WRITE(stdout, '(/5x, a, i8, a, i8/)' ) 'Number of (k,k+q) pairs on the Fermi surface: ', fermicount, ' out of ', nktotf
        ENDIF
        !
      ENDDO !smears
      !
    ENDDO ! itemp
    100 FORMAT(5x, 'Gaussian Broadening: ', f10.6, ' eV, ngauss=', i4)
    101 FORMAT(5x, 'DOS =', f10.6, ' states/spin/eV/Unit Cell at Ef=', f10.6, ' eV')
    102 FORMAT(5x, 'lambda___( ', i3, ' )=', f15.6, '   gamma___=', f15.6, ' meV', '   omega=', f12.4, ' meV')
    103 FORMAT(5x, 'lambda___( tot )=', f15.6)
    104 FORMAT(5x, 'lambda_tr( ',i3,' )=', f15.6, '   gamma_tr=', f15.6, ' meV', '   omega=', f12.4, ' meV')
    105 FORMAT(5x, 'lambda_tr( tot )=', f15.6)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE selfen_phon_q
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
    USE epwcom,        ONLY : nbndsub, fsthick, ngaussw, efermi_read, &
                              fermi_energy, degaussw, nel, meff, epsiheg, &
                              restart, restart_step, nstemp
    USE pwcom,         ONLY : ef
    USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, dmef, adapt_smearing, &
                              nkf, wqf, xkf, nkqtotf, efnew, nbndfst, nktotf,  &
                              gtemp, sigmar_all, sigmai_all, zi_all, lower_bnd
    USE constants_epw, ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, ci, eps6, eps8
    USE constants,     ONLY : pi
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE cell_base,     ONLY : omega, alat, bg
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE io_selfen,     ONLY : selfen_el_write
    USE poolgathering, ONLY : poolgather2
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
    REAL(KIND = DP) :: inv_eptemp
    !! Inverse of temperature defined for efficiency reasons
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
      inv_eptemp   = one / gtemp(itemp)
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
                  wgkq = wgauss(-ekq * inv_eptemp, -99)
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
                      dipole = REAL(      dmef(1, ibndmin - 1 + jbnd, ibndmin - 1 + ibnd, ikk) / 2.d0 *  &
                                    CONJG(dmef(1, ibndmin - 1 + jbnd, ibndmin - 1 + ibnd, ikk) / 2.d0) / &
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
    USE constants_epw, ONLY : eps6, eps10
    USE constants,     ONLY : pi
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
    FUNCTION dos_ef_seq(ngauss, degauss, ef, et, wk, nks, nbnd)
    !-----------------------------------------------------------------------
    !
    USE kinds, ONLY : DP
    USE mp,    ONLY : mp_sum
    USE constants_epw, ONLY : zero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ngauss
    !! Number of smearing
    INTEGER, INTENT(in) :: nbnd
    !! Total number of bands considered
    INTEGER, INTENT(in) :: nks
    !!  Number of kpoints
    REAL(KIND = DP), INTENT(in) :: et(nbnd, nks)
    !! Eigenenergies
    REAL(KIND = DP), INTENT(in) :: wk(nks)
    !! K-point weights
    REAL(KIND = DP), INTENT(in) :: ef
    !! Fermi level
    REAL(KIND = DP), INTENT(in) :: degauss
    !! Smearing value
    !
    REAL(KIND = DP) :: dos_ef_seq
    !! Output of the function
    !
    ! Local variables
    INTEGER :: ik
    !! K-point value
    INTEGER :: ibnd
    !! Band number
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! Fermi-Dirac function
    !
    ! Compute DOS at E_F (states per Ry per unit cell)
    !
    dos_ef_seq = zero
    DO ik = 1, nks
      DO ibnd = 1, nbnd
        dos_ef_seq = dos_ef_seq + wk(ik) * w0gauss((et(ibnd, ik) - ef) / degauss, ngauss) / degauss
      ENDDO
    ENDDO
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION dos_ef_seq
    !-----------------------------------------------------------------------
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
    USE epwcom,    ONLY : nbndsub, fsthick, ngaussw, degaussw, &
                          nsmear, delta_smear, efermi_read, fermi_energy
    USE pwcom,     ONLY : ef
    USE elph2,     ONLY : ibndmin, etf, wkf, xqf, wqf, nkqf, nktotf, &
                          nkf, xqf, nbndfst, efnew
    USE constants_epw, ONLY : ryd2ev, zero, one, two
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

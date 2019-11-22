  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE transport
  !----------------------------------------------------------------------
  !! 
  !! This module contains routines linked with electronic transport  
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE scattering_rate_q(iqq, iq, totq, ef0, efcb, first_cycle) 
    !-----------------------------------------------------------------------
    !!
    !! This routine computes the scattering rate (inv_tau)
    !!
    ! 
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE phcom,         ONLY : nmodes
    USE epwcom,        ONLY : nbndsub, fsthick, eps_acustic, degaussw, restart,      & 
                              nstemp, scattering_serta, scattering_0rta, shortrange, &
                              restart_step, restart_filq, vme, assume_metal
    USE pwcom,         ONLY : ef
    USE elph2,         ONLY : ibndmin, etf, nkqf, nkf, dmef, vmef, wf, wqf, & 
                              epf17, nkqtotf, inv_tau_all, inv_tau_allcb,    &
                              xqf, zi_allvb, zi_allcb, nbndfst, nktotf, transp_temp, &
                              lower_bnd
    USE constants_epw, ONLY : zero, one, two, pi, ryd2mev, kelvin2eV, ryd2ev,        & 
                              eps6, eps8, eps4
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : world_comm
    USE io_transport,  ONLY : scattering_write, tau_write, merge_read
    USE poolgathering, ONLY : poolgather2
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: first_cycle
    !! Use to determine weather this is the first cycle after restart 
    INTEGER, INTENT(in) :: iqq
    !! Q-point index from the selected q
    INTEGER, INTENT(in) :: iq
    !! Q-point index
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points within the fstichk window
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(in) :: efcb(nstemp)
    !! Second Fermi level for the temperature itemp. Could be unused (0).
    !
    ! Local variables
    INTEGER :: n
    !! Integer for the degenerate average over eigenstates  
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ikq
    !! Even k+q index to read etf
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: jbnd
    !! Local band index
    INTEGER :: imode
    !! Local mode index
    INTEGER :: itemp
    !! Index over temperature range
    INTEGER :: nqtotf_new
    !! Number of q-point in the new dataset
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: tmp
    !! Temporary variable to store real part of Sigma for the degenerate average
    REAL(KIND = DP) :: tmp2
    !! Temporary variable for zi_all
    REAL(KIND = DP) :: ekk2
    !! Temporary variable to the eigenenergies for the degenerate average  
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND = DP) :: ekq
    !! Energy relative to Fermi level: $$\varepsilon_{m\mathbf{k+q}}-\varepsilon_F$$
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared (g2 is Ry^2) 
    REAL(KIND = DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP) :: w0g1
    !! $$ \delta[\varepsilon_{nk} - \varepsilon_{mk+q} + \omega_{q}] $$ 
    REAL(KIND = DP) :: w0g2 
    !! $$ \delta[\varepsilon_{nk} - \varepsilon_{mk+q} - \omega_{q}] $$
    REAL(KIND = DP) :: inv_wq 
    !! Inverse phonon frequency. Defined for efficiency reasons.
    REAL(KIND = DP) :: inv_etemp
    !! Invese temperature inv_etemp = 1/etemp. Defined for efficiency reasons.
    REAL(KIND = DP) :: g2_tmp 
    !! Used to set component to 0 if the phonon freq. is too low. This is defined
    !! for efficiency reasons as if statement should be avoided in inner-most loops.
    REAL(KIND = DP) :: inv_degaussw
    !! 1.0/degaussw. Defined for efficiency reasons. 
    REAL(KIND = DP) :: wq
    !! Phonon frequency $$\omega_{q\nu}$$ on the fine grid.  
    REAL(KIND = DP) :: wgq
    !! Bose-Einstein occupation function $$n_{q\nu}$$
    REAL(KIND = DP) :: weight
    !! Self-energy factor 
    REAL(KIND = DP) :: fmkq
    !! Fermi-Dirac occupation function $$f_{m\mathbf{k+q}}$$
    REAL(KIND = DP) :: trans_prob
    !! Transition probability function
    REAL(KIND = DP) :: vkk(3, nbndfst)
    !! Electronic velocity $$v_{n\mathbf{k}}$$
    REAL(KIND = DP) :: vkq(3, nbndfst)
    !! Electronic velocity $$v_{m\mathbf{k+q}}$$
    REAL(KIND = DP) :: vel_factor(nbndfst, nbndfst)
    !! Velocity factor  $$ 1 - \frac{(v_{nk} \cdot v_{mk+q})}{ |v_{nk}|^2} $$
    REAL(KIND = DP) :: inv_tau_tmp(nbndfst)
    !! Temporary array to store the scattering rates
    REAL(KIND = DP) :: zi_tmp(nbndfst)
    !! Temporary array to store the zi
    REAL(KIND = DP), ALLOCATABLE :: inv_tau_all_new (:, :, :)
    !! New scattering rates to be merged
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    REAL(KIND = DP), EXTERNAL :: DDOT
    !! Dot product function
    REAL(KIND = DP), EXTERNAL :: efermig
    !! Function that returns the Fermi energy
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac 
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function  
    ! 
    IF (assume_metal) THEN
      CALL errore("scattering_rate_q", "metals not implemented.", 1)
    ENDIF
    CALL start_clock('SCAT')
    ! 
    IF (iqq == 1) THEN
      !
      WRITE(stdout, '(/5x,a)') REPEAT('=',67)
      WRITE(stdout, '(5x,"Scattering rate")')
      WRITE(stdout, '(5x,a/)') REPEAT('=',67)
      !
      IF (fsthick < 1.d3) &
        WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
        WRITE(stdout, '(5x,a,f10.6,a)' ) 'This is computed with respect to the fine Fermi level ',ef * ryd2ev, ' eV'
        WRITE(stdout, '(5x,a,f10.6,a,f10.6,a)' ) 'Only states between ',(ef - fsthick) * ryd2ev, ' eV and ', &
                (ef + fsthick) * ryd2ev, ' eV will be included'
        WRITE(stdout, '(5x,a/)')
      !
    ENDIF
    ! 
    ! In the case of a restart do not add the first step
    IF (first_cycle) THEN
      first_cycle = .FALSE.
    ELSE
      ! loop over temperatures
      DO itemp = 1, nstemp
        !
        etemp = transp_temp(itemp)
        !
        ! SP: Define the inverse so that we can efficiently multiply instead of dividing
        !
        inv_etemp = 1.0 / etemp
        inv_degaussw = 1.0 / degaussw
        !
        DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          IF (scattering_0rta) THEN 
            !vel_factor = 1 - (vk dot vkq) / |vk|^2  appears in Grimvall 8.20
            vel_factor(:, :) = zero
            IF (vme) THEN 
              DO ibnd = 1, nbndfst
                !
                ! vkk(3,nbnd) - velocity for k
                vkk(:, ibnd) = REAL(vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
                !
                DO jbnd = 1, nbndfst
                  !
                  ! vkq(3,nbnd) - velocity for k + q
                  vkq(:, jbnd) = REAL(vmef(:, ibndmin - 1 + jbnd, ibndmin - 1 + jbnd, ikq))
                  !
                  IF (ABS(vkk(1, ibnd)**2 + vkk(2, ibnd)**2 + vkk(3, ibnd)**2 ) > eps4) &
                    vel_factor(ibnd, jbnd) = DDOT(3, vkk(:, ibnd), 1, vkq(:, jbnd), 1) / &
                                             DDOT(3, vkk(:, ibnd), 1, vkk(:, ibnd), 1)
                ENDDO
              ENDDO
            ELSE
              DO ibnd = 1, nbndfst
                !
                ! vkk(3,nbnd) - velocity for k
                vkk(:, ibnd) = 2.0 * REAL(dmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
                !
                DO jbnd = 1, nbndfst
                  ! 
                  ! vkq(3,nbnd) - velocity for k + q
                  vkq(:, jbnd) = 2.0 * REAL(dmef(:, ibndmin - 1 + jbnd, ibndmin - 1 + jbnd, ikq))
                  !
                  IF (ABS(vkk(1, ibnd)**2 + vkk(2, ibnd)**2 + vkk(3, ibnd)**2 ) > eps4) &
                    vel_factor(ibnd, jbnd) = DDOT(3, vkk(:, ibnd), 1, vkq(:, jbnd), 1) / &
                                             DDOT(3, vkk(:, ibnd), 1, vkk(:, ibnd), 1)
                ENDDO  
              ENDDO
            ENDIF
            vel_factor(:, :) = one - vel_factor(:, :)
          ENDIF
          !
          ! We are not consistent with ef from ephwann_shuffle but it should not 
          ! matter if fstick is large enough.
          IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
              (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
            !
            DO imode = 1, nmodes
              !
              ! the phonon frequency and bose occupation
              wq = wf(imode, iq)
              !
              ! SP : Avoid if statement in inner loops
              ! the coupling from Gamma acoustic phonons is negligible
              IF (wq > eps_acustic) THEN
                g2_tmp = 1.0
                wgq = wgauss(-wq * inv_etemp, -99)
                wgq = wgq / (one - two * wgq)
                ! SP : Define the inverse for efficiency
                inv_wq =  1.0 / (two * wq) 
              ELSE
                g2_tmp = 0.0
                wgq = 0.0
                inv_wq = 0.0
              ENDIF
              !
              DO ibnd = 1, nbndfst
                !
                !  energy at k (relative to Ef)
                ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0(itemp)
                !
                DO jbnd = 1, nbndfst
                  !
                  !  energy and fermi occupation at k+q
                  ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0(itemp)
                  fmkq = wgauss(-ekq * inv_etemp, -99)
                  !
                  ! here we take into account the zero-point DSQRT(hbar/2M\omega)
                  ! with hbar = 1 and M already contained in the eigenmodes
                  ! g2 is Ry^2, wkf must already account for the spin factor
                  !
                  ! In case of q=\Gamma, then the short-range = the normal g. We therefore 
                  ! need to treat it like the normal g with ABS(g).
                  IF (shortrange .AND. (ABS(xqf(1, iq)) > eps8 .OR. ABS(xqf(2, iq)) > eps8 &
                     .OR. ABS(xqf(3, iq)) > eps8)) THEN
                    ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                    !     number, in which case its square will be a negative number. 
                    g2 = REAL((epf17(jbnd, ibnd, imode, ik)**two) * inv_wq * g2_tmp, KIND = DP)
                  ELSE
                    g2 = (ABS(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq * g2_tmp
                  ENDIF
                  !
                  ! delta[E_k - E_k+q + w_q] and delta[E_k - E_k+q - w_q]
                  w0g1 = w0gauss((ekk - ekq + wq) * inv_degaussw, 0) * inv_degaussw
                  w0g2 = w0gauss((ekk - ekq - wq) * inv_degaussw, 0) * inv_degaussw
                  !
                  ! transition probability 
                  ! (2 pi/hbar) * (k+q-point weight) * g2 * 
                  ! { [f(E_k+q) + n(w_q)] * delta[E_k - E_k+q + w_q] + 
                  !   [1 - f(E_k+q) + n(w_q)] * delta[E_k - E_k+q - w_q] } 
                  !
                  ! DBSP Just to try
                  trans_prob = pi * wqf(iq) * g2 * ((fmkq + wgq) * w0g1 + (one - fmkq + wgq) * w0g2)
                  !
                  IF (scattering_serta) THEN 
                    ! energy relaxation time approximation 
                    inv_tau_all(itemp, ibnd, ik + lower_bnd - 1) = inv_tau_all(itemp, ibnd, ik + lower_bnd - 1) + two * trans_prob
                  ELSEIF (scattering_0rta) THEN 
                    ! momentum relaxation time approximation
                    inv_tau_all(itemp, ibnd, ik + lower_bnd - 1) = inv_tau_all(itemp, ibnd, ik + lower_bnd - 1) &
                                           + two * trans_prob * vel_factor(ibnd, jbnd)
                  ENDIF
                  !
                  ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
                  !
                  weight = wqf(iq) * &
                          ((      fmkq + wgq ) * ((ekk - (ekq - wq))**two - degaussw**two) /       &
                                                 ((ekk - (ekq - wq))**two + degaussw**two)**two +  &
                           (one - fmkq + wgq ) * ((ekk - (ekq + wq))**two - degaussw**two) /       &
                                                 ((ekk - (ekq + wq))**two + degaussw**two)**two)
                  !
                  zi_allvb(itemp, ibnd, ik + lower_bnd - 1) = zi_allvb(itemp, ibnd, ik + lower_bnd - 1) + g2 * weight
                  ! 
                ENDDO !jbnd
                !
              ENDDO !ibnd
              !
              ! In this case we are also computing the scattering rate for another Fermi level position
              ! This is used to compute both the electron and hole mobility at the same time.  
              IF (ABS(efcb(itemp)) > eps4) THEN
                ! 
                DO ibnd = 1, nbndfst
                  !
                  !  energy at k (relative to Ef)
                  ekk = etf(ibndmin - 1 + ibnd, ikk) - efcb(itemp)
                  !
                  DO jbnd = 1, nbndfst
                    !
                    !  energy and fermi occupation at k+q
                    ekq = etf(ibndmin - 1 + jbnd, ikq) - efcb(itemp)
                    fmkq = wgauss(-ekq * inv_etemp, -99)
                    !
                    ! here we take into account the zero-point DSQRT(hbar/2M\omega)
                    ! with hbar = 1 and M already contained in the eigenmodes
                    ! g2 is Ry^2, wkf must already account for the spin factor
                    !
                    ! In case of q=\Gamma, then the short-range = the normal g. We therefore 
                    ! need to treat it like the normal g with ABS(g).
                    IF (shortrange .AND. (ABS(xqf(1, iq)) > eps8 .OR. ABS(xqf(2, iq)) > eps8 &
                       .OR. ABS(xqf(3, iq)) > eps8)) THEN
                      ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                      !     number, in which case its square will be a negative number. 
                      g2 = REAL((epf17(jbnd, ibnd, imode, ik)**two) * inv_wq * g2_tmp, KIND = DP)
                    ELSE
                      g2 = (ABS(epf17(jbnd, ibnd, imode, ik))**two) * inv_wq * g2_tmp
                    ENDIF
                    !
                    ! delta[E_k - E_k+q + w_q] and delta[E_k - E_k+q - w_q]
                    w0g1 = w0gauss((ekk - ekq + wq) * inv_degaussw, 0) * inv_degaussw
                    w0g2 = w0gauss((ekk - ekq - wq) * inv_degaussw, 0) * inv_degaussw
                    !
                    ! transition probability 
                    ! (2 pi/hbar) * (k+q-point weight) * g2 * 
                    ! { [f(E_k+q) + n(w_q)] * delta[E_k - E_k+q + w_q] + 
                    !   [1 - f(E_k+q) + n(w_q)] * delta[E_k - E_k+q - w_q] } 
                    !
                    trans_prob = pi * wqf(iq) * g2 * &
                                 ((fmkq + wgq) * w0g1 + (one - fmkq + wgq) * w0g2)
                    !
                    IF (scattering_serta) THEN
                      ! energy relaxation time approximation 
                      inv_tau_allcb(itemp, ibnd, ik + lower_bnd - 1) = &
                      inv_tau_allcb(itemp, ibnd, ik + lower_bnd - 1) + two * trans_prob
                      !
                    ELSEIF (scattering_0rta) THEN
                      ! momentum relaxation time approximation
                      inv_tau_allcb(itemp, ibnd, ik + lower_bnd - 1) = inv_tau_allcb(itemp, ibnd, ik + lower_bnd - 1) &
                                             + two * trans_prob * vel_factor(ibnd, jbnd)
                    ENDIF
                    !
                    ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
                    !
                    weight = wqf(iq) * &
                            ((      fmkq + wgq) * ((ekk - (ekq - wq))**two - degaussw**two) /       &
                                                  ((ekk - (ekq - wq))**two + degaussw**two)**two +  &
                             (one - fmkq + wgq) * ((ekk - (ekq + wq))**two - degaussw**two) /       &
                                                  ((ekk - (ekq + wq))**two + degaussw**two)**two)
                    !
                    zi_allcb(itemp, ibnd, ik + lower_bnd - 1) = zi_allcb(itemp, ibnd, ik + lower_bnd - 1) + g2 * weight
                    ! 
                  ENDDO !jbnd
                ENDDO !ibnd
              ENDIF ! ABS(efcb) < eps4
            ENDDO !imode
          ENDIF ! endif  fsthick
        ENDDO ! end loop on k
      ENDDO ! itemp
      !
      ! Creation of a restart point
      IF (restart) THEN
        IF (MOD(iqq, restart_step) == 0) THEN
          WRITE(stdout, '(a)' ) '     Creation of a restart point'
          ! 
          ! The mp_sum will aggreage the results on each k-points. 
          CALL mp_sum(inv_tau_all, world_comm)
          CALL mp_sum(zi_allvb,    world_comm)
          !
          IF (ABS(efcb(1)) > eps4) THEN
            ! 
            CALL mp_sum(inv_tau_allcb, world_comm) 
            CALL mp_sum(zi_allcb,      world_comm) 
            ! 
          ENDIF
          ! 
          IF (ABS(efcb(1)) > eps4) THEN
            CALL tau_write(iqq, totq, nktotf, .TRUE.)
          ELSE
            CALL tau_write(iqq, totq, nktotf, .FALSE.)
          ENDIF
          ! 
          ! Now show intermediate mobility with that amount of q-points
          CALL transport_coeffs(ef0, efcb)
        ENDIF
      ENDIF
    ENDIF ! first_cycle
    !
    ! The k points are distributed among pools: here we collect them
    !
    IF (iqq == totq) THEN
      !
      ! The total number of k points
      !
      ALLOCATE(etf_all(nbndsub, nkqtotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('scattering_rate_q', 'Error allocating etf_all', 1)
      !
      CALL mp_sum(inv_tau_all, world_comm)
      IF (ABS(efcb(1)) > eps4) CALL mp_sum(inv_tau_allcb, world_comm)
      CALL mp_sum(zi_allvb, world_comm)
      IF (ABS(efcb(1)) > eps4) CALL mp_sum(zi_allcb, world_comm)
      !
#if defined(__MPI)
      !
      ! Collect contributions from all pools (sum over k-points)
      ! this finishes the integral over the BZ (k)
      !
      CALL poolgather2(nbndsub, nkqtotf, nkqf, etf, etf_all)
#else
      !
      etf_all = etf
#endif
      !
      DO itemp = 1, nstemp  
        ! 
        etemp = transp_temp(itemp)
        WRITE(stdout, '(a,f8.3,a)' ) '     Temperature ', etemp * ryd2ev / kelvin2eV, ' K'
        !
        ! In case we read another q-file, merge the scattering here
        IF (restart_filq /= '') THEN
          ! 
          ALLOCATE(inv_tau_all_new(nstemp, nbndfst, nktotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('scattering_rate_q', 'Error allocating inv_tau_all_new', 1)
          inv_tau_all_new(:, :, :) = zero
          ! 
          CALL merge_read(nktotf, nqtotf_new, inv_tau_all_new) 
          ! 
          inv_tau_all(:, :, :) = (inv_tau_all(:, :, :) * totq &
                              + inv_tau_all_new(:, :, :) * nqtotf_new) / (totq + nqtotf_new)
          DEALLOCATE(inv_tau_all_new, STAT = ierr)
          IF (ierr /= 0) CALL errore('scattering_rate_q', 'Error deallocating inv_tau_all_new', 1)
          !
          WRITE(stdout, '(a)' ) '     '
          WRITE(stdout, '(a,i10,a)' ) '     Merge scattering for a total of ', totq + nqtotf_new, ' q-points'
          ! 
          CALL tau_write(iqq + nqtotf_new, totq + nqtotf_new, nktotf, .FALSE.)
          WRITE(stdout, '(a)' ) '     Write to restart file the sum'
          WRITE(stdout, '(a)' ) '     '
          !
          ! 
        ENDIF
        ! Average over degenerate eigenstates:
        WRITE(stdout, '(5x,"Average over degenerate eigenstates is performed")')
        ! 
        DO ik = 1, nktotf
          ikk = 2 * ik - 1
          ikq = ikk + 1
          ! 
          DO ibnd = 1, nbndfst
            ekk = etf_all(ibndmin - 1 + ibnd, ikk)
            n = 0
            tmp = 0.0_DP
            tmp2 = 0.0_DP
            DO jbnd = 1, nbndfst
              ekk2 = etf_all(ibndmin - 1 + jbnd, ikk)
              IF (ABS(ekk2 - ekk) < eps6) THEN
                n = n + 1
                tmp = tmp + inv_tau_all(itemp, jbnd, ik)
                tmp2 = tmp2 + zi_allvb(itemp, jbnd, ik)
              ENDIF
              ! 
            ENDDO ! jbnd
            inv_tau_tmp(ibnd) = tmp / FLOAT(n)
            zi_tmp(ibnd) = tmp2 / FLOAT(n)
            !
          ENDDO ! ibnd
          inv_tau_all(itemp, :, ik) = inv_tau_tmp(:)
          zi_allvb(itemp, :, ik) = zi_tmp(:)
          ! 
        ENDDO ! nkqtotf
        !
        IF (ABS(efcb(itemp)) > eps4) THEN 
          ! Average over degenerate eigenstates:
          WRITE(stdout, '(5x,"Average over degenerate eigenstates in CB is performed")')
          ! 
          DO ik = 1, nktotf
            ikk = 2 * ik - 1 
            ikq = ikk + 1 
            ! 
            DO ibnd = 1, nbndfst
              ekk = etf_all(ibndmin - 1 + ibnd, ikk)
              n = 0 
              tmp = 0.0_DP
              tmp2 = 0.0_DP
              DO jbnd = 1, nbndfst
                ekk2 = etf_all(ibndmin - 1 + jbnd, ikk)
                IF (ABS(ekk2 - ekk) < eps6) THEN
                  n = n + 1 
                  tmp =  tmp + inv_tau_allcb(itemp, jbnd, ik)
                  tmp2 =  tmp2 + zi_allcb(itemp, jbnd, ik)
                ENDIF
                ! 
              ENDDO ! jbnd
              inv_tau_tmp(ibnd) = tmp / FLOAT(n)
              zi_tmp(ibnd) = tmp2 / FLOAT(n)
              !
            ENDDO ! ibnd
            inv_tau_allcb(itemp, :, ik) = inv_tau_tmp(:)
            zi_allcb(itemp, :, ik) = zi_tmp(:)
            ! 
          ENDDO ! nkqtotf
        ENDIF
        ! 
        ! Output scattering rates here after looping over all q-points
        ! (with their contributions summed in inv_tau_all, etc.)
        CALL scattering_write(itemp, etemp, ef0, etf_all)
        !
      ENDDO !nstemp 
      !
      DEALLOCATE(etf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('scattering_rate_q', 'Error deallocating etf_all', 1)
      ! 
      ! Creation of a restart point at the end
      IF (restart) THEN
        WRITE(stdout, '(a)' ) '     Creation of the final restart point'
        ! 
        IF (ABS(efcb(1)) > eps4) THEN
          CALL tau_write(iqq, totq, nktotf, .TRUE.)
        ELSE
          CALL tau_write(iqq, totq, nktotf, .FALSE.)
        ENDIF
        ! 
      ENDIF ! restart
      ! 
    ENDIF ! iqq 
    !
    CALL stop_clock ('SCAT')
    ! DBSP
    !write(stdout,*),'iqq ',iqq
    !print*,shape(inv_tau_all)
    !write(stdout,*),'inv_tau_all(1,5:8,21) ',SUM(inv_tau_all(3,5:8,1))
    !write(stdout,*),'inv_tau_all(1,5:8,:) ',SUM(inv_tau_all(3,5:8,:))
    !write(stdout,*),'SUM(inv_tau_all) ',SUM(inv_tau_all(3,:,:))
    !write(stdout,*),'first_cycle ',first_cycle
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE scattering_rate_q
    !-----------------------------------------------------------------------
    !       
    !-----------------------------------------------------------------------
    SUBROUTINE transport_coeffs(ef0, efcb)
    !-----------------------------------------------------------------------
    !!
    !!  This routine computes the transport coefficients
    !!  SP - June 2018 - Update for symmetries in velocities when using homogeneous grids. 
    !!       This is currently commented out since it is ONLY needed if we are
    !!       interested in the off-diagonal mobility_\alpha\beta terms. 
    !!       At the moment we just want mobility_\alpha\alpha so it makes no difference.
    !!
    !-----------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE io_global,        ONLY : stdout, meta_ionode_id
    USE cell_base,        ONLY : alat, at, omega
    USE io_files,         ONLY : prefix 
    USE io_var,           ONLY : iufilsigma 
    USE epwcom,           ONLY : nbndsub, fsthick, system_2d, nstemp,              &
                                 int_mob, ncarrier, scatread, iterative_bte, vme, assume_metal
    USE pwcom,            ONLY : ef 
    USE elph2,            ONLY : ibndmin, etf, nkf, wkf, dmef, vmef,      & 
                                 inv_tau_all, nkqtotf, inv_tau_allcb, transp_temp, &
                                 zi_allvb, zi_allcb, map_rebal, nbndfst, nktotf
    USE constants_epw,    ONLY : zero, one, bohr2ang, ryd2ev, electron_SI,         &
                                 kelvin2eV, hbar, Ang2m, hbarJ, ang2cm, czero
    USE mp,               ONLY : mp_sum, mp_bcast
    USE mp_global,        ONLY : world_comm
    USE mp_world,         ONLY : mpime
    USE symm_base,        ONLY : s, t_rev, time_reversal, set_sym_bl, nrot
    USE cell_base,        ONLY : bg
    USE mp,               ONLY : mp_bcast
    USE epwcom,           ONLY : mp_mesh_k, nkf1, nkf2, nkf3
    USE constants_epw,    ONLY : eps6, eps4
    USE io_transport,     ONLY : scattering_read
    USE division,         ONLY : fkbounds
    USE grid,             ONLY : kpoint_grid_epw
    USE kinds_epw,        ONLY : SIK2
    USE poolgathering,    ONLY : poolgatherc4, poolgather2
    USE noncollin_module, ONLY : noncolin
    !
    IMPLICIT NONE
    ! 
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(in) :: efcb(nstemp)
    !! Second Fermi level for the temperature itemp (could be 0)
    !
    ! Local variables
    INTEGER :: i
    !! Cartesian direction index 
    INTEGER :: j
    !! Cartesian direction index 
    INTEGER :: ij
    !! Cartesian coupled index for matrix. 
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: lower_bnd
    !! Lower bounds index after k or q paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k or q paral
    INTEGER :: ikbz
    !! k-point index that run on the full BZ
    INTEGER :: nb
    !! Number of points in the BZ corresponding to a point in IBZ    
    INTEGER :: ierr
    !! Error status
    INTEGER :: bztoibz_tmp(nkf1 * nkf2 * nkf3)
    !! Temporary mapping
    INTEGER :: bztoibz(nkf1 * nkf2 * nkf3)
    !! BZ to IBZ mapping
    INTEGER(SIK2) :: s_bztoibz(nkf1 * nkf2 * nkf3)
    !! symmetry 
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND = DP) :: dfnk
    !! Derivative Fermi distribution $$-df_{nk}/dE_{nk}$$
    REAL(KIND = DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP) :: tau 
    !! Relaxation time
    REAL(KIND = DP) :: conv_factor1
    !! Conversion factor for the conductivity 
    REAL(KIND = DP) :: inv_cell 
    !! Inverse of the volume in [Bohr^{-3}]
    REAL(KIND = DP) :: carrier_density
    !! Carrier density [nb of carrier per unit cell]
    REAL(KIND = DP) :: carrier_density_prt
    !! Carrier density [nb of carrier per unit cell] in cm^-3 unit
    REAL(KIND = DP) :: fnk
    !! Fermi-Dirac occupation function 
    REAL(KIND = DP) :: mobility
    !! Sum of the diagonalized mobilities [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_xx
    !! Mobility along the xx axis after diagonalization [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_yy
    !! Mobility along the yy axis after diagonalization [cm^2/Vs] 
    REAL(KIND = DP) :: mobility_zz
    !! Mobility along the zz axis after diagonalization [cm^2/Vs] 
    REAL(KIND = DP) :: vkk(3, nbndfst)
    !! Electron velocity vector for a band. 
    REAL(KIND = DP) :: sigma(9, nstemp)
    !! Conductivity matrix in vector form
    REAL(KIND = DP) :: sigmaZ(9, nstemp)
    !! Conductivity matrix in vector form with Znk
    REAL(KIND = DP) :: sigma_m(3, 3, nstemp)
    !! Conductivity matrix
    REAL(KIND = DP) :: sigma_up(3, 3)
    !! Conductivity matrix in upper-triangle
    REAL(KIND = DP) :: sigma_eig(3)
    !! Eigenvalues from the diagonalized conductivity matrix
    REAL(KIND = DP) :: sigma_vect(3, 3)
    !! Eigenvectors from the diagonalized conductivity matrix
    REAL(KIND = DP) :: znk
    !! Real Znk from \lambda_nk (called zi_allvb or zi_allcb)
    REAL(KIND = DP) :: tdf_sigma(9)
    !! Temporary file
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac 
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function 
    REAL(KIND = DP), EXTERNAL :: efermig
    !! Function that returns the Fermi energy
    CHARACTER(LEN = 256) :: filsigma
    !! File for the conductivity  
    REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    COMPLEX(KIND = DP), ALLOCATABLE :: dmef_all(:, :, :, :)
    !! dipole matrix elements on the fine mesh among all pools
    COMPLEX(KIND = DP), ALLOCATABLE :: vmef_all(:, :, :, :)
    !! velocity matrix elements on the fine mesh among all pools
    REAL(KIND = DP), ALLOCATABLE :: tdf_sigma_m(:, :, :, :)
    !! transport distribution function
    REAL(KIND = DP), ALLOCATABLE :: wkf_all(:)
    !! k-point weight on the full grid across all pools
    !  SP - Uncomment to use symmetries on velocities
    REAL(KIND = DP) :: v_rot(3)
    !! Rotated velocity by the symmetry operation
    REAL(KIND = DP) :: vk_cart(3)
    !! veloctiy in cartesian coordinate
    REAL(KIND = DP) :: sa(3, 3)
    !! Rotation matrix
    REAL(KIND = DP) :: sb(3, 3)
    !! Rotation matrix
    REAL(KIND = DP) :: sr(3, 3)
    !! Rotation matrix
    ! 
    IF (assume_metal) THEN
      CALL errore("transport_coeffs", "metals not implemented.", 1)
    ENDIF
    CALL start_clock('MOB')
    !
    inv_cell = 1.0d0 / omega
    ! for 2d system need to divide by area (vacuum in z-direction)
    IF (system_2d) inv_cell = inv_cell * at(3, 3) * alat
    ! 
    ! We can read the scattering rate from files. 
    IF (scatread) THEN
      conv_factor1 = electron_SI / (hbar * bohr2ang * Ang2m)
      Sigma_m(:, :, :) = zero
      !
      ! Compute the Fermi level 
      DO itemp = 1, nstemp
        ! 
        etemp = transp_temp(itemp)
        ! 
        ! Lets gather the velocities from all pools
#if defined(__MPI)
        IF (vme) THEN 
          ALLOCATE(vmef_all(3, nbndsub, nbndsub, nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('transport_coeffs', 'Error allocating vmef_all', 1)
          vmef_all(:, :, :, :) = czero
          CALL poolgatherc4(3, nbndsub, nbndsub, nkqtotf, 2 * nkf, vmef, vmef_all)
        ELSE
          ALLOCATE(dmef_all(3, nbndsub, nbndsub, nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('transport_coeffs', 'Error allocating dmef_all', 1)
          dmef_all(:, :, :, :) = czero
          CALL poolgatherc4(3, nbndsub, nbndsub, nkqtotf, 2 * nkf, dmef, dmef_all)
        ENDIF
        ALLOCATE(wkf_all(nkqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('transport_coeffs', 'Error allocating wkf_all', 1)
        wkf_all(:) = zero
        CALL poolgather2(1, nkqtotf, 2 * nkf, wkf, wkf_all)
#else
        IF (vme) THEN
          ALLOCATE(vmef_all(3, nbndsub, nbndsub, nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('transport_coeffs', 'Error allocating vmef_all', 1)
          vmef_all = vmef
        ELSE
          ALLOCATE(dmef_all(3, nbndsub, nbndsub, nkqtotf), STAT = ierr)
          IF (ierr /= 0) CALL errore('transport_coeffs', 'Error allocating dmef_all', 1)
          dmef_all = dmef
        ENDIF
#endif     
        ALLOCATE(tdf_sigma_m(3, 3, nbndfst, nkqtotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('transport_coeffs', 'Error allocating tdf_sigma_m', 1)
        tdf_sigma_m(:, :, :, :) = zero 
        ! 
        ! In this case, the sum over q has already been done. It should therefore be ok 
        ! to do the mobility in sequential. Each cpu does the same thing below
        ALLOCATE(etf_all(nbndsub, nktotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('transport_coeffs', 'Error allocating etf_all', 1)
        !
        CALL scattering_read(etemp, ef0(itemp), etf_all, inv_tau_all)
        ! 
        ! This is hole mobility. ----------------------------------------------------
        IF (int_mob .OR. (ncarrier < -1E5)) THEN
          IF (itemp == 1) THEN
            WRITE(stdout, '(/5x,a)') REPEAT('=',67)
            WRITE(stdout, '(5x,"Temp [K]  Fermi [eV]  Hole density [cm^-3]  Hole mobility [cm^2/Vs]")')
            WRITE(stdout, '(5x,a/)') REPEAT('=',67)
          ENDIF
          !      
          DO ik = 1, nktotf 
            ikk = 2 * ik - 1
            ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
            IF (MINVAL(ABS(etf_all(:, ik) - ef)) < fsthick) THEN
              DO ibnd = 1, nbndfst
                ! This selects only valence bands for hole conduction
                IF (etf_all(ibndmin - 1 + ibnd, ik) < ef0(itemp)) THEN
                  IF (vme) THEN
                    vkk(:, ibnd) = REAL(vmef_all(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk)) 
                  ELSE
                    vkk(:, ibnd) = 2.0 * REAL(dmef_all (:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
                  ENDIF
                  tau = one / inv_tau_all(itemp, ibnd, ik)
                  ekk = etf_all(ibndmin - 1 + ibnd, ik) - ef0(itemp)
                  ! 
                  DO j = 1, 3
                    DO i = 1, 3
                      tdf_sigma_m(i, j, ibnd, ik) = vkk(i, ibnd) * vkk(j, ibnd) * tau
                    ENDDO
                  ENDDO
                  !
                  ! derivative Fermi distribution
                  dfnk = w0gauss(ekk / etemp, -99) / etemp
                  !
                  ! electrical conductivity matrix
                  Sigma_m(:, :, itemp) = Sigma_m(:, :, itemp) +  wkf_all(ikk) * dfnk * tdf_sigma_m(:, :, ibnd, ik)
                  !
                ENDIF ! valence bands
              ENDDO ! ibnd
            ENDIF ! fstick
          ENDDO ! ik
          ! 
          carrier_density = 0.0
          ! 
          DO ik = 1, nktotf
            ikk = 2 * ik - 1
            DO ibnd = 1, nbndfst
              ! This selects only valence bands for hole conduction
              IF (etf_all(ibndmin - 1 + ibnd, ik) < ef0(itemp)) THEN
                !  energy at k (relative to Ef)
                ekk = etf_all(ibndmin - 1 + ibnd, ik) - ef0(itemp)
                fnk = wgauss(-ekk / etemp, -99)
                ! The wkf(ikk) already include a factor 2
                carrier_density = carrier_density + wkf_all(ikk) * (1.0d0 - fnk)
              ENDIF
            ENDDO
          ENDDO
          ! 
          ! Diagonalize the conductivity matrix
          CALL rdiagh(3, Sigma_m(:, :, itemp), 3, sigma_eig, sigma_vect)
          ! 
          mobility_xx  = (sigma_eig(1) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility_yy  = (sigma_eig(2) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility_zz  = (sigma_eig(3) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility = (mobility_xx + mobility_yy + mobility_zz) / 3
          ! carrier_density in cm^-1
          carrier_density_prt = carrier_density * inv_cell * (bohr2ang * ang2cm)**(-3)
          WRITE(stdout, '(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV, &
                  ef0(itemp) * ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
          WRITE(stdout, '(45x, 1E18.6, a)') mobility_yy, '  y-axis'
          WRITE(stdout, '(45x, 1E18.6, a)') mobility_zz, '  z-axis'
          WRITE(stdout, '(45x, 1E18.6, a)') mobility, '     avg'
          !
        ENDIF ! int_mob .OR. (ncarrier < -1E5)
        ! 
        ! This is electron mobility. ----------------------------------------------------
        IF (int_mob .OR. (ncarrier > 1E5)) THEN
          IF (itemp == 1) THEN
            WRITE(stdout, '(/5x,a)') REPEAT('=',67)
            WRITE(stdout, '(5x,"Temp [K]  Fermi [eV]  Electron density [cm^-3]  Electron mobility [cm^2/Vs]")')
            WRITE(stdout, '(5x,a/)') REPEAT('=',67)
          ENDIF
          !      
          tdf_sigma_m(:, :, :, :) = zero
          !
          DO ik = 1, nktotf
            ikk = 2 * ik - 1
            ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
            IF (MINVAL(ABS(etf_all(:, ik) - ef)) < fsthick) THEN
              DO ibnd = 1, nbndfst
                ! This selects only conduction bands for electron conduction
                IF (etf_all(ibndmin - 1 + ibnd, ik) > ef0(itemp)) THEN
                  IF (vme) THEN 
                    vkk(:, ibnd) = REAL(vmef_all(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
                  ELSE 
                    vkk(:, ibnd) = 2.0 * REAL(dmef_all(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
                  ENDIF
                  tau = one / inv_tau_all(itemp, ibnd, ik)
                  ekk = etf_all(ibndmin - 1 + ibnd, ik) -  ef0(itemp)
                  ! 
                  DO j = 1, 3
                    DO i = 1, 3
                      tdf_sigma_m(i, j, ibnd, ik) = vkk(i, ibnd) * vkk(j, ibnd) * tau
                    ENDDO
                  ENDDO
                  !
                  ! derivative Fermi distribution
                  dfnk = w0gauss(ekk / etemp, -99) / etemp
                  !
                  ! electrical conductivity matrix
                  Sigma_m(:, :, itemp) = Sigma_m(:, :, itemp) + wkf_all(ikk) * dfnk * tdf_sigma_m(:, :, ibnd, ik)
                  !
                ENDIF ! valence bands
              ENDDO ! ibnd
            ENDIF ! fstick
          ENDDO ! ik
          ! 
          carrier_density = 0.0
          ! 
          DO ik = 1, nktotf
            ikk = 2 * ik - 1
            DO ibnd = 1, nbndfst
              ! This selects only conduction bands for electron conduction
              IF (etf_all(ibndmin - 1 + ibnd, ik) > ef0(itemp)) THEN
                !  energy at k (relative to Ef)
                ekk = etf_all(ibndmin - 1 + ibnd, ik) - ef0(itemp)
                fnk = wgauss(-ekk / etemp, -99)
                ! The wkf(ikk) already include a factor 2
                carrier_density = carrier_density + wkf_all(ikk) * fnk 
              ENDIF
            ENDDO
          ENDDO
          ! 
          ! Diagonalize the conductivity matrix
          CALL rdiagh(3, Sigma_m(:, :, itemp), 3, sigma_eig, sigma_vect)
          ! 
          mobility_xx = (sigma_eig(1) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility_yy = (sigma_eig(2) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility_zz = (sigma_eig(3) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility = (mobility_xx + mobility_yy + mobility_zz) / 3
          ! carrier_density in cm^-1
          carrier_density_prt = carrier_density * inv_cell * (bohr2ang * ang2cm)**(-3)
          WRITE(stdout, '(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV, &
                  ef0(itemp) * ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
          WRITE(stdout, '(45x, 1E18.6, a)') mobility_yy, '  y-axis'
          WRITE(stdout, '(45x, 1E18.6, a)') mobility_zz, '  z-axis'
          WRITE(stdout, '(45x, 1E18.6, a)') mobility, '     avg'
          !
        ENDIF ! int_mob .OR. (ncarrier > 1E5)
        ! 
        IF (vme) THEN
          DEALLOCATE(vmef_all, STAT = ierr)
          IF (ierr /= 0) CALL errore('transport_coeffs', 'Error deallocating vmef_all', 1)
        ELSE
          DEALLOCATE(dmef_all, STAT = ierr)
          IF (ierr /= 0) CALL errore('transport_coeffs', 'Error deallocating dmef_all', 1)
        ENDIF
        DEALLOCATE(tdf_sigma_m, STAT = ierr)
        IF (ierr /= 0) CALL errore('transport_coeffs', 'Error deallocating tdf_sigma_m', 1)
        DEALLOCATE(etf_all, STAT = ierr)
        IF (ierr /= 0) CALL errore('transport_coeffs', 'Error deallocating etf_all', 1)
      ENDDO ! itemp
      !
    ELSE ! Case without reading the scattering rates from files.
      !
      !  SP - Uncomment to use symmetries on velocities
      IF (mp_mesh_k) THEN
        bztoibz(:) = 0
        s_bztoibz(:) = 0
        ! 
        CALL set_sym_bl()
        ! What we get from this call is bztoibz
        CALL kpoint_grid_epw(nrot, time_reversal, .FALSE., s, t_rev, nkf1, nkf2, nkf3, bztoibz, s_bztoibz)
        ! 
        IF (iterative_bte) THEN
          ! Now we have to remap the points because the IBZ k-points have been
          ! changed to to load balancing. 
          bztoibz_tmp(:) = 0
          DO ikbz = 1, nkf1 * nkf2 * nkf3
            bztoibz_tmp(ikbz) = map_rebal(bztoibz(ikbz))
          ENDDO
          bztoibz(:) = bztoibz_tmp(:) 
        ENDIF
        !
      ENDIF
      !
      ! This is hole mobility. In the case of intrinsic mobilities we can do both
      ! electron and hole mobility because the Fermi level is the same. This is not
      ! the case for doped mobilities.
      ! 
      ! find the bounds of k-dependent arrays in the parallel case in each pool
      CALL fkbounds(nktotf, lower_bnd, upper_bnd)
      ! 
      IF (int_mob .OR. (ncarrier < -1E5)) THEN
        ! 
        DO itemp = 1, nstemp
          !
          etemp = transp_temp(itemp)
          !
          IF (itemp == 1) THEN 
            !
            ! tdf_sigma_ij(ibnd,ik) = v_i(ik,ibnd) * v_j(ik,ibnd) * tau(ik,ibnd)
            ! i,j - cartesian components and ij combined (i,j) index
            ! 1 = (1,1) = xx, 2 = (1,2) = xy, 3 = (1,3) = xz
            ! 4 = (2,1) = yx, 5 = (2,2) = yy, 6 = (2,3) = yz
            ! 7 = (3,1) = zx, 8 = (3,2) = zy, 9 = (3,3) = zz
            ! this can be reduced to 6 if we take into account symmetry xy=yx, ...
            tdf_sigma(:) = zero
            sigma(:, :)  = zero
            sigmaZ(:, :) = zero
            !
          ENDIF
          !
          DO ik = 1, nkf
            !
            ikk = 2 * ik - 1
            !
            ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
            IF (MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) THEN
              !
              DO ibnd = 1, nbndfst
                !
                ! This selects only valence bands for hole conduction
                IF (etf(ibndmin - 1 + ibnd, ikk) < ef0(itemp)) THEN 
                  !
                  ! vkk(3,nbnd) - velocity for k
                  tdf_sigma(:) = zero
                  IF (vme) THEN
                    ! vmef is in units of Ryd * bohr
                    vkk(:, ibnd) = REAL(vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
                  ELSE 
                    ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
                    ! 1/m  = 2 in Rydberg atomic units
                    ! dmef is in units of 1/a.u. (where a.u. is bohr)
                    ! v_(k,i) is in units of Ryd * a.u.
                    vkk(:, ibnd) = 2.0 * REAL(dmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
                  ENDIF
                  ! Use symmetries on k-point (from Homogeneous grid only)
                  IF (mp_mesh_k) THEN
                    !
                    vk_cart(:) = vkk(:, ibnd)
                    ! 
                    ! Loop on full BZ 
                    nb = 0
                    DO ikbz = 1, nkf1 * nkf2 * nkf3
                      ! If the k-point from the full BZ is related by a symmetry operation 
                      ! to the current k-point, then take it.  
                      IF (bztoibz(ikbz) == ik + lower_bnd - 1) THEN
                        nb = nb + 1
                        ! Transform the symmetry matrix from Crystal to cartesian
                        sa(:, :) = DBLE(s(:, :, s_bztoibz(ikbz)))
                        sb       = MATMUL(bg, sa)
                        sr(:, :) = MATMUL(at, TRANSPOSE(sb))
                        CALL DGEMV('n', 3, 3, 1.d0, sr, 3, vk_cart(:), 1, 0.d0, v_rot(:), 1)
                        ij = 0
                        DO j = 1, 3
                          DO i = 1, 3
                            ij = ij + 1
                            ! The factor two in the weight at the end is to account for spin
                            IF (noncolin) THEN
                              tdf_sigma(ij) = tdf_sigma(ij) + (v_rot(i) * v_rot(j)) * 1.0 / (nkf1 * nkf2 * nkf3)
                            ELSE
                              tdf_sigma(ij) = tdf_sigma(ij) + (v_rot(i) * v_rot(j)) * 2.0 / (nkf1 * nkf2 * nkf3)
                            ENDIF
                          ENDDO
                        ENDDO
                      ENDIF
                    ENDDO ! ikbz 
                    IF (noncolin) THEN
                      IF (ABS(nb * 1.0 / (nkf1 * nkf2 * nkf3) - wkf(ikk)) > eps6) THEN
                        CALL errore('transport', ' The number of kpoint in the IBZ is not equal to the weight', 1)
                      ENDIF
                    ELSE
                      IF (ABS(nb * 2.0 / (nkf1 * nkf2 * nkf3) - wkf(ikk)) > eps6) THEN
                        CALL errore('transport', ' The number of kpoint in the IBZ is not equal to the weight', 1)
                      ENDIF
                    ENDIF 
                  ! withtout symmetries
                  ELSE
                    ! 
                    ij = 0
                    DO j = 1, 3
                      DO i = 1, 3
                        ij = ij + 1
                        tdf_sigma(ij) = vkk(i, ibnd) * vkk(j, ibnd) * wkf(ikk)
                      ENDDO
                    ENDDO
                  ENDIF ! mp_mesh_k
                  !
                  !  energy at k (relative to Ef)
                  ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0(itemp)
                  tau = one / inv_tau_all(itemp, ibnd, ik + lower_bnd - 1)
                  ! derivative Fermi distribution
                  ! (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T) 
                  dfnk = w0gauss(ekk / etemp, -99) / etemp
                  ! electrical conductivity
                  sigma(:, itemp) = Sigma(:, itemp) + dfnk * tdf_sigma(:) * tau
                  !
                  ! Now do the same but with Znk multiplied
                  ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
                  znk = one / (one + zi_allvb(itemp, ibnd, ik + lower_bnd - 1))
                  tau = one / (Znk * inv_tau_all(itemp, ibnd, ik + lower_bnd - 1))
                  sigmaz(:, itemp) = sigmaz(:, itemp) + dfnk * tdf_sigma(:) * tau
                  !
                ENDIF
              ENDDO ! ibnd
            ENDIF ! endif  fsthick
          ENDDO ! end loop on k
          !
          ! The k points are distributed among pools: here we collect them
          !
          CALL mp_sum(sigma(:, itemp), world_comm)
          CALL mp_sum(sigmaz(:, itemp), world_comm)
          !
        ENDDO ! nstemp
        !
        IF (mpime == meta_ionode_id) THEN
          filsigma = TRIM(prefix) // '_elcond_h'
          OPEN(iufilsigma, FILE = filsigma, FORM = 'formatted')
          WRITE(iufilsigma, '(a)') "# Electrical conductivity in 1/(Ohm * m)"
          WRITE(iufilsigma, '(a)') "#         Ef(eV)         Temp(K)        Sigma_xx        Sigma_xy        Sigma_xz" // & 
                                                    "       Sigma_yx         Sigma_yy        Sigma_yz " // &
                                                    "        Sigma_xz        Sigma_yz        Sigma_zz"
        ENDIF
        !
        conv_factor1 = electron_SI / (hbar * bohr2ang * Ang2m)
        !
        WRITE(stdout, '(/5x,a)') REPEAT('=',67)
        WRITE(stdout, '(5x,"Temp [K]  Fermi [eV]  Hole density [cm^-3]  Hole mobility [cm^2/Vs]")')
        WRITE(stdout, '(5x,a/)') REPEAT('=',67)
        ! 
        DO itemp = 1, nstemp
          etemp = transp_temp(itemp)
          ! sigma in units of 1/(a.u.) is converted to 1/(Ohm * m)
          IF (mpime ==  meta_ionode_id) THEN
            WRITE(iufilsigma, '(11E16.8)') ef0(itemp) * ryd2ev, etemp * ryd2ev / kelvin2eV, &
                                          conv_factor1 * sigma(:, itemp) * inv_cell
          ENDIF
          carrier_density = 0.0
          ! 
          DO ik = 1, nkf
            ikk = 2 * ik - 1
            DO ibnd = 1, nbndfst
              ! This selects only valence bands for hole conduction
              IF (etf(ibndmin - 1 + ibnd, ikk) < ef0(itemp)) THEN
                !  energy at k (relative to Ef)
                ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0(itemp)      
                fnk = wgauss(-ekk / etemp, -99)
                ! The wkf(ikk) already include a factor 2
                carrier_density = carrier_density + wkf(ikk) * (1.0d0 - fnk) 
              ENDIF
            ENDDO
          ENDDO 
          ! 
          CALL mp_sum(carrier_density, world_comm)
          !
          ! Diagonalize the conductivity matrix
          ! 1 = (1,1) = xx, 2 = (1,2) = xy, 3 = (1,3) = xz
          ! 4 = (2,1) = yx, 5 = (2,2) = yy, 6 = (2,3) = yz
          ! 7 = (3,1) = zx, 8 = (3,2) = zy, 9 = (3,3) = zz
          sigma_up(:, :) = zero
          sigma_up(1, 1) = sigma(1, itemp)
          sigma_up(1, 2) = sigma(2, itemp)
          sigma_up(1, 3) = sigma(3, itemp)
          sigma_up(2, 1) = sigma(4, itemp)
          sigma_up(2, 2) = sigma(5, itemp)
          sigma_up(2, 3) = sigma(6, itemp)
          sigma_up(3, 1) = sigma(7, itemp)
          sigma_up(3, 2) = sigma(8, itemp)
          sigma_up(3, 3) = sigma(9, itemp)
          ! 
          CALL rdiagh(3, sigma_up, 3, sigma_eig, sigma_vect)
          ! 
          !Sigma_diag = (Sigma(1,itemp)+Sigma(5,itemp)+Sigma(9,itemp))/3
          !Sigma_offdiag = (Sigma(2,itemp)+Sigma(3,itemp)+Sigma(4,itemp)+&
          !                 Sigma(6,itemp)+Sigma(7,itemp)+Sigma(8,itemp))/6
          mobility_xx = (sigma_eig(1) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility_yy = (sigma_eig(2) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility_zz = (sigma_eig(3) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility = (mobility_xx + mobility_yy + mobility_zz) / 3
          ! carrier_density in cm^-1
          carrier_density_prt = carrier_density * inv_cell * (bohr2ang * ang2cm)**(-3)
          WRITE(stdout, '(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV, &
                  ef0(itemp) * ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
          WRITE(stdout, '(45x, 1E18.6, a)') mobility_yy, '  y-axis'
          WRITE(stdout, '(45x, 1E18.6, a)') mobility_zz, '  z-axis'
          WRITE(stdout, '(45x, 1E18.6, a)') mobility, '     avg' 
          ! 
  !        ! Now do Znk ----------------------------------------------------------
  !        sigma_up(:, :) = zero
  !        sigma_up(1,1) = SigmaZ(1,itemp)
  !        sigma_up(1,2) = SigmaZ(2,itemp)
  !        sigma_up(1,3) = SigmaZ(3,itemp)
  !        sigma_up(2,1) = SigmaZ(4,itemp)
  !        sigma_up(2,2) = SigmaZ(5,itemp)
  !        sigma_up(2,3) = SigmaZ(6,itemp)
  !        sigma_up(3,1) = SigmaZ(7,itemp)
  !        sigma_up(3,2) = SigmaZ(8,itemp)
  !        sigma_up(3,3) = SigmaZ(9,itemp)
  !        CALL rdiagh(3,sigma_up,3,sigma_eig,sigma_vect)
  !        mobility_xx  = ( sigma_eig(1) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
  !        mobility_yy  = ( sigma_eig(2) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
  !        mobility_zz  = ( sigma_eig(3) * electron_SI * ( bohr2ang * ang2cm  )**2)  /( carrier_density * hbarJ)
  !        mobility = (mobility_xx+mobility_yy+mobility_zz)/3
  !        ! carrier_density in cm^-1
  ! DBSP - Z-factor
  !        WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev /kelvin2eV, &
  !                ef0(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis [Z]'
  !        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis [Z]'
  !        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis [Z]'
  !        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg [Z]'
          ! 
        ENDDO ! nstemp
        !
        IF (mpime == meta_ionode_id) CLOSE(iufilsigma)
        !
      ENDIF ! Hole mob
      ! 
      ! Now the electron conduction and mobilities
      ! 
      IF (int_mob .OR. (ncarrier > 1E5)) THEN
        DO itemp = 1, nstemp
          etemp = transp_temp(itemp)
          IF (itemp == 1) THEN
            tdf_sigma(:) = zero
            sigma(:, :)  = zero
            sigmaZ(:, :) = zero
          ENDIF
          DO ik = 1, nkf
            ikk = 2 * ik - 1
            IF (MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) THEN
              IF (ABS(efcb(itemp)) < eps4) THEN  
                DO ibnd = 1, nbndfst
                  ! This selects only cond bands for electron conduction
                  IF (etf(ibndmin - 1 + ibnd, ikk) > ef0(itemp)) THEN
                    tdf_sigma(:) = zero
                    IF (vme) THEN
                      ! vmef is in units of Ryd * bohr
                      vkk(:, ibnd) = REAL(vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
                    ELSE
                      ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
                      ! 1/m  = 2 in Rydberg atomic units
                      ! dmef is in units of 1/a.u. (where a.u. is bohr)
                      ! v_(k,i) is in units of Ryd * a.u.
                      vkk(:, ibnd) = 2.0 * REAL(dmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
                    ENDIF
                    IF (mp_mesh_k) THEN
                      !
                      vk_cart(:) = vkk(:, ibnd)
                      ! 
                      ! Loop on full BZ 
                      nb = 0
                      DO ikbz = 1, nkf1 * nkf2 * nkf3
                        ! If the k-point from the full BZ is related by a symmetry operation 
                        ! to the current k-point, then take it.  
                        IF (bztoibz(ikbz) == ik + lower_bnd - 1) THEN
                          nb = nb + 1
                          ! Transform the symmetry matrix from Crystal to cartesian
                          sa(:, :) = DBLE(s(:, :, s_bztoibz(ikbz)))
                          sb       = MATMUL(bg, sa)
                          sr(:, :) = MATMUL(at, TRANSPOSE(sb))
                          CALL DGEMV('n', 3, 3, 1.d0, sr, 3, vk_cart(:), 1, 0.d0, v_rot(:), 1)
                          ij = 0
                          DO j = 1, 3
                            DO i = 1, 3
                              ij = ij + 1
                              ! The factor two in the weight at the end is to account for spin
                              IF (noncolin) THEN
                                tdf_sigma(ij) = tdf_sigma(ij) + (v_rot(i) * v_rot(j)) * 1.0 / (nkf1 * nkf2 * nkf3)
                              ELSE
                                tdf_sigma(ij) = tdf_sigma(ij) + (v_rot(i) * v_rot(j)) * 2.0 / (nkf1 * nkf2 * nkf3)
                              ENDIF
                            ENDDO
                          ENDDO
                        ENDIF
                      ENDDO ! ikbz 
                      IF (noncolin) THEN
                        IF (ABS(nb * 1.0 / (nkf1 * nkf2 * nkf3) - wkf(ikk)) > eps6) THEN
                          CALL errore('transport', ' The number of kpoint in the IBZ is not equal to the weight', 1)
                        ENDIF
                      ELSE
                        IF (ABS(nb * 2.0 / (nkf1 * nkf2 * nkf3) - wkf(ikk)) > eps6) THEN
                          CALL errore('transport', ' The number of kpoint in the IBZ is not equal to the weight', 1)
                        ENDIF
                      ENDIF
                    ! withtout symmetries
                    ELSE
                      ! 
                      ij = 0
                      DO j = 1, 3
                        DO i = 1, 3
                          ij = ij + 1
                          tdf_sigma(ij) = vkk(i, ibnd) * vkk(j, ibnd) * wkf(ikk)
                        ENDDO
                      ENDDO
                    ENDIF ! mp_mesh_k
                    !
                    ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0(itemp)
                    tau = one / inv_tau_all(itemp, ibnd, ik + lower_bnd - 1)
                    dfnk = w0gauss(ekk / etemp, -99) / etemp
                    sigma(:, itemp) = Sigma(:, itemp) + dfnk * tdf_sigma(:) * tau
                    !
                    ! Now do the same but with Znk multiplied
                    ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
                    znk = one / (one + zi_allvb(itemp, ibnd, ik + lower_bnd - 1))
                    tau = one / (Znk * inv_tau_all(itemp, ibnd, ik + lower_bnd - 1))
                    sigmaz(:, itemp) = sigmaz(:, itemp) + dfnk * tdf_sigma(:) * tau
                  ENDIF
                ENDDO 
              ELSE ! In this case we have 2 Fermi levels
                DO ibnd = 1, nbndfst
                  ! This selects only cond bands for hole conduction
                  IF (etf(ibndmin - 1 + ibnd, ikk) > efcb(itemp)) THEN
                    ! 
                    !  SP - Uncomment to use symmetries on velocities
                    tdf_sigma(:) = zero
                    IF (vme) THEN
                      vkk(:, ibnd) = REAL(vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
                    ELSE
                      vkk(:, ibnd) = 2.0 * REAL(dmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
                    ENDIF
                    ! 
                    IF (mp_mesh_k) THEN
                      !
                      vk_cart(:) = vkk(:, ibnd)
                      ! 
                      ! Loop on full BZ 
                      nb = 0
                      DO ikbz = 1, nkf1 * nkf2 * nkf3
                        ! If the k-point from the full BZ is related by a symmetry operation 
                        ! to the current k-point, then take it.  
                        IF (bztoibz(ikbz) == ik + lower_bnd - 1) THEN
                          nb = nb + 1
                          ! Transform the symmetry matrix from Crystal to cartesian
                          sa(:, :) = DBLE(s(:, :, s_bztoibz(ikbz)))
                          sb       = MATMUL(bg, sa)
                          sr(:, :) = MATMUL(at, TRANSPOSE(sb))
                          CALL DGEMV('n', 3, 3, 1.d0, sr, 3, vk_cart(:), 1, 0.d0, v_rot(:), 1)
                          ij = 0
                          DO j = 1, 3
                            DO i = 1, 3
                              ij = ij + 1
                              ! The factor two in the weight at the end is to account for spin
                              IF (noncolin) THEN
                                tdf_sigma(ij) = tdf_sigma(ij) + (v_rot(i) * v_rot(j)) * 1.0 / (nkf1 * nkf2 * nkf3)
                              ELSE
                                tdf_sigma(ij) = tdf_sigma(ij) + (v_rot(i) * v_rot(j)) * 2.0 / (nkf1 * nkf2 * nkf3)
                              ENDIF
                            ENDDO
                          ENDDO
                        ENDIF
                      ENDDO ! ikbz 
                      IF (noncolin) THEN
                        IF (ABS(nb * 1.0 / (nkf1 * nkf2 * nkf3) - wkf(ikk)) > eps6) THEN
                          CALL errore('transport', ' The number of kpoint in the IBZ is not equal to the weight', 1)
                        ENDIF
                      ELSE
                        IF (ABS(nb * 2.0 / (nkf1 * nkf2 * nkf3) - wkf(ikk)) > eps6) THEN
                          CALL errore ('transport', ' The number of kpoint in the IBZ is not equal to the weight', 1)
                        ENDIF
                      ENDIF
                    ! withtout symmetries
                    ELSE
                      ! 
                      ij = 0
                      DO j = 1, 3
                        DO i = 1, 3
                          ij = ij + 1
                          tdf_sigma(ij) = vkk(i, ibnd) * vkk(j, ibnd) * wkf(ikk)
                        ENDDO
                      ENDDO
                      !
                    ENDIF ! mp_mesh_k
                    ekk = etf(ibndmin - 1 + ibnd, ikk) - efcb(itemp)
                    tau = one / inv_tau_allcb(itemp, ibnd, ik + lower_bnd - 1)
                    dfnk = w0gauss(ekk / etemp, -99) / etemp
                    sigma(:, itemp) = sigma(:, itemp) +  dfnk * tdf_sigma(:) * tau                      
                    ! 
                    ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
                    znk = one / (one + zi_allcb(itemp, ibnd, ik + lower_bnd - 1))
                    tau = one / (znk * inv_tau_allcb(itemp, ibnd, ik + lower_bnd - 1))
                    ij = 0
                    DO j = 1, 3
                      DO i = 1, 3
                        ij = ij + 1
                        tdf_sigma(ij) = vkk(i, ibnd) * vkk(j, ibnd) * tau
                      ENDDO
                    ENDDO
                    sigmaZ(:, itemp) = sigmaZ(:, itemp) + wkf(ikk) * dfnk * tdf_sigma(:)                  
                  ENDIF 
                ENDDO ! ibnd
              ENDIF ! etcb
            ENDIF ! endif  fsthick
          ENDDO ! end loop on k
          CALL mp_sum(sigma(:, itemp), world_comm)
          CALL mp_sum(sigmaZ(:, itemp), world_comm)
          ! 
        ENDDO ! nstemp
        IF (mpime == meta_ionode_id) THEN
          filsigma = TRIM(prefix) // '_elcond_e'
          OPEN(iufilsigma, FILE = filsigma, FORM = 'formatted')
          WRITE(iufilsigma, '(a)') "# Electrical conductivity in 1/(Ohm * m)"
          WRITE(iufilsigma, '(a)') "#         Ef(eV)         Temp(K)        Sigma_xx        Sigma_xy        Sigma_xz" // &
                                                    "       Sigma_yx         Sigma_yy        Sigma_yz " // &
                                                    "        Sigma_xz        Sigma_yz        Sigma_zz"
        ENDIF
        !
        conv_factor1 = electron_SI / ( hbar * bohr2ang * Ang2m )
        WRITE(stdout, '(/5x,a)') REPEAT('=',67)
        WRITE(stdout, '(5x,"Temp [K]  Fermi [eV]  Elec density [cm^-3]  Elec mobility [cm^2/Vs]")')
        WRITE(stdout, '(5x,a/)') REPEAT('=',67)
        DO itemp = 1, nstemp
          etemp = transp_temp(itemp)
          IF (mpime == meta_ionode_id) THEN
            ! sigma in units of 1/(a.u.) is converted to 1/(Ohm * m)
            IF (ABS(efcb(itemp)) < eps4) THEN 
              WRITE(iufilsigma, '(11E16.8)') ef0(itemp) * ryd2ev, etemp * ryd2ev / kelvin2eV, &
                                             conv_factor1 * sigma(:, itemp) * inv_cell
            ELSE
              WRITE(iufilsigma, '(11E16.8)') efcb(itemp) * ryd2ev, etemp * ryd2ev / kelvin2eV, &
                                             conv_factor1 * sigma(:, itemp) * inv_cell
            ENDIF
          ENDIF
          carrier_density = 0.0
          ! 
          DO ik = 1, nkf
            DO ibnd = 1, nbndfst
              ikk = 2 * ik - 1
              ! This selects only conduction bands for electron conduction
              IF (ABS(efcb(itemp)) < eps4) THEN 
                IF (etf(ibndmin - 1 + ibnd, ikk) > ef0(itemp)) THEN
                  ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0(itemp)
                  fnk = wgauss( -ekk / etemp, -99)
                  ! The wkf(ikk) already include a factor 2
                  carrier_density = carrier_density + wkf(ikk) * fnk
                ENDIF
              ELSE
                IF (etf(ibndmin - 1 + ibnd, ikk) > efcb(itemp)) THEN
                  ekk = etf(ibndmin - 1 + ibnd, ikk) - efcb(itemp)
                  fnk = wgauss(-ekk / etemp, -99)
                  ! The wkf(ikk) already include a factor 2
                  carrier_density = carrier_density + wkf(ikk) * fnk
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          CALL mp_sum(carrier_density, world_comm)
          ! Diagonalize the conductivity matrix
          ! 1 = (1,1) = xx, 2 = (1,2) = xy, 3 = (1,3) = xz
          ! 4 = (2,1) = yx, 5 = (2,2) = yy, 6 = (2,3) = yz
          ! 7 = (3,1) = zx, 8 = (3,2) = zy, 9 = (3,3) = zz
          sigma_up(:, :) = zero
          sigma_up(1, 1) = sigma(1, itemp)
          sigma_up(1, 2) = sigma(2, itemp)
          sigma_up(1, 3) = sigma(3, itemp)
          sigma_up(2, 1) = sigma(4, itemp)
          sigma_up(2, 2) = sigma(5, itemp)
          sigma_up(2, 3) = sigma(6, itemp)
          sigma_up(3, 1) = sigma(7, itemp)
          sigma_up(3, 2) = sigma(8, itemp)
          sigma_up(3, 3) = sigma(9, itemp)
          ! 
          CALL rdiagh(3, sigma_up, 3, sigma_eig, sigma_vect)
          ! 
          mobility_xx = (sigma_eig(1) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility_yy = (sigma_eig(2) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility_zz = (sigma_eig(3) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility = (mobility_xx + mobility_yy + mobility_zz) / 3
          !
          ! Carrier_density in cm^-1
          carrier_density_prt = carrier_density * inv_cell * (bohr2ang * ang2cm)**(-3)
          IF (ABS(efcb(itemp)) < eps4) THEN
            WRITE(stdout, '(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV,&
                                                     ef0(itemp) * ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
          ELSE
            WRITE(stdout, '(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV,&
                                                     efcb(itemp) * ryd2ev, carrier_density_prt, mobility_xx, '  x-axis'
          ENDIF
          WRITE(stdout, '(45x, 1E18.6, a)') mobility_yy, '  y-axis'
          WRITE(stdout, '(45x, 1E18.6, a)') mobility_zz, '  z-axis'
          WRITE(stdout, '(45x, 1E18.6, a)') mobility, '     avg'
          ! Issue warning if the material is anisotropic
         ! IF (Sigma_offdiag > 0.1*Sigma_diag) THEN
         !   WRITE(stdout,'(5x,a,1f10.5,a)') 'Warning: Sigma_offdiag = ',(Sigma_offdiag*100)/Sigma_diag, '% of Sigma_diag'
         ! ENDIF
          ! Now do the mobility with Znk factor ----------------------------------------------------------
          sigma_up(:, :) = zero
          sigma_up(1, 1) = sigmaZ(1, itemp)
          sigma_up(1, 2) = sigmaZ(2, itemp)
          sigma_up(1, 3) = sigmaZ(3, itemp)
          sigma_up(2, 1) = sigmaZ(4, itemp)
          sigma_up(2, 2) = sigmaZ(5, itemp)
          sigma_up(2, 3) = sigmaZ(6, itemp)
          sigma_up(3, 1) = sigmaZ(7, itemp)
          sigma_up(3, 2) = sigmaZ(8, itemp)
          sigma_up(3, 3) = sigmaZ(9, itemp)
          CALL rdiagh(3, sigma_up, 3, sigma_eig, sigma_vect)
          mobility_xx = (sigma_eig(1) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility_yy = (sigma_eig(2) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility_zz = (sigma_eig(3) * electron_SI * (bohr2ang * ang2cm)**2) / (carrier_density * hbarJ)
          mobility = (mobility_xx + mobility_yy + mobility_zz) / 3
          !
  ! DBSP - Z-factor
  !        IF (ABS(efcb(itemp)) < eps4) THEN
  !          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV,&
  !                                                   ef0(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis [Z]'
  !        ELSE
  !          WRITE(stdout,'(5x, 1f8.3, 1f12.4, 1E19.6, 1E19.6, a)') etemp * ryd2ev / kelvin2eV,&
  !                                                   efcb(itemp)*ryd2ev, carrier_density_prt, mobility_xx, '  x-axis [Z]'
  !        ENDIF
  !        WRITE(stdout,'(45x, 1E18.6, a)') mobility_yy, '  y-axis [Z]'
  !        WRITE(stdout,'(45x, 1E18.6, a)') mobility_zz, '  z-axis [Z]'
  !        WRITE(stdout,'(45x, 1E18.6, a)') mobility, '     avg [Z]'
          ! 
        ENDDO ! nstemp
        WRITE(stdout,'(5x)')
        WRITE(stdout,'(5x,"Note: Mobility are sorted by ascending values and might not correspond")')
        WRITE(stdout,'(5x,"                                         to the expected (x,y,z) axis.")')
        WRITE(stdout,'(5x)')
        !
        IF (mpime == meta_ionode_id) CLOSE(iufilsigma)
        ! 
      ENDIF ! Electron mobilities
    ENDIF ! scatread
    !
    CALL stop_clock('MOB')
    ! 
    WRITE(stdout,  * ) '    Total time so far'
    CALL print_clock('SCAT')
    CALL print_clock('MOB')
    WRITE(stdout, '(5x)')
    ! 
    RETURN
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE transport_coeffs
    !--------------------------------------------------------------------------
    ! 
  !--------------------------------------------------------------------------
  END MODULE transport
  !--------------------------------------------------------------------------

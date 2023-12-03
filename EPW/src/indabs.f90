  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE indabs
  !----------------------------------------------------------------------
  !!
  !! This module contains the routines related to phonon assisted absorption
  !! 12/03/2018 Kyle and E. Kioupakis: First implementation
  !! 08/04/2018 S. Ponce: Cleaning
  !! 17/09/2019 S. Ponce: Modularization and cleaning
  !! 07/2021 X. Zhang: Added free carrier absorption related subroutines
  !! 04/2022 X. Zhang: Merged with impurity matrix elements implementation
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE indabs_main(iq, totq, first_cycle, iq_restart)
    !-----------------------------------------------------------------------
    !!
    !! Main routine for phonon assisted absorption
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iuindabs
    USE modes,         ONLY : nmodes
    USE epwcom,        ONLY : nstemp, fsthick, degaussw, &
                              eps_acustic, efermi_read, fermi_energy,&
                              vme, omegamin, omegamax, omegastep, carrier, &
                              nomega, neta, restart, restart_step, &
                              ii_g, ii_n, ii_lscreen
    USE elph2,         ONLY : etf, ibndmin, nkf, epf17, wkf, nqtotf, wf, wqf, &
                              sigmar_all, efnew, gtemp, &
                              dmef, omegap, epsilon2_abs, epsilon2_abs_lorenz, vmef, &
                              nbndfst, nktotf, ef0_fca, partion, &
                              epsilon2_abs_all, epsilon2_abs_lorenz_all, epstf_therm, &
                              epsilon2_abs_imp, epsilon2_abs_lorenz_imp, eimpf17
    USE constants_epw, ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi, ci, eps6, czero, &
                              bohr2ang, ang2cm
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_world,      ONLY : mpime
    USE mp_global,     ONLY : inter_pool_comm
    USE cell_base,     ONLY : omega
    USE io_indabs,     ONLY : indabs_write
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! Q-point index
    INTEGER, INTENT(in) :: iq_restart
    !! Restart q points, determine if initialization
    LOGICAL, INTENT(inout) :: first_cycle
    !! Use to determine weather this is the first cycle after restart
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points from the selecq.fmt grid.
    !
    ! Local variables
    CHARACTER(LEN = 256) :: nameF
    !! Name of the file
    CHARACTER(LEN = 10) :: c
    !! Number of eta values, in string format
    CHARACTER(LEN = 256) :: format_string
    !! Format string
    CHARACTER(LEN = 20) :: tp
    !! Temperature, in string format
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
    INTEGER :: nksqtotf
    !! Total number of k+q points
    INTEGER :: iw
    !! Index for frequency
!    INTEGER :: nomega
    !! Number of points on the photon energy axis
    INTEGER :: mbnd
    !! Index for summation over intermediate bands
    INTEGER :: ipol
    !! Polarization direction
    INTEGER :: m
    !! Counter on denominator imaginary broadening values
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: ierr
    !! Error status
!    INTEGER, PARAMETER :: neta = 9
    !! Broadening parameter
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekmk
    !! Eigen energy on the fine grid relative to the Fermi level for the intermediate band
    REAL(KIND = DP) :: ekmq
    !! Eigen energy of k+q on the fine grid relative to the Fermi level for the intermediate band
    REAL(KIND = DP) :: wq(nmodes), nqv(nmodes)
    !! Phonon frequencies and phonon occupations on the fine grid
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: wgkk, wgkq
    !! Fermi-Dirac occupation factor $f_{nk+q}(T)$, $f_{nk}(T)$
    REAL(KIND = DP) :: weighta, weighte
    !!- delta function for absorption, emission
    REAL(KIND = DP) :: inv_eptemp0
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse of the smearing for efficiency reasons
    REAL(KIND = DP) :: pfac
    !! Occupation prefactors
    REAL(KIND = DP) :: pface
    !! Occupation prefactors
    REAL(KIND = DP) :: cfac
    !! Absorption prefactor
    REAL(KIND = DP) :: n_imp_au(nstemp)
    !! Density of charged impurity in a.u.
    REAL(KIND = DP) :: inveps
    !! Inverse of thermal eps
    REAL(KIND = DP) :: eta(9) = (/ 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5 /) / ryd2eV
    !! Imaginary broadening of matrix element denominators
    REAL(KIND = DP) :: etemp_fca
    !! Temperature for fermi level calculation
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! Function to compute the Density of States at the Fermi level
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    COMPLEX(KIND = DP) :: vkk(3, nbndfst, nbndfst)
    !!- Velocity matrix elements at k, k+q
    COMPLEX(KIND = DP) :: vkq(3, nbndfst, nbndfst)
    !!- Velocity matrix elements at k, k+q
    COMPLEX(KIND = DP) :: s1a(3), s1e(3), s2a(3), s2e(3)
    !! Transition probability function
    COMPLEX(KIND = DP) :: s1imp(3), s2imp(3)
    !! Charged impurity transition
    COMPLEX(KIND = DP) :: eimpf(nbndfst, nbndfst)
    !! Electron-charged-impurity matrix elements
    COMPLEX(KIND = DP) :: epf(nbndfst, nbndfst, nmodes)
    !! Generalized matrix elements for phonon-assisted absorption
    !
    ! SP: Define the inverse so that we can efficiently multiply instead of dividing
    !
    inv_degaussw = 1.0 / degaussw
    !
    ! Done outside
!    nomega = INT((omegamax - omegamin) / omegastep) + 1
    !
    ! 300 K
    ! Epsilon2 prefactor for velocity matrix elements, including factor of 2 for spin (weights for k-points are divided by 2 to be normalized to 1)
    ! C = 8*pi^2*e^2 = 8*pi^2*2 approx 157.9136704
    !
    cfac = 16.d0 * pi**2
    !
    IF (iq == iq_restart) THEN
      !
      ALLOCATE(omegap(nomega), STAT = ierr)
      IF (ierr /= 0) CALL errore('indabs', 'Error allocating omegap', 1)
      !
      IF (carrier .and. ii_g) THEN
        ALLOCATE(epsilon2_abs_imp(3, nomega, neta, nstemp), STAT = ierr)
        IF (ierr /= 0) CALL errore('indabs', 'Error allocating epsilon2_abs_imp', 1)
        ALLOCATE(epsilon2_abs_lorenz_imp(3, nomega, neta, nstemp), STAT = ierr)
        IF (ierr /= 0) CALL errore('indabs', 'Error allocating epsilon2_abs_lorenz_imp', 1)
      ENDIF
      !
      IF (iq_restart == 1) THEN
        epsilon2_abs_all = 0.d0
        epsilon2_abs_lorenz_all = 0.d0
      ENDIF
      epsilon2_abs = 0.d0
      epsilon2_abs_lorenz = 0.d0
      DO iw = 1, nomega
        omegap(iw) = omegamin + (iw - 1) * omegastep
      ENDDO
      IF (iq_restart == 1) THEN
        CALL dirabs()
      ENDIF
      WRITE(stdout, '(/5x,a/)') REPEAT('=',67)
      WRITE(stdout, '(5x,"Phonon-assisted absorption")')
      WRITE(stdout, '(5x,a/)') REPEAT('=',67)
      !
      IF (fsthick < 1.d3) WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
      WRITE(stdout, '(/5x,a)') 'The following temperatures are calculated:'
      DO itemp = 1, nstemp
        WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Temperature T = ', gtemp(itemp) * ryd2ev, ' eV'
      ENDDO
      !
      IF (carrier) THEN
        IF (ii_g) THEN
        WRITE(stdout, '(5x,"Impurity-assisted optics calculated for the given charge density.")')
          epsilon2_abs_imp = 0.d0
          epsilon2_abs_lorenz_imp = 0.d0
        ENDIF
        CALL conduc_fca(ef0_fca)
      ENDIF
    ENDIF
    !
    ! The total number of k points
    !
    nksqtotf = nktotf ! odd-even for k,k+q
    !
    IF (ii_g) THEN
      DO itemp = 1, nstemp
        n_imp_au(itemp) = partion(itemp) * ii_n * (bohr2ang * ang2cm)**(3.d0)
      ENDDO
    ENDIF
    !
    DO itemp = 1, nstemp
      IF (first_cycle .and. itemp == nstemp) THEN
        first_cycle = .false.
      ELSE
        !
        IF (carrier) THEN
          !
          ef0 = ef0_fca(itemp)
        ELSEIF (efermi_read) THEN
          !
          ef0 = fermi_energy
        ELSE
          !
          ef0 = efnew
        ENDIF
        !
        inv_eptemp0 = 1.0 / gtemp(itemp)
        !
        IF (carrier .and. ii_g .and. ii_lscreen) THEN
          inveps = 1.0 / epstf_therm(itemp)
        ENDIF
        !
        DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          IF (carrier .and. ii_g) THEN
            IF (ii_lscreen) THEN
              eimpf(:, :) = eimpf17(:, :, ik) * inveps
            ELSE
              eimpf(:, :) = eimpf17(:, :, ik)
            ENDIF
          ENDIF
          !
          DO imode = 1, nmodes
            !
            ! the phonon frequency at this q and nu
            wq(imode) = wf(imode, iq)
            !
            epf(:, :, imode) = epf17(:, :, imode,ik)
            IF (wq(imode) > eps_acustic) THEN
              nqv(imode) = wgauss(-wq(imode) / gtemp(itemp), -99)
              nqv(imode) = nqv(imode) / (one - two * nqv(imode))
            ENDIF
          ENDDO
          !
          ! RM - vme version should be checked
          IF (vme == 'wannier') THEN
            DO ibnd = 1, nbndfst
              DO jbnd = 1, nbndfst
                ! vmef is in units of Ryd * bohr
                vkk(:, ibnd, jbnd) = vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + jbnd, ikk)
                vkq(:, ibnd, jbnd) = vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + jbnd, ikq)
              ENDDO
            ENDDO
          ELSE
            DO ibnd = 1, nbndfst
              DO jbnd = 1, nbndfst
                ! Dme's already corrected for GW corrections in wan2bloch.f90
                vkk(:, ibnd, jbnd) = dmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + jbnd, ikk)
                vkq(:, ibnd, jbnd) = dmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + jbnd, ikq)
              ENDDO
            ENDDO
          ENDIF
          !
          DO ibnd = 1, nbndfst
            !  the energy of the electron at k (relative to Ef)
            ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0
            !
            IF (ABS(ekk) < fsthick) THEN
              !
              wgkk = wgauss(-ekk * inv_eptemp0, -99)
              !
              DO jbnd = 1, nbndfst
                !
                ! The fermi occupation for k+q
                ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0
                !
                IF (ABS(ekq) < fsthick .AND. ekq < ekk + wq(nmodes) + omegamax + 6.0 * degaussw) THEN
                  !
                  wgkq = wgauss(-ekq * inv_eptemp0, -99)
                  !
                  IF (ekq - ekk - wq(nmodes) - omegamax > 6.0 * degaussw) CYCLE
                  IF (ekq - ekk + wq(nmodes) - omegamin < - 6.0 * degaussw) CYCLE
                  !
                  DO imode = 1, nmodes
                    !
                    IF (wq(imode) > eps_acustic) THEN
                      !
                      DO m = 1, neta
                        s1a = czero
                        s1e = czero
                        s2a = czero
                        s2e = czero
                        !
                        DO mbnd = 1, nbndfst
                          !
                          ! The energy of the electron at k (relative to Ef)
                          ekmk = etf(ibndmin - 1 + mbnd, ikk) - ef0
                          ! The energy of the electron at k+q (relative to Ef)
                          ekmq = etf(ibndmin - 1 + mbnd, ikq) - ef0
                          !
                          s1a(:) = s1a(:) + epf(jbnd, mbnd,imode) * vkk(:, mbnd, ibnd) / &
                                   (ekmk  - ekq + wq(imode) + ci * eta(m))
                          s1e(:) = s1e(:) + epf(jbnd, mbnd,imode) * vkk(:, mbnd, ibnd) / &
                                   (ekmk  - ekq - wq(imode) + ci * eta(m))
                          s2a(:) =  s2a(:) + epf(mbnd, ibnd,imode) * vkq(:, jbnd, mbnd) / &
                                   (ekmq  - ekk - wq(imode)+ ci * eta(m))
                          s2e(:) =  s2e(:) + epf(mbnd, ibnd,imode) * vkq(:, jbnd, mbnd) / &
                                   (ekmq  - ekk + wq(imode)+ ci * eta(m))
                        ENDDO
                        !
                        pfac  =  nqv(imode)      * wgkk * (one - wgkq) - (nqv(imode) + one) * (one - wgkk) * wgkq
                        pface = (nqv(imode) + one) * wgkk * (one - wgkq) -  nqv(imode) * (one - wgkk) * wgkq
                        !
                        DO iw = 1, nomega
                          !
                          IF (ABS(ekq - ekk - wq(imode) - omegap(iw)) > 6.0 * degaussw .AND. &
                              ABS(ekq - ekk + wq(imode) - omegap(iw)) > 6.0 * degaussw) CYCLE
                          !
                          weighte = w0gauss((ekq - ekk - omegap(iw) + wq(imode)) / degaussw, 0) / degaussw
                          weighta = w0gauss((ekq - ekk - omegap(iw) - wq(imode)) / degaussw, 0) / degaussw
                          !
                          DO ipol = 1, 3
                            epsilon2_abs(ipol, iw, m, itemp) = epsilon2_abs(ipol, iw, m, itemp) + (wkf(ikk) / 2.0) * wqf(iq) * &
                                 cfac / omegap(iw)**2 * pfac  * weighta * ABS(s1a(ipol) + s2a(ipol))**2 / (2 * wq(imode) * omega)
                            epsilon2_abs(ipol, iw, m, itemp) = epsilon2_abs(ipol, iw, m, itemp) + (wkf(ikk) / 2.0) * wqf(iq) * &
                                 cfac / omegap(iw)**2 * pface * weighte * ABS(s1e(ipol) + s2e(ipol))**2 / (2 * wq(imode) * omega)
                            epsilon2_abs_lorenz(ipol, iw, m, itemp) = epsilon2_abs_lorenz(ipol, iw, m, itemp) + &
                                  (wkf(ikk) / 2.0) * wqf(iq) * &
                                 cfac / omegap(iw)**2 * pfac  * ABS(s1a(ipol) + s2a(ipol))**2 / (2 * wq(imode) * omega) * &
                                 (degaussw / (degaussw**2 + (ekq - ekk - omegap(iw) - wq(imode))**2)) / pi
                            epsilon2_abs_lorenz(ipol, iw, m, itemp) = epsilon2_abs_lorenz(ipol, iw, m, itemp)  + &
                                 (wkf(ikk) / 2.0) * wqf(iq) * &
                                 cfac / omegap(iw)**2 * pface * ABS(s1e(ipol) + s2e(ipol))**2 / (2 * wq(imode) * omega) * &
                                 (degaussw / (degaussw**2 + (ekq - ekk - omegap(iw) + wq(imode))**2 )) / pi
                          ENDDO ! ipol
                        ENDDO ! iw
                      ENDDO ! neta
                    ENDIF ! if wq > acoustic
                  ENDDO ! imode
                  !
                  IF (carrier .and. ii_g) THEN
                    !
                    DO m = 1, neta
                      !
                      s1imp = czero
                      s2imp = czero
                      !
                      DO mbnd = 1, nbndfst
                        !
                        ! The energy of the electron at k (relative to Ef)
                        ekmk = etf(ibndmin - 1 + mbnd, ikk) - ef0
                        ! The energy of the electron at k+q (relative to Ef)
                        ekmq = etf(ibndmin - 1 + mbnd, ikq) - ef0
                        !
                        s1imp(:) = s1imp(:) + eimpf(jbnd, mbnd) * vkk(:, mbnd, ibnd) / &
                                 (ekmk  - ekq + ci * eta(m))
                        s2imp(:) = s2imp(:) + eimpf(mbnd, ibnd) * vkq(:, jbnd, mbnd) / &
                                 (ekmq  - ekk + ci * eta(m))
                      ENDDO
                      !
                      pfac  =  wgkk - wgkq
                      !
                      DO iw = 1, nomega
                        !
                        IF (ABS(ekq - ekk  - omegap(iw)) > 6.0 * degaussw) CYCLE
                        !
                        weighta = w0gauss((ekq - ekk - omegap(iw)) / degaussw, 0) / degaussw
                        !
                        DO ipol = 1, 3
                          epsilon2_abs_imp(ipol, iw, m, itemp) = epsilon2_abs_imp(ipol, iw, m, itemp) + (wkf(ikk) / 2.0) * &
                               wqf(iq) * cfac / omegap(iw)**2 * pfac  * weighta * &
                               ABS(s1imp(ipol) + s2imp(ipol))**2 * n_imp_au(itemp)
                          epsilon2_abs_lorenz_imp(ipol, iw, m, itemp) = epsilon2_abs_lorenz_imp(ipol, iw, m, itemp) + &
                               (wkf(ikk) / 2.0) * wqf(iq) * &
                               cfac / omegap(iw)**2 * pfac  * ABS(s1imp(ipol) + s2imp(ipol))**2 * n_imp_au(itemp) * &
                               (degaussw / (degaussw**2 + (ekq - ekk - omegap(iw))**2)) / pi
                        ENDDO ! ipol
                      ENDDO ! iw
                    ENDDO ! neta
                  ENDIF ! carrier&imp
                  !
                ENDIF ! endif  ekq in fsthick
              ENDDO ! jbnd
            ENDIF  ! endif  ekk in fsthick
          ENDDO ! ibnd
        ENDDO ! ik
        IF (restart) THEN
          IF (MOD(iq, restart_step) == 0 .and. itemp == nstemp) THEN
            WRITE(stdout, '(5x, a, i10)' ) 'Creation of a restart point at ', iq
            CALL mp_sum(epsilon2_abs, inter_pool_comm)
            CALL mp_sum(epsilon2_abs_lorenz, inter_pool_comm)
            CALL mp_barrier(inter_pool_comm)
            epsilon2_abs_all = epsilon2_abs_all + epsilon2_abs
            epsilon2_abs_lorenz_all = epsilon2_abs_lorenz_all + epsilon2_abs_lorenz
            CALL indabs_write(iq, totq, epsilon2_abs_all, epsilon2_abs_lorenz_all)
            epsilon2_abs = 0.d0
            epsilon2_abs_lorenz = 0.d0
          ENDIF
        ENDIF
      ENDIF ! Skip first step in restart
    ENDDO ! itemp
    !
    ! The k points are distributed among pools: here we collect them
    !
    IF (iq == nqtotf) THEN
      !
#if defined(__MPI)
      !
      ! Note that poolgather2 works with the doubled grid (k and k+q)
      !
      CALL mp_barrier(inter_pool_comm)
      CALL mp_sum(epsilon2_abs, inter_pool_comm)
      CALL mp_sum(epsilon2_abs_lorenz, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      IF (carrier .and. ii_g) THEN
        CALL mp_barrier(inter_pool_comm)
        CALL mp_sum(epsilon2_abs_imp, inter_pool_comm)
        CALL mp_sum(epsilon2_abs_lorenz_imp, inter_pool_comm)
        CALL mp_barrier(inter_pool_comm)
      ENDIF
      !
#endif
      !
      !Add everything again to all
      epsilon2_abs_all = epsilon2_abs_all + epsilon2_abs
      epsilon2_abs_lorenz_all = epsilon2_abs_lorenz_all + epsilon2_abs_lorenz
      ! Output to stdout
      c = 'X'
      WRITE(c,"(i0)") neta
      format_string = "(5x,f15.6," // TRIM(c) // "E22.14)"
      !
      WRITE(stdout, '(5x,a)')
      WRITE(stdout, '(5x,a)') 'Phonon-assisted absorption versus energy'
      WRITE(stdout, '(5x,a,4f15.6)') 'Broadenings: ', eta(1:4)
      WRITE(stdout, '(5x,a,5f15.6)') '   ', eta(5:9)
      WRITE(stdout, '(5x,a)')
      WRITE(stdout, '(5x,a)') 'For the first Broadening and Temperature we have:'
      ! For test-farm checking purposes, only show m=1
      WRITE(stdout, '(5x,a)') 'Photon energy (eV), Imaginary dielectric function along x,y,z'
      DO iw = 1, nomega
        WRITE(stdout, '(5x,f15.6,3E22.14)') omegap(iw) * ryd2ev, (epsilon2_abs_all(ipol, iw, 1, 1), ipol = 1, 3)
      ENDDO
      WRITE(stdout, '(5x,a)')
      WRITE(stdout, '(5x,a)') 'Values with other broadenings for temperature X are reported in the files epsilon2_indabs_X.dat'
      WRITE(stdout, '(5x,a)')
      !
      ! Output to file
      DO itemp = 1,nstemp
        IF (mpime == ionode_id) THEN
          WRITE(c,"(i0)") neta + 1
          WRITE(tp,"(f8.1)") gtemp(itemp) * ryd2ev / kelvin2eV
          format_string = "("//TRIM(c) // "E22.14)"
          nameF = 'epsilon2_indabs_' // trim(adjustl(tp)) // 'K.dat'
          OPEN(UNIT = iuindabs, FILE = nameF)
          WRITE(iuindabs, '(a)') '# Phonon-assisted absorption versus energy'
          WRITE(iuindabs, '(a)') '# Photon energy (eV), Directionally-averaged imaginary dielectric function along x,y,z'
          DO iw = 1, nomega
            WRITE(iuindabs, format_string) omegap(iw) * ryd2ev, (SUM(epsilon2_abs_all(:, iw, m, itemp)) / 3.0d0, m = 1, neta)
          ENDDO
          CLOSE(iuindabs)
          !
          nameF = 'epsilon2_indabs_lorenz' // trim(adjustl(tp)) // 'K.dat'
          OPEN(UNIT = iuindabs, FILE = nameF)
          WRITE(iuindabs, '(a)') '# Phonon-assisted absorption versus energy'
          WRITE(iuindabs, '(a)') '# Photon energy (eV), Directionally-averaged imaginary dielectric function along x,y,z'
          DO iw = 1, nomega
            WRITE(iuindabs, format_string) omegap(iw) * ryd2ev, (SUM(epsilon2_abs_lorenz_all(:, iw, m, itemp)) / 3.0d0, m = 1, neta)
          ENDDO
          CLOSE(iuindabs)
          !
          IF (carrier .and. ii_g) THEN
            nameF = 'epsilon2_indabs_imp_' // trim(adjustl(tp)) // 'K.dat'
            OPEN(UNIT = iuindabs, FILE = nameF)
            WRITE(iuindabs, '(a)') '# Charged-impurity-assisted absorption versus energy'
            WRITE(iuindabs, '(a)') '# Photon energy (eV), Directionally-averaged imaginary dielectric function along x,y,z'
            DO iw = 1, nomega
              WRITE(iuindabs, format_string) omegap(iw) * ryd2ev, (SUM(epsilon2_abs_imp(:, iw, m, itemp)) / 3.0d0, m = 1, neta)
            ENDDO
            CLOSE(iuindabs)
            !
            nameF = 'epsilon2_indabs_lorenz_imp_' // trim(adjustl(tp)) // 'K.dat'
            OPEN(UNIT = iuindabs, FILE = nameF)
            WRITE(iuindabs, '(a)') '# Charged-impurity-assisted absorption versus energy'
            WRITE(iuindabs, '(a)') '# Photon energy (eV), Directionally-averaged imaginary dielectric function along x,y,z'
            DO iw = 1, nomega
              WRITE(iuindabs, format_string) omegap(iw) * ryd2ev, &
                      (SUM(epsilon2_abs_lorenz_imp(:, iw, m, itemp)) / 3.0d0, m = 1, neta)
            ENDDO
            CLOSE(iuindabs)
            !
          ENDIF! carrier&ii_g
        ENDIF
      ENDDO
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE indabs_main
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE dirabs()
    !-----------------------------------------------------------------------
    !! Xiao Zhang 03/2021
    !! This routine calculates the direct part of the imaginary dielectric
    !! function.
    !! Only independent particles scheme is implemented, as simple as
    !! Eq. (27) in Rohlfing, M and Louie, S. G. PRB 62.8 (2000)
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iudirabs
    USE modes,         ONLY : nmodes
    USE epwcom,        ONLY : nbndsub, nstemp, fsthick, degaussw, &
                              eps_acustic, efermi_read, fermi_energy,&
                              vme, omegamin, omegamax, omegastep, carrier, &
                              nomega, neta, lindabs
    USE elph2,         ONLY : etf, ibndmin, nkf, epf17, wkf, nqtotf, wf, wqf, &
                              sigmar_all, efnew, gtemp, &
                              dmef, omegap, epsilon2_abs_dir, epsilon2_abs_lorenz_dir, vmef, &
                              nbndfst, nktotf, ef0_fca
    USE constants_epw, ONLY : kelvin2eV, ryd2mev, one, ryd2ev, two, zero, pi, ci, eps6, czero
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE cell_base,     ONLY : omega
    USE mp_world,      ONLY : mpime
    !
    IMPLICIT NONE
    !
    ! Local variables
    CHARACTER(LEN = 256) :: nameF
    !! Name of the file
    CHARACTER(LEN = 10) :: c
    !! Number of eta values, in string format
    CHARACTER(LEN = 256) :: format_string
    !! Format string
    CHARACTER(LEN = 20) :: tp
    !! Temperature, in string format
    INTEGER :: ik
    !! Counter on k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    INTEGER :: iw
    !! Index for frequency
    INTEGER :: ipol
    !! Counter on polarization
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: ekki
    !! Energy of the fine grid relative to the Fermi level, band i
    REAL(KIND = DP) :: ekkj
    !! Energy of the fine grid relative to the Fermi level, band j
    REAL(KIND = DP) :: ef0
    !! Fermi level
    REAL(KIND = DP) :: wgki, wgkj
    !! Fermi-Dirac occupation factor f_{i,k} and f_{j,k}
    REAL(KIND = DP) :: weighta
    !! delta function for absorption
    REAL(KIND = DP) :: inv_degaussw
    !! Inverse of the smearing for efficiency reasons
    REAL(KIND = DP) :: inv_temp
    !! Inverse of the temperature
    REAL(KIND = DP) :: pfac
    !! Occupation factor
    REAL(KIND = DP) :: cfac
    !! Absorption prefactor
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! Function to compute the Density of States at the Fermi level
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    COMPLEX(KIND = DP) :: vkk(3, nbndsub, nbndsub)
    !!- Velocity matrix elements at k
    COMPLEX(KIND = DP) :: optmat(3)
    !! Transition probability function
    !
    !! Define inverse so that multiply is more efficient
    !
    inv_degaussw = 1.0 / degaussw
    ! Epsilon2 prefactor for velocity matrix elements, including factor of 2 for spin (weights for k-points are divided by 2 to be normalized to 1)
    ! C = 8*pi^2*e^2 = 8*pi^2*2 approx 157.9136704
    !
    cfac = 16.d0 * pi**2
    !
    WRITE(stdout, '(/5x,a/)') REPEAT('=',67)
    WRITE(stdout, '(5x,"Direct absorption with independent particle approximation")')
    WRITE(stdout, '(5x,a/)') REPEAT('=',67)
    !
    IF (fsthick < 1.d3) WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
    WRITE(stdout, '(/5x,a)') 'The following temperatures are calculated:'
    DO itemp = 1, nstemp
      WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Temperature T = ', gtemp(itemp) * ryd2ev, ' eV'
    ENDDO
!    IF (.NOT. lindabs) THEN
!      ALLOCATE(omegap(nomega), STAT = ierr)
!      IF (ierr /= 0) CALL errore('dirabs', 'Error allocating omegap', 1)
!      DO iw = 1, nomega
!        omegap(iw) = omegamin + (iw - 1) * omegastep
!      ENDDO
    ALLOCATE(epsilon2_abs_dir(3, nomega, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('dirabs', 'Error allocating epsilon2_abs_dir', 1)
    ALLOCATE(epsilon2_abs_lorenz_dir(3, nomega, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('dirabs', 'Error allocating epsilon2_abs_lorenz_dir', 1)
    epsilon2_abs_dir = 0.d0
    epsilon2_abs_lorenz_dir = 0.d0
    !
    DO itemp = 1, nstemp
      IF (carrier) THEN
        !
        ef0 = ef0_fca(itemp)
      ELSEIF (efermi_read) THEN
        !
        ef0 = fermi_energy
      ELSE
        !
        ef0 = efnew
      ENDIF
      !
      inv_temp = 1.0 / gtemp(itemp)
      !
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        IF (vme == 'wannier') THEN
          DO ibnd = 1, nbndsub
            DO jbnd = 1, nbndsub
              ! vmef is in units of Ryd * bohr
              vkk(:, ibnd, jbnd) = vmef(:, ibnd, jbnd, ikk)
            ENDDO
          ENDDO
         ELSE
          DO ibnd = 1, nbndsub
            DO jbnd = 1, nbndsub
              ! Dme's already corrected for GW corrections in wan2bloch.f90
              vkk(:, ibnd, jbnd) = dmef(:, ibnd, jbnd, ikk)
            ENDDO
          ENDDO
        ENDIF
        DO ibnd = 1, nbndsub
          !  the energy of the electron at k (relative to Ef)
          ekki = etf(ibnd, ikk) - ef0
          IF (ABS(ekki) < fsthick * 4.0) THEN
            !
            ! Occupation factor
            wgki = wgauss(-ekki * inv_temp, -99)
            !
            DO jbnd = 1, nbndsub
              !
              ekkj = etf(jbnd, ikk) - ef0
              !
              IF (ABS(ekkj) < fsthick * 4.0 .AND. ekkj < ekki + omegamax + 6.0 * degaussw) THEN
                !
                wgkj = wgauss(-ekkj * inv_temp, -99)
                !
                IF (ekkj - ekki - omegamax > 6.0 * degaussw) CYCLE
                IF (ekkj - ekki - omegamin < - 6.0 * degaussw) CYCLE
                !
                optmat(:) = vkk(:,jbnd,ibnd)
                pfac = wgki - wgkj
                DO iw = 1, nomega
                  IF (ABS(ekkj - ekki  - omegap(iw)) > 6.0 * degaussw) CYCLE
                  !
                  weighta = w0gauss((ekkj - ekki - omegap(iw)) / degaussw, 0) / degaussw
                  !
                  DO ipol = 1, 3
                    epsilon2_abs_dir(ipol, iw, itemp) = epsilon2_abs_dir(ipol, iw, itemp) + (wkf(ikk) / 2.0) * cfac / &
                        omegap(iw) ** 2 * pfac * weighta * ABS(optmat(ipol)) ** 2 / omega
                    epsilon2_abs_lorenz_dir(ipol, iw, itemp) = epsilon2_abs_lorenz_dir(ipol, iw, itemp) + (wkf(ikk) / 2.0) * &
                        cfac / omegap(iw) ** 2 * pfac * ABS(optmat(ipol)) ** 2  / omega * &
                        (degaussw / (degaussw**2 + (ekkj - ekki - omegap(iw))**2)) / pi
                  ENDDO ! ipol
                ENDDO ! iw
              ENDIF ! ekkj
            ENDDO !jbnd
          ENDIF ! ekki
        ENDDO !ibnd
      ENDDO !ik
    ENDDO ! itemp
    !
    ! Sum the k-points across pool
    !
#if defined(__MPI)
    !
    CALL mp_barrier(inter_pool_comm)
    CALL mp_sum(epsilon2_abs_dir, inter_pool_comm)
    CALL mp_sum(epsilon2_abs_lorenz_dir, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
#endif
    !
    ! Write to stdout
    !
    WRITE(stdout, '(5x,a)')
    WRITE(stdout, '(5x,a)') 'Direct absorption versus energy'
    WRITE(stdout, '(5x,a)') 'For the first temperature we have:'
    WRITE(stdout, '(5x,a)') 'Photon energy (eV), Imaginary dielectric function along x,y,z'
    DO iw = 1, nomega
      WRITE(stdout, '(5x,f15.6,3E22.14)') omegap(iw) * ryd2ev, (epsilon2_abs_dir(ipol, iw, 1), ipol = 1, 3)
    ENDDO
    WRITE(stdout, '(5x,a)')
    WRITE(stdout, '(5x,a)') 'Values with other for temperature X are reported in the files epsilon2_dirabs_X.dat'
    WRITE(stdout, '(5x,a)')
    !
    !Output to file
    !
    DO itemp = 1, nstemp
      IF (mpime == ionode_id) THEN
        WRITE(tp,"(f8.1)") gtemp(itemp) * ryd2ev / kelvin2eV
        nameF = 'epsilon2_dirabs_' // trim(adjustl(tp)) // 'K.dat'
        OPEN(UNIT = iudirabs, FILE = nameF)
        WRITE(iudirabs, '(a)') '# Direct absorption versus energy'
        WRITE(iudirabs, '(a)') '# Photon energy (eV), Imaginary dielectric function along x,y,z,average'
        DO iw = 1, nomega
          WRITE(iudirabs, '(5x,f15.6,4E22.14)') omegap(iw) * ryd2ev, (epsilon2_abs_dir(ipol, iw, itemp), ipol = 1, 3), &
                  SUM(epsilon2_abs_dir(:, iw, itemp)) / 3.0d0
        ENDDO
        CLOSE(iudirabs)
        !
        nameF = 'epsilon2_dirabs_lorenz' // trim(adjustl(tp)) // 'K.dat'
        OPEN(UNIT = iudirabs, FILE = nameF)
        WRITE(iudirabs, '(a)') '# Direct absorption versus energy'
        WRITE(iudirabs, '(a)') '# Photon energy (eV), Imaginary dielectric function along x,y,z,average'
        DO iw = 1, nomega
          WRITE(iudirabs, '(5x,f15.6,4E22.14)') omegap(iw) * ryd2ev, (epsilon2_abs_lorenz_dir(ipol, iw, itemp), ipol = 1, 3), &
                  SUM(epsilon2_abs_lorenz_dir(:, iw, itemp)) / 3.0d0
        ENDDO
        CLOSE(iudirabs)
      ENDIF
    ENDDO
    ! Deallocate
    DEALLOCATE(epsilon2_abs_dir, STAT = ierr)
    IF (ierr /= 0) CALL errore('dirabs', 'Error deallocating epsilon2_abs_dir', 1)
    DEALLOCATE(epsilon2_abs_lorenz_dir, STAT = ierr)
    IF (ierr /= 0) CALL errore('dirabs', 'Error deallocating epsilon2_abs_lorenz_dir', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE dirabs
    !-----------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE renorm_eig(ikk, ikq, nrr_k, dims, ndegen_k, irvec_k, irvec_r, cufkk, cufkq, cfac, cfacq)
    !--------------------------------------------------------------------------
    !!
    !! This routine computes the renormalization of the eigenenergies to be applied
    !! in case one read external eigenvalues.
    !! The implementation follows Eq. 30 of  Phys. Rev. B 62, 4927 (2000)
    !! Samuel Ponce, Kyle and Emmanouil Kioupakis
    !!
    USE kinds,         ONLY : DP
    USE elph2,         ONLY : xkfd, chw, chw_ks, etf_ks, etf, vmef, nkqf
    USE epwcom,        ONLY : use_ws, nbndsub
    USE constants_epw, ONLY : eps40, ryd2mev, twopi, zero, eps6, ci, czero
    USE wan2bloch,     ONLY : hamwan2bloch, vmewan2bloch
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ikk
    !! Current k point on that core (ikk = 2 * ik + 1)
    INTEGER, INTENT(in) :: ikq
    !! k+q point on that core
    INTEGER, INTENT(in) :: nrr_k
    !! Number of WS points for electrons
    INTEGER, INTENT(in) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(in) :: ndegen_k(nrr_k, dims, dims)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, INTENT(in) :: irvec_k(3, nrr_k)
    !! Integer components of the ir-th Wigner-Seitz grid point in the basis of the lattice vectors for electrons
    REAL(KIND = DP), INTENT(in) :: irvec_r(3, nrr_k)
    !!  Wigner-Size supercell vectors, store in real instead of integer
    COMPLEX(KIND = DP), INTENT(inout) :: cufkk(nbndsub, nbndsub)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP), INTENT(inout) :: cufkq(nbndsub, nbndsub)
    !! the same, for points k+q
    COMPLEX(KIND = DP)  :: cufkkd(nbndsub,nbndsub)
    !! Rotation matrix, shifted mesh, points k
    COMPLEX(KIND = DP)  :: cufkqd(nbndsub,nbndsub)
    !! Rotation matrix, shifted mesh, k+q
    COMPLEX(KIND = DP), INTENT(in) :: cfac(nrr_k, dims, dims)
    !! Exponential factor
    COMPLEX(KIND = DP), INTENT(in) :: cfacq(nrr_k, dims, dims)
    !! Exponential factor
    !
    ! Local variables
    INTEGER :: icounter
    !! Integer counter for displaced points
    INTEGER :: ibnd
    !! Band index
    INTEGER :: jbnd
    !! Band index
    INTEGER :: ir
    !! Counter for WS loop
    INTEGER :: iw
    !! Counter on bands when use_ws == .TRUE.
    INTEGER :: iw2
    !! Counter on bands when use_ws == .TRUE.
    REAL(KIND = DP) :: rdotk(nrr_k)
    !! $r\cdot k$
    REAL(KIND = DP) :: rdotk2(nrr_k)
    !! $r\cdot k$
    REAL(KIND = DP) :: etfd(nbndsub, nkqf, 6)
    !! interpolated eigenvalues (nbnd, nkqf) eigenvalues for shifted grid in the case of eig_read
    REAL(KIND = DP) :: etfd_ks(nbndsub, nkqf, 6)
    !! interpolated eigenvalues (nbnd, nkqf) KS eigenvalues for shifted grid in the case of eig_read
    COMPLEX(KIND = DP) :: cfacd(nrr_k, dims, dims, 6)
    !! Used to store $e^{2\pi r \cdot k}$ exponential of displaced vector
    COMPLEX(KIND = DP) :: cfacqd(nrr_k, dims, dims, 6)
    !! Used to store $e^{2\pi r \cdot k+q}$ exponential of dispaced vector
    !
    cfacd(:, :, :, :) = czero
    cfacqd(:, :, :, :)= czero
    DO icounter = 1, 6
      CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkfd(:, ikk, icounter), 1, 0.0_DP, rdotk, 1)
      CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkfd(:, ikq, icounter), 1, 0.0_DP, rdotk2, 1)
      IF (use_ws) THEN
        DO iw = 1, dims
          DO iw2 = 1, dims
            DO ir = 1, nrr_k
              IF (ndegen_k(ir, iw2, iw) > 0) THEN
                cfacd(ir, iw2, iw, icounter)  = EXP(ci * rdotk(ir))  / ndegen_k(ir, iw2, iw)
                cfacqd(ir, iw2, iw, icounter) = EXP(ci * rdotk2(ir)) / ndegen_k(ir, iw2, iw)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        cfacd(:, 1, 1, icounter)  = EXP(ci * rdotk(:)) / ndegen_k(:, 1, 1)
        cfacqd(:, 1, 1, icounter) = EXP(ci * rdotk2(:)) / ndegen_k(:, 1, 1)
      ENDIF
      !
      CALL hamwan2bloch(nbndsub, nrr_k, cufkkd, etfd(:, ikk, icounter), chw, cfacd, dims)
      CALL hamwan2bloch(nbndsub, nrr_k, cufkqd, etfd(:, ikq, icounter), chw, cfacqd, dims)
      CALL hamwan2bloch(nbndsub, nrr_k, cufkkd, etfd_ks(:, ikk, icounter), chw_ks, cfacd, dims)
      CALL hamwan2bloch(nbndsub, nrr_k, cufkqd, etfd_ks(:, ikq, icounter), chw_ks, cfacqd, dims)
    ENDDO ! icounter
    ! -----------------------------------------------------------------------------------------
    CALL vmewan2bloch(nbndsub, nrr_k, irvec_k, cufkk, vmef(:, :, :, ikk), &
                      etf(:, ikk), etf_ks(:, ikk), chw_ks, cfac, dims)
    CALL vmewan2bloch(nbndsub, nrr_k, irvec_k, cufkq, vmef(:, :, :, ikq), &
                      etf(:, ikq), etf_ks(:, ikq), chw_ks, cfacq, dims)
    ! To Satisfy Phys. Rev. B 62, 4927-4944 (2000) , Eq. (30)
    DO ibnd = 1, nbndsub
      DO jbnd = 1, nbndsub
        IF (ABS(etfd_ks(ibnd, ikk, 1) - etfd_ks(jbnd, ikk, 2)) > eps6) THEN
          vmef(1, ibnd, jbnd, ikk) = vmef(1, ibnd, jbnd, ikk) * (etfd(ibnd, ikk, 1) - etfd(jbnd, ikk, 2)) / &
              (etfd_ks(ibnd, ikk, 1) - etfd_ks(jbnd, ikk, 2))
        ENDIF
        IF (ABS(etfd_ks(ibnd, ikk, 3) - etfd_ks(jbnd, ikk, 4)) > eps6) THEN
          vmef(2, ibnd, jbnd, ikk) = vmef(2, ibnd, jbnd, ikk) * (etfd(ibnd, ikk, 3) - etfd(jbnd, ikk, 4)) / &
              (etfd_ks(ibnd, ikk, 3) - etfd_ks(jbnd, ikk, 4))
        ENDIF
        IF (ABS(etfd_ks(ibnd, ikk, 5) - etfd_ks(jbnd, ikk, 6)) > eps6) THEN
          vmef(3, ibnd, jbnd, ikk) = vmef(3, ibnd, jbnd, ikk) * (etfd(ibnd, ikk, 5) - etfd(jbnd, ikk, 6)) / &
              (etfd_ks(ibnd, ikk, 5) - etfd_ks(jbnd, ikk, 6))
        ENDIF
        IF (ABS(etfd_ks(ibnd, ikq, 1) - etfd_ks(jbnd, ikq, 2)) > eps6) THEN
          vmef(1, ibnd, jbnd, ikq) = vmef(1, ibnd, jbnd, ikq) * (etfd(ibnd, ikq, 1) - etfd(jbnd, ikq, 2)) / &
              (etfd_ks(ibnd, ikq, 1) - etfd_ks(jbnd, ikq, 2))
        ENDIF
        IF (ABS(etfd_ks(ibnd, ikq, 3) - etfd_ks(jbnd, ikq, 4)) > eps6) THEN
          vmef(2, ibnd, jbnd, ikq) = vmef(2, ibnd, jbnd, ikq) * (etfd(ibnd, ikq, 3) - etfd(jbnd, ikq, 4)) / &
              (etfd_ks(ibnd, ikq, 3) - etfd_ks(jbnd, ikq, 4))
        ENDIF
        IF (ABS(etfd_ks(ibnd, ikq, 5) - etfd_ks(jbnd, ikq, 6)) > eps6) THEN
          vmef(3, ibnd, jbnd, ikq) = vmef(3, ibnd, jbnd, ikq) * (etfd(ibnd, ikq, 5) - etfd(jbnd, ikq, 6) ) / &
              (etfd_ks(ibnd, ikq, 5) - etfd_ks(jbnd, ikq, 6))
        ENDIF
      ENDDO
    ENDDO
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE renorm_eig
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE fermi_carrier_indabs(itemp, etemp_fca, ef0_fca, ctype)
    !-----------------------------------------------------------------------
    !! Xiao Zhang: Implemented 03/2021
    !! This is a slightly modified version of subroutine fermicarrier
    !! in Utility.f90. This is used to calculate the fermi energy associated
    !! with a certain carrier density in when the user want to calculate
    !! indirect optical absorption for doped semiconductor
    !!
    !-----------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : omega, alat, at
    USE io_global, ONLY : stdout
    USE elph2,     ONLY : etf, nkf, wkf, efnew, nkqf, evbm, ecbm
    USE constants_epw, ONLY : ryd2ev, bohr2ang, ang2cm, eps5, kelvin2eV, &
                              zero, eps80, eps6
    USE noncollin_module, ONLY : noncolin
    USE pwcom,     ONLY : nelec
    USE epwcom,    ONLY : int_mob, nbndsub, ncarrier, nstemp, fermi_energy, &
                          system_2d, carrier, efermi_read, assume_metal, ngaussw
    USE klist_epw, ONLY : isk_dummy
    USE mp,        ONLY : mp_barrier, mp_sum, mp_max, mp_min
    USE mp_global, ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Temperature index
    INTEGER, INTENT(out) :: ctype
    !! Calculation type: -1 = hole, +1 = electron and 0 = both.
    REAL(KIND = DP), INTENT(in) :: etemp_fca
    !! Temperature in kBT [Ry] unit.
    REAL(KIND = DP), INTENT(inout) :: ef0_fca(nstemp)
    !! Fermi level for the temperature itemp
!    REAL(KIND = DP), INTENT(inout) :: efcb(nstemp)
    !! Second fermi level for the temperature itemp
    REAL(KIND = DP), EXTERNAL :: efermig
    !! External function to calculate the fermi energy
    !
    ! Local variables
    INTEGER :: i
    !! Index for the bisection iteration
    INTEGER :: ik
    !! k-point index per pool
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: ivbm
    !! Index of the VBM
    INTEGER :: icbm
    !! Index of the CBM
    INTEGER, PARAMETER :: maxiter = 500 ! 300
    !! Maximum interation
    REAL(KIND = DP) :: fermi
    !! Fermi level returned
    REAL(KIND = DP) :: fermicb
    !! Fermi level returned for second Fermi level
    REAL(KIND = DP) :: fnk
    !! Fermi-Diract occupation
    REAL(KIND = DP) :: ks_exp(nbndsub, nkf)
    !! Exponential of the eigenvalues divided by kBT
    REAL(KIND = DP) :: ks_expcb(nbndsub, nkf)
    !! Exponential of the eigenvalues divided by kBT for CB
    REAL(KIND = DP) :: fermi_exp
    !! Fermi level in exponential format
    REAL(KIND = DP) :: rel_err
    !! Relative error
    REAL(KIND = DP) :: factor
    !! Factor that goes from number of carrier per unit cell to number of
    !! carrier per cm^-3
    REAL(KIND = DP) :: arg
    !! Argument of the exponential
    REAL(KIND = DP) :: inv_cell
    !! Inverse of the volume in [Bohr^{-3}]
    REAL(KIND = DP) :: ef_tmp
    !! Energy of the current Fermi level for the bisection method
    REAL(KIND = DP) :: elw
    !! Energy lower bound for the bisection method
    REAL(KIND = DP) :: eup
    !! Energy upper bound for the bisection method
    REAL(KIND = DP) :: hole_density
    !! Hole carrier density
    REAL(KIND = DP) :: electron_density
    !! Electron carrier density
    REAL(KIND = DP) :: intrinsic_density
    !! Intrinsic carrier density
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP) :: ekk
    !! Eigen energy on the fine grid relative to the Fermi level
    !
    IF (assume_metal) THEN
      !! set conduction band chemical potential to 0 since it is irrelevent
 !     efcb(itemp) = zero
      ef0_fca(itemp) = efermig(etf, nbndsub, nkqf, nelec, wkf, etemp_fca, ngaussw, 0, isk_dummy)
      RETURN
    ENDIF
    !
    IF (ncarrier < 0.d0) THEN
      ctype = -1
    ELSE
      ctype = 1
    ENDIF
    !
    ef_tmp  = zero
    fermi   = zero
    fermicb = zero
    !
    IF (system_2d == 'no') THEN
      inv_cell = 1.0d0 / omega
    ELSE 
      ! for 2d system need to divide by area (vacuum in z-direction)
      inv_cell = ( 1.0d0 / omega ) * at(3, 3) * alat
    ENDIF
    ! vbm index
    IF (noncolin) THEN
      ivbm = FLOOR(nelec / 1.0d0)
    ELSE
      ivbm = FLOOR(nelec / 2.0d0)
    ENDIF
    icbm = ivbm + 1 ! Nb of bands
    !
    ! Initialization value. Should be large enough ...
    evbm = -10000d0
    ecbm = 10000d0 ! In Ry
    !
    WRITE(stdout, '(5x, "calculating band extrema...")')
    DO ik = 1, nkf
      ikk = 2 * ik - 1
      DO ibnd = 1, nbndsub
        IF (ibnd < ivbm + 1) THEN
          IF (etf(ibnd, ikk) > evbm) THEN
            evbm = etf(ibnd, ikk)
          ENDIF
        ENDIF
       ! Find cbm index
        IF (ibnd > ivbm) THEN
          IF (etf(ibnd, ikk) < ecbm) THEN
            ecbm = etf(ibnd, ikk)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !
    ! Find max and min across pools
    !
    CALL mp_max(evbm, inter_pool_comm)
    CALL mp_min(ecbm, inter_pool_comm)
    !
    IF (itemp == 1) THEN
      WRITE(stdout, '(5x, "Valence band maximum    = ", f10.6, " eV")') evbm * ryd2ev
      WRITE(stdout, '(5x, "Conduction band minimum = ", f10.6, " eV")') ecbm * ryd2ev
    ENDIF
    !
    ! We first calculate intrinsic carrier density
    ! Using this we can determine if the user input make sense
    WRITE(stdout, '(5x, "Calculating intrinsic carrier density")')
    ef_tmp = (ecbm + evbm) / 2.d0
    eup = ecbm
    elw = evbm
    factor = inv_cell * (bohr2ang * ang2cm)**(-3.d0)
    Do i = 1, maxiter
      !
      ef_tmp = (eup + elw) / 2.d0
      hole_density = zero
      electron_density = zero
      !
      DO ik = 1, nkf
        ikk = 2 * ik -1
!        WRITE(stdout, '(5x, i, i)') icbm, nbndsub
        DO ibnd = icbm, nbndsub
          ekk = etf(ibnd, ikk) - ef_tmp
          fnk = wgauss(-ekk / etemp_fca, -99)
          ! The wkf(ikk) already include a factor 2
          electron_density = electron_density + wkf(ikk) * fnk * factor
        ENDDO ! ibnd
      ENDDO ! ik
      DO ik = 1, nkf
        ikk = 2 * ik -1
        DO ibnd = 1, ivbm
          ekk = etf(ibnd, ikk) - ef_tmp
          fnk = wgauss(-ekk / etemp_fca, -99)
          ! The wkf(ikk) already include a factor 2
          hole_density = hole_density + wkf(ikk) * (1.0d0 - fnk) * factor
        ENDDO ! ibnd
      ENDDO ! ik
      CALL mp_barrier(inter_pool_comm)
      CALL mp_sum(hole_density, inter_pool_comm)
      CALL mp_sum(electron_density, inter_pool_comm)
      IF (ABS(hole_density) < eps80) THEN
        rel_err = -1.d0
      ELSE
        rel_err = (hole_density - electron_density) / hole_density
      ENDIF
      IF (ABS(rel_err) < eps5) THEN
        intrinsic_density = hole_density
        EXIT
      ELSEIF ((rel_err) > eps5) THEN
        elw = ef_tmp
      ELSE
        eup = ef_tmp
      ENDIF
      intrinsic_density = hole_density
    ENDDO ! maxiter
    IF (i == maxiter) THEN
      WRITE(stdout, '(5x, a)') "Too many iterations when calculating intrinsic density"
      WRITE(stdout, '(5x, a, f8.5)') "Relative error of hole density versus electron density: ", rel_err
      WRITE(stdout, '(5x, a)') "Something likely wrong unless we are dealing with a metallic system"
    ENDIF
    !
    WRITE(stdout, '(5x, "Intrinsic density = ", E18.6, "Cm^-3")' ) intrinsic_density
    WRITE(stdout, '(/5x, "calculating fermi level...")')
    IF (ABS(ncarrier) < intrinsic_density) THEN
      WRITE(stdout, '(5x, a)') 'ncarrier not given, or smaller than intrinsic density'
      WRITE(stdout, '(5x, a)') 'Setting the fermi level to mig-gap'
      fermi = (evbm + ecbm) / 2.d0
    ELSEIF ( ncarrier > intrinsic_density ) THEN
      ! assuming free electron density
      eup = 10000d0 + ecbm
      elw = evbm - 10000d0
      ef_tmp = (ecbm + evbm) / 2.d0
      DO i = 1, maxiter
        electron_density = zero
        DO ik = 1, nkf
          ikk = 2 * ik -1
          DO ibnd = icbm, nbndsub
            ekk = etf(ibnd, ikk) - ef_tmp
            fnk = wgauss(-ekk / etemp_fca, -99)
            ! The wkf(ikk) already include a factor 2
            electron_density = electron_density + wkf(ikk) * fnk * factor
          ENDDO ! ibnd
        ENDDO ! ik
        CALL mp_barrier(inter_pool_comm)
        CALL mp_sum(electron_density, inter_pool_comm)
        ! In this case ncarrier is a negative number
        rel_err = (electron_density - ncarrier) / electron_density
        !
        IF (ABS(rel_err) < eps6) THEN
          fermi = ef_tmp
          EXIT
        ELSEIF ((rel_err) > eps6) THEN
          eup = ef_tmp
        ELSE
          elw = ef_tmp
        ENDIF
        IF (MOD(i, 10) == 0) THEN
          WRITE(stdout, '(5x, a, i5, a, f8.3)') 'Iteration #', i, " Fermi level = ", ef_tmp
        ENDIF
        ef_tmp = (eup + elw) / 2.0d0
      ENDDO ! maxiter
      fermi = ef_tmp
    ELSEIF ( ncarrier < - intrinsic_density) THEN
      ! assuming free hole density
      eup = 10000d0 + ecbm
      elw = evbm - 10000d0
      ef_tmp = (ecbm + evbm) / 2.d0
      DO i = 1, maxiter
        hole_density = zero
        DO ik = 1, nkf
          ikk = 2 * ik -1
          DO ibnd = 1, ivbm
            ekk = etf(ibnd, ikk) - ef_tmp
            fnk = wgauss(-ekk / etemp_fca, -99)
            ! The wkf(ikk) already include a factor 2
            hole_density = hole_density + wkf(ikk) * (1.0d0 - fnk) * factor
          ENDDO ! ibnd
        ENDDO ! ik
        CALL mp_barrier(inter_pool_comm)
        CALL mp_sum(hole_density, inter_pool_comm)
        rel_err = (hole_density - ABS(ncarrier)) / hole_density
        !
        IF (ABS(rel_err) < eps6) THEN
          fermi = ef_tmp
          EXIT
        ELSEIF ((rel_err) > eps6) THEN
          elw = ef_tmp
        ELSE
          eup = ef_tmp
        ENDIF
        IF (MOD(i, 10) == 0) THEN
          WRITE(stdout, '(5x, a, i5, a, f8.3)') 'Iteration #', i, " Fermi level = ", ef_tmp
        ENDIF
        ef_tmp = (eup + elw) / 2.0d0
      ENDDO ! maxiter
      fermi = ef_tmp
    ENDIF !ncarrier
    IF (i == maxiter) THEN
      WRITE(stdout, '(5x, "Warning: too many iterations in bisection"/ &
                    5x, "ef_tmp = ", f10.6)' ) fermi * ryd2ev
    ENDIF
    WRITE(stdout, '(/5x, "Temperature ", f8.3, " K")' ) etemp_fca * ryd2ev / kelvin2eV
    IF (ncarrier > intrinsic_density) THEN
      ef0_fca(itemp) = fermi
      WRITE(stdout, '(5x, "Electron density = ", E18.6, "Cm^-3")' ) electron_density
      WRITE(stdout, '(5x, "Calculated Fermi level = ", f10.5, " eV")' )  ef0_fca(itemp) * ryd2ev
    ELSEIF (ncarrier < - intrinsic_density) THEN
      ef0_fca(itemp) = fermi
      WRITE(stdout, '(5x, "Hole density = ", E18.6, "Cm^-3")' ) hole_density
      WRITE(stdout, '(5x, "Calculated Fermi level = ", f10.5, " eV")' )  ef0_fca(itemp) * ryd2ev
    ELSE
      ef0_fca(itemp) = fermi
      WRITE(stdout, '(5x, "Calculated Fermi level = ", f10.5, " eV")' )  ef0_fca(itemp) * ryd2ev
    ENDIF
    !
    RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE fermi_carrier_indabs
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE conduc_fca(ef0_fca)
    !----------------------------------------------------------------------
    !!
    !! This routine calculates the electrical conductivity
    !! within constant relaxation time approxiamtion
    !! The aim is to calculate resistive term within EPW in a
    !! consistent manner
    !!
    !! Xiao Zhang, 06/2021 First implementation
    !----------------------------------------------------------------------
    USE kinds,             ONLY : DP
    USE io_global,         ONLY : stdout, ionode_id
    USE cell_base,         ONLY : alat, at, omega, bg
    USE symm_base,         ONLY : s
    USE epwcom,            ONLY : nbndsub, fsthick, system_2d, nstemp, assume_metal, &
                                  vme, mp_mesh_k, nkf1, nkf2, nkf3, omegamin, omegamax, &
                                  omegastep, nomega, sigma_ref
    USE elph2,             ONLY : ibndmin, etf, nkf, wkf, vmef, dmef, bztoibz,  &
                                  nkqtotf, gtemp, nbndfst, nktotf, nkqf, s_bztoibz, &
                                  omegap
    USE constants_epw,     ONLY : zero, one, bohr2ang, ryd2ev, ang2cm, czero, &
                                  kelvin2eV, hbar, Ang2m, hbarJ, eps6, eps4, pi, &
                                  ryd2mev, meV2invps
    USE constants,         ONLY : electron_si
    USE mp,                ONLY : mp_sum, mp_bcast
    USE mp_global,         ONLY : world_comm
    USE mp_world,          ONLY : mpime
    USE poolgathering,     ONLY : poolgatherc4, poolgather2
    USE division,          ONLY : fkbounds
    USE grid,              ONLY : kpoint_grid_epw
    USE symm_base,         ONLY : s
    USE noncollin_module,  ONLY : noncolin
    USE pwcom,             ONLY : ef
    USE io_var,            ONLY : iuindabs
    !
    IMPLICIT NONE
    !
    REAL (KIND = DP), INTENT(in) :: ef0_fca(nstemp)
    !! Fermi level for temperature itemp and the given carrier density.
    !
    ! Local variables
    CHARACTER(LEN = 256) :: format_string
    !! Format string
    CHARACTER(LEN = 256) :: nameF
    !! Name of the file
    CHARACTER(LEN = 20) :: tp
    !! Temperature, in string format
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
    INTEGER :: iww
    !! Frequency point
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: efcalc
    !! Fermi level
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
    REAL(KIND = DP) :: sigma_ref_au
    !! Reference conductivity from user input, a.u.
    REAL(KIND = DP) :: sigma_calc
    !! Calculated average conductivity
    REAL(KIND = DP) :: inv_cell
    !! Inverse of the volume in [Bohr^{-3}]
    REAL(KIND = DP) :: vkk(3, nbndfst)
    !! Electron velocity vector for a band.
    REAL(KIND = DP) :: sigma(9, nstemp)
    !! Conductivity matrix in vector form
    REAL(KIND = DP) :: sigma_m(3, 3)
    !! Conductivity matrix
    REAL(KIND = DP) :: sigma_eig(3)
    !! Eigenvalues from the diagonalized conductivity matrix
    REAL(KIND = DP) :: sigma_vect(3, 3)
    !! Eigenvectors from the diagonalized conductivity matrix
    REAL(KIND = DP) :: tdf_sigma(9)
    !! Temporary file
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    REAL(KIND = DP) :: epsilon2_resistive(nomega)
    !! Resistive spectra
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
    IF (system_2d == 'no') THEN
      inv_cell = 1.0d0 / omega
    ELSE 
      ! for 2d system need to divide by area (vacuum in z-direction)
      inv_cell = ( 1.0d0 / omega ) * at(3, 3) * alat
    ENDIF
    !
    conv_factor1 = electron_si / (hbar * bohr2ang * Ang2m)
    !
    CALL fkbounds(nktotf, lower_bnd, upper_bnd)
    !
    DO itemp = 1, nstemp
      !
      etemp = gtemp(itemp)
      efcalc = ef0_fca(itemp)
      !
      IF (itemp == 1) THEN
        WRITE(stdout, '(/5x, a)') 'Calculate conducitivity within RTA for the given density'
        tdf_sigma(:) = zero
        sigma(:, :)  = zero
        !
      ENDIF
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        !
        IF (MINVAL(ABS(etf(:, ik) - ef)) < fsthick) THEN
          !
          DO ibnd = 1, nbndfst
            !
            tdf_sigma(:) = zero
            !
            IF (vme == 'wannier') THEN
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
            ekk = etf(ibndmin - 1 + ibnd, ikk) - efcalc
            !
            ! derivative Fermi distribution
            dfnk = w0gauss(ekk / etemp, -99) / etemp
            !
            ! electrical conductivity matrix
            sigma(:, itemp) = sigma(:, itemp) + dfnk * tdf_sigma(:)! * tau
            !
          ENDDO! ibnd
        ENDIF!fsthick
      ENDDO! ik
      !
      !Sum over all pool to gather results
      !
      CALL mp_sum(sigma(:, itemp), world_comm)
      !
      sigma_m(:, :) = zero
      sigma_m(1, 1) = sigma(1, itemp)
      sigma_m(1, 2) = sigma(2, itemp)
      sigma_m(1, 3) = sigma(3, itemp)
      sigma_m(2, 1) = sigma(4, itemp)
      sigma_m(2, 2) = sigma(5, itemp)
      sigma_m(2, 3) = sigma(6, itemp)
      sigma_m(3, 1) = sigma(7, itemp)
      sigma_m(3, 2) = sigma(8, itemp)
      sigma_m(3, 3) = sigma(9, itemp)
      ! Diagonalize the conductivity matrix
      CALL rdiagh(3, sigma_m(:, :), 3, sigma_eig, sigma_vect)
      !
      !
      sigma_calc = SUM(sigma_eig) / 3.d0 * inv_cell
      sigma_ref_au = sigma_ref / conv_factor1
      tau = sigma_ref_au / sigma_calc
      WRITE(stdout, '(5x, a, 3E16.7)') 'Calculated constant relaxation time: ', tau / (ryd2mev * meV2invps)
      WRITE(stdout, '(5x, a, 3E16.7)') 'Conductivity xx, yy, zz: ', conv_factor1 * sigma_eig(:) * &
                                        inv_cell * tau
      WRITE(stdout, '(5x, a)') 'Calculate and Write the resistive contribution (Drude term).'
      !
      !4*pi*sigma/(w*(1+w^2*tau^2)), an additional factor of two comes from e^2
      !
      DO iww = 1, nomega
        epsilon2_resistive(iww) = 8.d0 * pi * sigma_ref_au / (omegap(iww) * (1 + omegap(iww)**2.d0 * tau**2.d0))
      ENDDO
      !
      IF (mpime == ionode_id) THEN
        WRITE(tp,"(f8.1)") gtemp(itemp) * ryd2ev / kelvin2eV
        nameF = 'epsilon2_indabs_resis_' // trim(adjustl(tp)) // 'K.dat'
        OPEN(UNIT = iuindabs, FILE = nameF)
        WRITE(iuindabs, '(a)') '# Resistive contribution versus energy'
        WRITE(iuindabs, '(a)') '# Photon energy (eV), Directionally-averaged imaginary dielectric function along x,y,z'
        DO iww = 1, nomega
          WRITE(iuindabs, '(2E22.14)') omegap(iww) * ryd2ev, epsilon2_resistive(iww)
        ENDDO
        CLOSE(iuindabs)
      ENDIF
      !
    ENDDO ! itemp
    !-----------------------------------------------------------------------
    END SUBROUTINE conduc_fca
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE prepare_indabs()
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iuindabs
    USE modes,         ONLY : nmodes
    USE epwcom,        ONLY : nstemp, omegamin, omegamax, omegastep,     &
                              nomega, neta, restart

    USE elph2,         ONLY : epsilon2_abs, epsilon2_abs_lorenz,         &
                              epsilon2_abs_all, epsilon2_abs_lorenz_all
    !
    IMPLICIT NONE
    !
    INTEGER  :: ierr
    !! error in allocation
    !                    
    ! Calculate the number of frequency points
    nomega = INT((omegamax - omegamin) / omegastep) + 1
    ALLOCATE(epsilon2_abs(3, nomega, neta, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('prepare_indabs', 'Error allocating epsilon2_abs', 1)
    ALLOCATE(epsilon2_abs_lorenz(3, nomega, neta, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('prepare_indabs', 'Error allocating epsilon2_abs_lorenz', 1)
    ALLOCATE(epsilon2_abs_all(3, nomega, neta, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('prepare_indabs', 'Error allocating epsilon2_abs_all', 1)
    ALLOCATE(epsilon2_abs_lorenz_all(3, nomega, neta, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('prepare_indabs', 'Error allocating & 
                              &epsilon2_abs_lorenz_all', 1)
    !-------------------------------------------------------------------------
    END SUBROUTINE prepare_indabs 
    !-------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  END MODULE indabs
  !-----------------------------------------------------------------------


